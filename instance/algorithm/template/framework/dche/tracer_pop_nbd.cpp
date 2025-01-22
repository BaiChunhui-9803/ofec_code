#include "tracer_pop_nbd.h"
#include "../../../../../core/global.h"

namespace ofec {
	void TracerPopNBD::addInputParameters() {
		m_input_parameters.add("model file", new FileName(m_model_file_name,
			"instance/algorithm/template/framework/edhe",
			"Models (*.pt)", "nea2_40_50.pt"));
	}

	void TracerPopNBD::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		try {
			m_img_seg_mdl = torch::jit::load(g_working_directory +
				"/instance/algorithm/template/framework/edhe/" + m_model_file_name);
#ifdef OFEC_PLAYBACK
			std::cout << "The model has been successfully loaded." << std::endl;
#endif // OFEC_PLAYBACK
		}
		catch (const c10::Error &e) {
			std::cout << "Failed to load the model." << std::endl;
			std::cout << e.what() << std::endl;
			throw Exception("Error at EDHE::initialize_(): failed to load the model.");
		}
		size_t last_delimiter_pos = m_model_file_name.rfind('_');
		size_t second_last_delimieter_pso = m_model_file_name.rfind('_', last_delimiter_pos - 1);
		m_nbd_resolution = std::stoi(m_model_file_name.substr(last_delimiter_pos + 1));
		//m_max_num_iters = std::stoi(m_model_file_name.substr(second_last_delimieter_pso + 1,
			//last_delimiter_pos - second_last_delimieter_pso));
		m_max_num_iters = env->problem()->numberVariables() * 20;
		clearEvoPop();
	}

	void TracerPopNBD::clearEvoPop() {
		m_evo_pop.clear();
		m_evo_nbd.clear();
		m_outlier.clear();
		m_evo_nbd2.clear();
		m_sols.clear();
	}

	void TracerPopNBD::identifyOutliers(std::list<const Solution<>*> &candidates, Environment *env) {
		int num_iters = m_evo_pop.size();
		Real max_nbd = 0;
		for (size_t iter = 0; iter < num_iters; ++iter) {
			Real max_nbd_cur_iter = *std::max_element(m_evo_nbd[iter].begin(), m_evo_nbd[iter].end());
			if (max_nbd_cur_iter > max_nbd)
				max_nbd = max_nbd_cur_iter;
		}
		max_nbd *= 1.1;
		//Real max_freq = 0;
		std::vector<std::vector<Real>> freq_heat_map_data(num_iters, std::vector<Real>(m_nbd_resolution, 0));
		for (size_t iter = 0; iter < num_iters; ++iter) {
			for (Real nbd : m_evo_nbd[iter]) {
				freq_heat_map_data[iter][m_nbd_resolution * nbd / max_nbd]++;
				//if (freq > max_freq)
				//	max_freq = freq;
			}
		}

		auto input_tensor = torch::empty({ 1, 1, m_nbd_resolution, num_iters });
		for (int x = 0; x < m_nbd_resolution; ++x) {
			for (int y = 0; y < num_iters; ++y) {
				int gs = 2 / acos(-1.0) * atan(freq_heat_map_data[y][m_nbd_resolution - 1 - x]) * 255;
				input_tensor[0][0][x][y] = gs;
			}
		}
		std::vector<torch::jit::IValue> inputs({ input_tensor });
		auto output_tensor = m_img_seg_mdl.forward(inputs).toTensor();

		std::vector<std::vector<bool>> is_outlier(num_iters, std::vector<bool>(m_nbd_resolution));
		for (int y = 0; y < num_iters; ++y) {
			for (int x = 0; x < m_nbd_resolution; ++x) {
				if (output_tensor[0][0][x][y].item<float>() > 0.5)
					is_outlier[y][m_nbd_resolution - 1 - x] = true;
				else
					is_outlier[y][m_nbd_resolution - 1 - x] = false;
			}
		}

		for (size_t i = 0; i < m_evo_nbd.size(); ++i) {
			for (size_t j = 0; j < m_evo_nbd[i].size(); ++j) {
				int y = i;
				int x = m_nbd_resolution * m_evo_nbd[i][j] / max_nbd;
				m_outlier[i][j] = is_outlier[y][x];
			}
		}

		std::vector<std::vector<const SolutionBase*>> seeds_each_iter(num_iters);
		for (size_t iter = 0; iter < num_iters; ++iter) {
			for (size_t i = 0; i < m_evo_nbd2[iter].size(); ++i) {
				if (m_evo_nbd2[iter][i] == -1) {
					seeds_each_iter[iter].push_back(m_evo_pop[iter][i]);
				}
				else {
					int x = m_nbd_resolution * m_evo_nbd2[iter][i] / max_nbd;
					if (is_outlier[iter][x]) {
						seeds_each_iter[iter].push_back(m_evo_pop[iter][i]);
					}
				}
			}
		}
		std::vector<const SolutionBase*> solutions;
		for (size_t iter = 0; iter < num_iters; ++iter) {
			if (solutions.size() >= seeds_each_iter[iter].size()) {
				for (auto seed : seeds_each_iter[iter]) {
					Real min_dis = std::numeric_limits<Real>::max();
					size_t nearest = solutions.size();
					for (size_t i = 0; i < solutions.size(); ++i) {
						Real dis = seed->variableDistance(*solutions[i], env);
						if (dis < min_dis) {
							min_dis = dis;
							nearest = i;
						}
					}
					if (dominate(*seed, *solutions[nearest], env->problem()->optimizeMode())) {
						solutions[nearest] = seed;
					}
				}
			}
			else {
				std::vector<bool> mapped(seeds_each_iter[iter].size(), false);
				for (auto &seed : solutions) {
					Real min_dis = std::numeric_limits<Real>::max();
					size_t nearest = seeds_each_iter[iter].size();
					for (size_t i = 0; i < seeds_each_iter[iter].size(); ++i) {
						Real dis = seed->variableDistance(*seeds_each_iter[iter][i], env);
						if (dis < min_dis) {
							min_dis = dis;
							nearest = i;
						}
					}
					mapped[nearest] = true;
					if (dominate(*seeds_each_iter[iter][nearest], *seed, env->problem()->optimizeMode()))
						seed = seeds_each_iter[iter][nearest];
				}
				for (size_t i = 0; i < seeds_each_iter[iter].size(); ++i) {
					if (!mapped[i])
						solutions.push_back(seeds_each_iter[iter][i]);
				}
			}
		}
		for (auto s : solutions) {
			candidates.push_back(dynamic_cast<const Solution<>*>(s));
		}
	}
}
