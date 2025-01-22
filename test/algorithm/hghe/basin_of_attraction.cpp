/********* Begin Register Information **********
{ "description": "Estimation of basins of attraction based on space-paritioning" }
*********** End Register Information **********/

#include <fstream>
#include <sstream>
#include <iomanip>
#include "../../../utility/catch.hpp"
#include "../../../run/include_problem.h"
#include "../../../core/global.h"
#include "../../../core/problem/solution.h"
#include "../../../utility/clustering/nbc.h"
#include "../../../utility/kd-tree/kdtree_space.h"
#include "../../../utility/nondominated_sorting/filter_sort.h"
#include "../../../instance/algorithm/continuous/single_objective/multi_modal/de_nrand_1/de_nrand_1_pop.h"
#include "../../../utility/kd-tree/data_adaptor.h"

using namespace ofec;

void monotonicSearching(
	const std::vector<Real> &obj,
	const std::vector<std::list<size_t>> &neighbor,
	size_t center,
	std::set<size_t> &cluster,
	std::vector<size_t> &num_clustered,
	bool is_minimize)
{
	for (auto i : neighbor[center]) {
		if ((is_minimize && obj[i] >= obj[center])
			|| (!is_minimize && obj[i] <= obj[center])) {
			if (cluster.count(i) == 0) {
				cluster.emplace(i);
				num_clustered[i]++;
				monotonicSearching(obj, neighbor, i, cluster, num_clustered, is_minimize);
			}
		}
	}
}

void groupSubspace(
	const std::vector<Real> &obj,
	const std::vector<std::list<size_t>> &neighbor,
	std::vector<std::list<size_t>> &aff_bsn,
	bool is_minimize,
	size_t &num_bsns)
{
	std::multimap<Real, size_t> sort_fitness;
	std::vector<size_t> num_clustered(obj.size(), 0);
	std::list<std::set<size_t>> clusters;
	for (size_t i = 0; i < obj.size(); ++i)
		sort_fitness.emplace(obj[i], i);
	if (is_minimize) {
		for (auto it = sort_fitness.begin(); it != sort_fitness.end(); ++it) {
			if (num_clustered[it->second] == 0) {
				std::set<size_t> cluster;
				cluster.emplace(it->second);
				num_clustered[it->second]++;
				monotonicSearching(obj, neighbor, it->second, cluster, num_clustered, is_minimize);
				clusters.emplace_back(cluster);
			}
		}
	}
	else {
		for (auto it = sort_fitness.rbegin(); it != sort_fitness.rend(); ++it) {
			if (num_clustered[it->second] == 0) {
				std::set<size_t> cluster;
				cluster.emplace(it->second);
				num_clustered[it->second]++;
				monotonicSearching(obj, neighbor, it->second, cluster, num_clustered, is_minimize);
				clusters.emplace_back(cluster);
			}
		}
	}
	num_bsns = 0;
	aff_bsn.clear();
	aff_bsn.resize(obj.size());
	for (auto &cluster : clusters) {
		for (auto i : cluster)
			aff_bsn[i].push_back(num_bsns);
		num_bsns++;
	}
}

void generateInfoBoA(const std::string &pro_name, int num_vars, size_t num_div) {
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;

	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);
	bool is_minimize = pro->optimizeMode(0) == OptimizeMode::kMinimize;
	size_t num_ssp, id_ssp, num_bsns;
	Real gap_x1, gap_x2, l_x1, l_x2;
	num_ssp = num_div * num_div;
	auto &range_x1 = CAST_CONOP(pro)->range(0);
	auto &range_x2 = CAST_CONOP(pro)->range(1);
	gap_x1 = (range_x1.second - range_x1.first) / num_div;
	gap_x2 = (range_x2.second - range_x2.first) / num_div;

	std::vector<Real> obj(num_ssp, (is_minimize ? std::numeric_limits<Real>::max() : -std::numeric_limits<Real>::max()));
	Solution<> s(1, 0, num_vars);
	for (size_t i = 0; i < num_div; ++i) {
		l_x1 = range_x1.first + i * gap_x1;
		for (size_t j = 0; j < num_div; ++j) {
			l_x2 = range_x2.first + j * gap_x2;
			id_ssp = i * num_div + j;
			for (size_t k = 0; k < 100; ++k) {
				s.variable()[0] = rnd->uniform.nextNonStd(l_x1, l_x1 + gap_x1);
				s.variable()[1] = rnd->uniform.nextNonStd(l_x2, l_x2 + gap_x2);
				s.evaluate(pro, -1, false);
				if (is_minimize) {
					if (s.objective(0) < obj[id_ssp])
						obj[id_ssp] = s.objective(0);
				}
				else if (s.objective(0) > obj[id_ssp])
					obj[id_ssp] = s.objective(0);;
			}
		}
	}
	std::vector<std::list<size_t>> aff_bsn;
	if (num_div < 3) {
		aff_bsn.assign(num_ssp, { 0 });
	}
	else {
		std::vector<std::list<size_t>> neighbor(num_ssp);
		for (size_t i = 0; i < 1; ++i) {
			for (size_t j = 0; j < 1; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp + 1, id_ssp + num_div, id_ssp + num_div + 1 };
			}
			for (size_t j = 1; j < num_div - 1; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp - 1, id_ssp + 1, id_ssp + num_div - 1,id_ssp + num_div, id_ssp + num_div + 1 };
			}
			for (size_t j = num_div - 1; j < num_div; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp - 1, id_ssp + num_div - 1, id_ssp + num_div };
			}
		}
		for (size_t i = 1; i < num_div - 1; ++i) {
			for (size_t j = 0; j < 1; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp - num_div, id_ssp - num_div + 1, id_ssp + 1, id_ssp + num_div, id_ssp + num_div + 1 };
			}
			for (size_t j = 1; j < num_div - 1; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp - num_div - 1, id_ssp - num_div, id_ssp - num_div + 1, id_ssp - 1, id_ssp + 1, id_ssp + num_div - 1,id_ssp + num_div, id_ssp + num_div + 1 };
			}
			for (size_t j = num_div - 1; j < num_div; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp - num_div - 1, id_ssp - num_div,  id_ssp - 1, id_ssp + num_div - 1, id_ssp + num_div };
			}
		}
		for (size_t i = num_div - 1; i < num_div; ++i) {
			for (size_t j = 0; j < 1; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp - num_div, id_ssp - num_div + 1, id_ssp + 1 };
			}
			for (size_t j = 1; j < num_div - 1; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp - num_div - 1, id_ssp - num_div, id_ssp - num_div + 1, id_ssp - 1, id_ssp + 1 };
			}
			for (size_t j = num_div - 1; j < num_div; ++j) {
				id_ssp = i * num_div + j;
				neighbor[id_ssp] = { id_ssp - num_div - 1, id_ssp - num_div,  id_ssp - 1 };
			}
		}
		groupSubspace(obj, neighbor, aff_bsn, is_minimize, num_bsns);
	}

	std::stringstream out_sstream;
	out_sstream << "# This file stores the information of basin of attractions of problem\n";
	out_sstream << "# problem name: " << pro_name << "\n"
		<< "# number of decision variables: " << num_vars << "\n"
		<< "# number of divisions for each variable: " << num_div << "\n"
		<< "# total number of subspaces: " << num_ssp << "\n"
		<< "# total number of basins of attraction: " << num_bsns << "\n\n"
		<< "SSp obj BoA\n";
	for (size_t i = 0; i < num_ssp; ++i) {
		out_sstream << i << " " << obj[i] << " ";
		if (aff_bsn[i].empty())
			throw Exception("Should be affiliated to at least one basin of attraction");
		while (aff_bsn[i].size() > 1) {
			out_sstream << aff_bsn[i].front() << ",";
			aff_bsn[i].pop_front();
		}
		out_sstream << aff_bsn[i].front();
		out_sstream << "\n";
	}

	auto file_name = g_working_directory + "/result/boa/" + pro_name + "(" + std::to_string(num_vars) + "D).dat";
	std::ofstream out_fstream(file_name);
	out_fstream << out_sstream.str();
	out_fstream.close();

	DEL_PARAM(id_param);
	DEL_PRO(pro);
	DEL_RND(rnd);
}

void generateSamples(Problem *pro, size_t num_div) {
	std::vector<std::unique_ptr<Solution<>>> tmp_samples;

	size_t number_objectives = pro->numberObjectives();
	size_t num_cons = pro->numberConstraints();
	size_t num_vars = pro->numberVariables();

	std::string file_path = g_working_directory + "/result/samples/"
		+ pro->name() + '(' + std::to_string(num_vars) + "D)"
		+ '_' + std::to_string(num_div) + "_use_rank_as_obj.dat";
	std::ifstream in_file(file_path);
	if (in_file.fail()) {
		auto boundary = CAST_CONOP(pro)->boundary();
		if (num_vars == 2) {
			tmp_samples.resize(num_div * num_div);
			for (auto &sol : tmp_samples)
				sol.reset(new Solution<>(number_objectives, num_cons, num_vars));
			ofec::Real x1;
			for (size_t i = 0; i < num_div; ++i) {
				x1 = boundary[0].first + i * (boundary[0].second - boundary[0].first) / (num_div - 1);
				for (size_t j = 0; j < num_div; ++j)
					tmp_samples[j * num_div + i]->variable()[0] = x1;
			}
			ofec::Real x2;
			for (size_t i = 0; i < num_div; ++i) {
				x2 = boundary[1].first + i * (boundary[1].second - boundary[1].first) / (num_div - 1);
				for (size_t j = 0; j < num_div; ++j)
					tmp_samples[i * num_div + j]->variable()[1] = x2;
			}
		}
		if (!tmp_samples.empty()) {
			for (auto &sol : tmp_samples)
				sol->evaluate(pro, -1, false);
			std::vector<std::vector<Real>*> data(tmp_samples.size());
			for (size_t i = 0; i < tmp_samples.size(); ++i)
				data[i] = &(tmp_samples[i]->objective());
			std::vector<int> rank;
			nd_sort::filterSort<Real>(data, rank, pro->optimizeMode());

			std::stringstream out_stream;
			int flag_blank_line = 0;
			out_stream << std::fixed << std::setprecision(3);
			for (size_t i = 0; i < tmp_samples.size(); ++i) {
				for (size_t j = 0; j < num_vars; ++j)
					out_stream << std::setw(10) << tmp_samples[i]->variable()[j] << ' ';
				out_stream << std::setw(10) << rank[i] << '\n';
				if (num_vars == 2 && ++flag_blank_line == num_div) {
					out_stream << '\n';
					flag_blank_line = 0;
				}
			}
			std::ofstream out_file(file_path);
			out_file << out_stream.str();
		}
	}
}

void splitSubspace(size_t id_ssp,
	nanoflann::KDTreeSpace<Real> &ssp_tree,
	std::vector<std::list<size_t>> &neighbor)
{
	size_t id_ssp_new = ssp_tree.splitRegion(id_ssp);

	neighbor.push_back(neighbor[id_ssp]);
	auto &adj_ssp = neighbor[id_ssp];
	for (auto it = adj_ssp.begin(); it != adj_ssp.end();) {
		if (!ssp_tree.checkAdjacency(*it, id_ssp)) {
			neighbor[*it].remove(id_ssp);
			it = adj_ssp.erase(it);
		}
		else
			it++;
	}
	adj_ssp.push_back(id_ssp_new);
	auto &adj_ssp_new = neighbor[id_ssp_new];
	for (auto it = adj_ssp_new.begin(); it != adj_ssp_new.end();) {
		if (ssp_tree.checkAdjacency(*it, id_ssp_new)) {
			neighbor[*it].push_back(id_ssp_new);
			it++;
		}
		else
			it = adj_ssp_new.erase(it);
	}
	adj_ssp_new.push_back(id_ssp);
}

void splitBySols(
	size_t num_vars,
	const std::vector<std::vector<Real>> &centers,
	nanoflann::KDTreeSpace<Real> &ssp_tree,
	std::vector<std::list<size_t>> &neighbor)
{
	std::vector<std::vector<Real>> max_size(centers.size(), std::vector<Real>(num_vars, -1));
	std::vector<std::vector<std::vector<Real>>> dist_each_dim(num_vars);
	for (size_t j = 0; j < num_vars; ++j)
		dist_each_dim[j].resize(centers.size(), std::vector<Real>(centers.size()));
	for (size_t j = 0; j < num_vars; ++j) {
		for (size_t i = 0; i < centers.size(); ++i) {
			for (size_t k = i + 1; k < centers.size(); ++k) {
				dist_each_dim[j][i][k] = dist_each_dim[j][k][i] =
					fabs(centers[i][j] - centers[k][j]);
			}
		}
	}
	size_t dim_max;
	for (size_t i = 0; i < centers.size(); ++i) {
		for (size_t k = i + 1; k < centers.size(); ++k) {
			dim_max = 0;
			for (size_t j = 1; j < num_vars; ++j) {
				if (dist_each_dim[j][i][k] > dist_each_dim[dim_max][i][k])
					dim_max = j;
			}
			if (max_size[i][dim_max] == -1 || max_size[i][dim_max] > dist_each_dim[dim_max][i][k])
				max_size[i][dim_max] = dist_each_dim[dim_max][i][k];
			if (max_size[k][dim_max] == -1 || max_size[k][dim_max] > dist_each_dim[dim_max][i][k])
				max_size[k][dim_max] = dist_each_dim[dim_max][i][k];
		}
	}
	for (size_t i = 0; i < max_size.size(); ++i) {
		for (size_t j = 0; j < num_vars; ++j) {
			if (max_size[i][j] != -1)
				max_size[i][j] = max_size[i][j] / 2;
		}
	}
	for (size_t i = 0; i < centers.size(); ++i) {
		bool split;
		while (true) {
			split = false;
			auto id_ssp = ssp_tree.getRegionIdx(centers[i]);
			auto &box = ssp_tree.getBox(id_ssp);
			for (size_t j = 0; j < num_vars; ++j) {
				if (max_size[i][j] != -1 && (box[j].second - box[j].first) >= max_size[i][j]) {
					split = true;
					break;
				}
			}
			if (split)
				splitSubspace(id_ssp, ssp_tree, neighbor);
			else
				break;
		};
	}
}

void testGSP(Problem *pro, Random *rnd, const std::vector<std::vector<Real>> &centers) {

	if (centers.size() == 1) {
		std::cout << 1 << ' ';
		return;
	}

	bool is_minimize = pro->optimizeMode(0) == OptimizeMode::kMinimize;
	size_t num_vars = pro->numberVariables();
	size_t num_div = 3;
	size_t num_ssp;

	if (num_vars > 2) {
		std::cout << "--" << ' ';
		return;
	}

	while (true) {
		size_t id_ssp, num_bsns;
		if (num_vars == 1) {
			Real gap_x1, l_x1;
			num_ssp = num_div;
			auto &range_x1 = CAST_CONOP(pro)->range(0);
			gap_x1 = (range_x1.second - range_x1.first) / num_div;
			std::vector<Real> obj(num_ssp, (is_minimize ? std::numeric_limits<Real>::max() : -std::numeric_limits<Real>::max()));
			Solution<> s(1, 0, num_vars);
			for (size_t i = 0; i < num_div; ++i) {
				l_x1 = range_x1.first + i * gap_x1;
				id_ssp = i;
				for (size_t k = 0; k < 100; ++k) {
					s.variable()[0] = rnd->uniform.nextNonStd(l_x1, l_x1 + gap_x1);
					s.evaluate(pro, -1, false);
					if (is_minimize) {
						if (s.objective(0) < obj[id_ssp])
							obj[id_ssp] = s.objective(0);
					}
					else if (s.objective(0) > obj[id_ssp])
						obj[id_ssp] = s.objective(0);;
				}
			}
			std::vector<std::list<size_t>> neighbor(num_ssp), aff_bsn;
			for (size_t i = 0; i < 1; ++i) {
				id_ssp = i;
				neighbor[id_ssp] = { id_ssp + 1 };
			}
			for (size_t i = 1; i < num_div - 1; ++i) {
				id_ssp = i;
				neighbor[id_ssp] = { id_ssp - 1, id_ssp + 1 };
			}
			for (size_t i = num_div - 1; i < num_div; ++i) {
				id_ssp = i;
				neighbor[id_ssp] = { id_ssp - 1 };
			}
			groupSubspace(obj, neighbor, aff_bsn, is_minimize, num_bsns);

			std::vector<size_t> num_center(num_bsns, 0);
			size_t row;
			for (auto &c : centers) {
				row = floor((c[0] - range_x1.first) / gap_x1);
				if (row == num_div)
					id_ssp = row - 1;
				else
					id_ssp = row;
				for (auto b : aff_bsn[id_ssp])
					num_center[b]++;
			}
			bool separated = true;
			for (size_t b = 0; b < num_bsns; ++b) {
				if (num_center[b] > 1) {
					separated = false;
					break;
				}
			}
			if (separated)
				break;
			else
				num_div++;
		}
		else if (num_vars == 2) {
			Real gap_x1, gap_x2, l_x1, l_x2;
			num_ssp = num_div * num_div;
			auto &range_x1 = CAST_CONOP(pro)->range(0);
			auto &range_x2 = CAST_CONOP(pro)->range(1);
			gap_x1 = (range_x1.second - range_x1.first) / num_div;
			gap_x2 = (range_x2.second - range_x2.first) / num_div;

			std::vector<Real> obj(num_ssp, (is_minimize ? std::numeric_limits<Real>::max() : -std::numeric_limits<Real>::max()));
			Solution<> s(1, 0, num_vars);
			for (size_t i = 0; i < num_div; ++i) {
				l_x1 = range_x1.first + i * gap_x1;
				for (size_t j = 0; j < num_div; ++j) {
					l_x2 = range_x2.first + j * gap_x2;
					id_ssp = i * num_div + j;
					for (size_t k = 0; k < 100; ++k) {
						s.variable()[0] = rnd->uniform.nextNonStd(l_x1, l_x1 + gap_x1);
						s.variable()[1] = rnd->uniform.nextNonStd(l_x2, l_x2 + gap_x2);
						s.evaluate(pro, -1, false);
						if (is_minimize) {
							if (s.objective(0) < obj[id_ssp])
								obj[id_ssp] = s.objective(0);
						}
						else if (s.objective(0) > obj[id_ssp])
							obj[id_ssp] = s.objective(0);;
					}
				}
			}

			std::vector<std::list<size_t>> neighbor(num_ssp), aff_bsn;
			for (size_t i = 0; i < 1; ++i) {
				for (size_t j = 0; j < 1; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp + 1, id_ssp + num_div, id_ssp + num_div + 1 };
				}
				for (size_t j = 1; j < num_div - 1; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp - 1, id_ssp + 1, id_ssp + num_div - 1,id_ssp + num_div, id_ssp + num_div + 1 };
				}
				for (size_t j = num_div - 1; j < num_div; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp - 1, id_ssp + num_div - 1, id_ssp + num_div };
				}
			}
			for (size_t i = 1; i < num_div - 1; ++i) {
				for (size_t j = 0; j < 1; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp - num_div, id_ssp - num_div + 1, id_ssp + 1, id_ssp + num_div, id_ssp + num_div + 1 };
				}
				for (size_t j = 1; j < num_div - 1; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp - num_div - 1, id_ssp - num_div, id_ssp - num_div + 1, id_ssp - 1, id_ssp + 1, id_ssp + num_div - 1,id_ssp + num_div, id_ssp + num_div + 1 };
				}
				for (size_t j = num_div - 1; j < num_div; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp - num_div - 1, id_ssp - num_div,  id_ssp - 1, id_ssp + num_div - 1, id_ssp + num_div };
				}
			}
			for (size_t i = num_div - 1; i < num_div; ++i) {
				for (size_t j = 0; j < 1; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp - num_div, id_ssp - num_div + 1, id_ssp + 1 };
				}
				for (size_t j = 1; j < num_div - 1; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp - num_div - 1, id_ssp - num_div, id_ssp - num_div + 1, id_ssp - 1, id_ssp + 1 };
				}
				for (size_t j = num_div - 1; j < num_div; ++j) {
					id_ssp = i * num_div + j;
					neighbor[id_ssp] = { id_ssp - num_div - 1, id_ssp - num_div,  id_ssp - 1 };
				}
			}

			groupSubspace(obj, neighbor, aff_bsn, is_minimize, num_bsns);

			std::vector<size_t> num_center(num_bsns, 0);
			size_t row, col;
			for (auto &c : centers) {
				row = floor((c[0] - range_x1.first) / gap_x1);
				col = floor((c[1] - range_x2.first) / gap_x2);
				id_ssp = row * num_div + col;
				for (auto b : aff_bsn[id_ssp])
					num_center[b]++;
			}
			bool separated = true;
			for (size_t b = 0; b < num_bsns; ++b) {
				if (num_center[b] > 1) {
					separated = false;
					break;
				}
			}
			if (separated)
				break;
			else
				num_div++;
		}
	}

	std::cout << num_ssp << ' ';
}

void testBSP(Problem *pro, Random *rnd, const std::vector<std::vector<Real>> &centers) {

	bool is_minimize = pro->optimizeMode(0) == OptimizeMode::kMinimize;
	size_t num_vars = pro->numberVariables();

	nanoflann::KDTreeSpace<Real> ssp_tree;
	ssp_tree.setInitBox(CAST_CONOP(pro)->boundary());
	ssp_tree.inputRatioData({ 1.0 });
	ssp_tree.buildIndex();

	std::vector<std::list<size_t>> neighbor(1), aff_bsn;
	splitBySols(num_vars, centers, ssp_tree, neighbor);
	std::vector<Real> obj(ssp_tree.size(), (is_minimize ? std::numeric_limits<Real>::max() : -std::numeric_limits<Real>::max()));
	Solution<> s(1, 0, num_vars);
	for (size_t id_ssp = 0; id_ssp < obj.size(); ++id_ssp) {
		auto &box = ssp_tree.getBox(id_ssp);
		for (size_t k = 0; k < 100; ++k) {
			for (size_t j = 0; j < num_vars; ++j)
				s.variable()[j] = rnd->uniform.nextNonStd(box[j].first, box[j].second);
			s.evaluate(pro, -1, false);
			if (is_minimize) {
				if (s.objective(0) < obj[id_ssp])
					obj[id_ssp] = s.objective(0);
			}
			else if (s.objective(0) > obj[id_ssp])
				obj[id_ssp] = s.objective(0);;
		}
	}
	for (size_t id_opt = 0; id_opt < pro->optima().number(); ++id_opt) {
		for (size_t j = 0; j < num_vars; ++j)
			s.variable()[j] = pro->optima().variable(id_opt)[j];
		s.evaluate(pro, -1, false);
		size_t id_ssp = ssp_tree.getRegionIdx(s.variable().vect());
		if (is_minimize) {
			if (s.objective(0) < obj[id_ssp])
				obj[id_ssp] = s.objective(0);
		}
		else if (s.objective(0) > obj[id_ssp])
			obj[id_ssp] = s.objective(0);;
	}
	//std::vector<std::vector<Real>> vec_obj(obj.size(), std::vector<Real>(1));
	//for (size_t i = 0; i < obj.size(); ++i)
	//	vec_obj[i][0] = obj[i];
	//std::vector<std::vector<Real> *> data(vec_obj.size());
	//for (size_t i = 0; i < vec_obj.size(); ++i)
	//	data[i] = &(vec_obj[i]);
	//std::vector<int> rank;
	//int num_ranks = nd_sort::filterSort<Real>(data, rank, pro->optimizeMode());
	size_t num_bsns;
	groupSubspace(obj, neighbor, aff_bsn, is_minimize, num_bsns);

	std::vector<size_t> num_center(num_bsns, 0);
	for (auto &c : centers) {
		size_t id_ssp = ssp_tree.getRegionIdx(c);
		for (auto b : aff_bsn[id_ssp])
			num_center[b]++;
	}
	bool separated = true;
	for (size_t b = 0; b < num_bsns; ++b) {
		if (num_center[b] > 1) {
			separated = false;
			break;
		}
	}

	std::cout << (separated ? "true" : "false") << " " << ssp_tree.size() << " ";
}

void compareGridAndKDTree(const std::string &pro_name, int num_vars) {
	ParameterMap v;
	v["problem name"] = pro_name;
	v["number of variables"] = num_vars;
	int id_param = std::make_shared<const ParameterMap>(v);
	Problem *pro = Problem::generateByFactory(id_param, 0.5);
	pro->initialize();
	Random *rnd = ADD_RND(0.5);

	//generateSamples(pro, 30);

	std::vector<std::vector<Real>> centers;
	for (size_t i = 0; i < pro->optima().number(); ++i)
		centers.push_back(pro->optima().variable(i).vect());

	std::cout << pro_name << ' ' << num_vars << ' ' << centers.size() << ' ';

	testBSP(pro, rnd, centers);

	testGSP(pro, rnd, centers);

	std::cout << std::endl;

	DEL_PARAM(id_param);
	DEL_PRO(pro);
	DEL_RND(rnd);
}

TEST_CASE("Generate BoAs", "[generate][BoA]") {
	generateInfoBoA("Classic_sphere", 2, 5);
	generateInfoBoA("Classic_six_hump_camel_back", 2, 100);
	generateInfoBoA("Classic_Vincent", 2, 200);
	generateInfoBoA("Classic_Rastrigin_modified", 2, 100);
	generateInfoBoA("Classic_Himmenblau", 2, 100);
}

TEST_CASE("Compared grid and k-d tree", "[BSP][GSP]") {
	compareGridAndKDTree("MMOP_CEC2013_F01", 1);
	compareGridAndKDTree("MMOP_CEC2013_F02", 1);
	compareGridAndKDTree("MMOP_CEC2013_F03", 1);
	compareGridAndKDTree("MMOP_CEC2013_F04", 2);
	compareGridAndKDTree("MMOP_CEC2013_F05", 2);
	compareGridAndKDTree("MMOP_CEC2013_F06", 2);
	compareGridAndKDTree("MMOP_CEC2013_F07", 2);
	compareGridAndKDTree("MMOP_CEC2013_F06", 3);
	compareGridAndKDTree("MMOP_CEC2013_F07", 3);
	compareGridAndKDTree("MMOP_CEC2013_F08", 2);
	compareGridAndKDTree("MMOP_CEC2013_F09", 2);
	compareGridAndKDTree("MMOP_CEC2013_F10", 2);
	compareGridAndKDTree("MMOP_CEC2013_F11", 2);
	compareGridAndKDTree("MMOP_CEC2013_F11", 3);
	compareGridAndKDTree("MMOP_CEC2013_F12", 3);
	compareGridAndKDTree("MMOP_CEC2013_F11", 5);
	compareGridAndKDTree("MMOP_CEC2013_F12", 5);
	compareGridAndKDTree("MMOP_CEC2013_F11", 10);
	compareGridAndKDTree("MMOP_CEC2013_F12", 10);
	compareGridAndKDTree("MMOP_CEC2013_F12", 20);
}
