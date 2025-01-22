#include "hghe_de.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../utility/clustering/nbc.h"
#include "../../../../../../utility/nondominated_sorting/fast_sort.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	void HGHE_DE::addInputParameters() {
		m_input_parameters.add("scaling factor", new RangedReal(m_scaling_factor, 0, 2, 0.5));
		m_input_parameters.add("crossover rate", new RangedReal(m_crossover_rate, 0, 1, 0.6));
		m_input_parameters.add("size of exploitive subpopulation", new RangedSizeT(
			m_size_subpop_exploit, 5, 50, 15));
		m_input_parameters.add("number of solutions for exploration", new RangedSizeT(
			m_num_sols_explore, 1, 200, 10));
		m_input_parameters.add("potential prediction strategy", new Enumeration(
			m_ptnl_pred_strat, 
			{ "mean", "random", "Pareto", "obj first", "frequecy first", "number of optima" },
			PtnlPredStrat::kMean));
	}
	
	void HGHE_DE::initialize_(Environment *env) {
		HGHEC::initialize_(env);
		m_subpops_exploit.clear();
		m_sols_explore.clear();
	}

	void HGHE_DE::run_(Environment *env) {
		initHills(env);
		while (!terminating()) {
			updateHills(env);
			if (classifyPopulation(env) & kTerminate)
				break;
			updatePotentialExplore(env);
			if (exploreHills(env) & kTerminate)
				break;
			if (exploitHills(env) & kTerminate)
				break;
#ifdef OFEC_DATUM_MULTI_POP_H
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(1 + m_subpops_exploit.size());
			for (size_t i = 0; i < m_num_sols_explore; ++i) {
				g_multi_pop.pops[0].push_back(m_sols_explore[i].get());
			}
			for (size_t k = 0; k < m_subpops_exploit.size(); ++k) {
				for (size_t i = 0; i < m_subpops_exploit[k].size(); ++i) {
					g_multi_pop.pops[k + 1].push_back(&m_subpops_exploit[k][i]);
				}
			}
			datumUpdated(env, g_multi_pop);
#endif
		}
	}

	int HGHE_DE::classifyPopulation(Environment *env) {
		/* Find the hill region to which each subpopulation belongs */
		std::vector<std::list<size_t>> subpop_in_each_hill(m_hills.size());
		for (size_t k = 0; k < m_subpops_exploit.size(); ++k) {
			if (m_subpops_exploit[k].size() == 0)
				continue;
			std::vector<size_t> freq_each_hill(m_hills.size(), 0);
			for (size_t i = 0; i < m_subpops_exploit[k].size(); ++i) {
				for (Hill *hill : theHillsLocated(m_subpops_exploit[k][i])) {
					freq_each_hill[m_ptr_to_id[hill]]++;
				}
			}
			size_t id_max = std::max_element(freq_each_hill.begin(), freq_each_hill.end()) - freq_each_hill.begin();
			subpop_in_each_hill[id_max].push_back(k);
		}
		/* Merge subpopulations belonging to the same hill region */
		decltype(m_subpops_exploit) tmp_subpops_exploit;
		m_subpops_exploit.moveTo(tmp_subpops_exploit);
		m_subpops_exploit.resize(subpop_in_each_hill.size());
		for (size_t k = 0; k < m_subpops_exploit.size(); ++k) {
			m_subpops_exploit[k].setParameter(m_crossover_rate, m_scaling_factor);
			for (size_t tmp_k : subpop_in_each_hill[k]) {
				m_subpops_exploit[k].append(tmp_subpops_exploit[tmp_k]);
			}
		}
		/* According to the hill region to which each subpopulation belongs,
		removeSubspace the Solutions that are out of the region and make up for the shortage of Solutions */
		for (size_t k = 0; k < m_subpops_exploit.size(); ++k) {
			for (auto it = m_subpops_exploit[k].begin(); it != m_subpops_exploit[k].end();) {
				if (isSolInHill(**it, m_id_to_ptr[k]))
					it++; 
				else
					it = m_subpops_exploit[k].remove(it);
			}
			while (m_subpops_exploit[k].size() > m_size_subpop_exploit) {
				auto it = m_subpops_exploit[k].begin(), it_worst = it;
				while (++it != m_subpops_exploit[k].end()) {
					if (dominate(**it_worst, **it, env->problem()->optimizeMode()))
						it_worst = it;
				}
				m_subpops_exploit[k].remove(it_worst);
			}
			while (m_subpops_exploit[k].size() < m_size_subpop_exploit) {
				auto new_ind = std::make_unique<IndividualDE>(env);
				randomSolInHill(*new_ind, m_id_to_ptr[k]);
				int rf = new_ind->evaluate(env);
				if (rf != kNormalEval)
					return rf;
				archiveSolution(*new_ind, TaskEval::kExplore, env);
				m_subpops_exploit[k].append(new_ind);
			}
		}
		return kNormalEval;
	}

	void HGHE_DE::updatePotentialExplore(Environment *env) {
		if (m_hills.size() > 1) {
			if (m_ptnl_pred_strat == PtnlPredStrat::kMean) {
				for (auto &hill : m_hills)
					hill->potential_explore = 1.0;
			}
			else if (m_ptnl_pred_strat == PtnlPredStrat::kRandom) {
				for (auto &hill : m_hills)
					hill->potential_explore = m_random->uniform.next();
			}
			else if (m_ptnl_pred_strat == PtnlPredStrat::kNumOpts) { /* god mode, only for test */
				auto optima = CAST_CONOP(env->problem())->optima();
				for (auto &hill : m_hills) {
					hill->potential_explore = 0;
					for (size_t i = 0; i < optima->numberSolutions(); ++i) {
						if (isVarInHill(optima->solution(i).variable(), hill.get())) {
							hill->potential_explore++;
						}
					}
				}
				bool all_zero = true;
				for (auto &hill : m_hills) {
					if (hill->potential_explore < 2)
						hill->potential_explore = 0;
					else {
						all_zero = false;
						hill->potential_explore--;
					}
				}
				if (all_zero) {
					for (auto &hill : m_hills)
						hill->potential_explore = 1.0;
				}
			}
			else {
				std::vector<std::vector<Real>> attributes(m_hills.size(), std::vector<Real>(2));
				for (auto &hill : m_hills) {
					const SolutionType *best_sol = nullptr;
					attributes[m_ptr_to_id[hill.get()]][1] = 0;										// exploration frequency
					for (Subspace* subspace : hill->subspaces()) {
						auto cur_best = subspace->best_sol;
						if (cur_best && (!best_sol || dominate(*cur_best, *best_sol, env->problem()->optimizeMode())))
							best_sol = cur_best;
						attributes[m_ptr_to_id[hill.get()]][1] += subspace->his_explore.size();
					}
					if (best_sol)
						attributes[m_ptr_to_id[hill.get()]][0] = best_sol->objective(0);		// best objective
					else
						throw Exception("A hill has no solutions evaluted");
				}

				Real min_obj, max_obj, min_freq, max_freq;
				min_obj = max_obj = attributes[0][0];
				min_freq = max_freq = attributes[0][1];
				for (size_t i = 1; i < attributes.size(); ++i) {
					if (min_obj > attributes[i][0])
						min_obj = attributes[i][0];
					if (max_obj < attributes[i][0])
						min_obj = attributes[i][0];
					if (min_freq > attributes[i][1])
						min_freq = attributes[i][1];
					if (max_freq < attributes[i][1])
						max_freq = attributes[i][1];
				}

				Real obj_step = 0.01;
				for (size_t i = 0; i < attributes.size(); ++i) {
					attributes[i][0] = mapReal(attributes[i][0], min_obj, max_obj, 0.0, 1.0);
					attributes[i][0] = int(attributes[i][0] / obj_step) * obj_step;
					attributes[i][1] = mapReal(attributes[i][1], min_freq, max_freq, 0.0, 1.0);
				}

				std::vector<int> rank(attributes.size());
				int num_rank = 0;

				if (m_ptnl_pred_strat == PtnlPredStrat::kObjFirst) {
					std::vector<Real> data(attributes.size());
					for (size_t i = 0; i < attributes.size(); ++i)
						data[i] = attributes[i][0];
					std::vector<int> sequence;
					mergeSort(data, attributes.size(), sequence, env->problem()->optimizeMode()[0] == OptimizeMode::kMinimize);
					for (size_t i = 0; i < attributes.size(); ++i)
						rank[sequence[i]] = i;
					num_rank = attributes.size();
				}
				else if (m_ptnl_pred_strat == PtnlPredStrat::kFreqFirst) {
					std::vector<Real> data(attributes.size());
					for (size_t i = 0; i < attributes.size(); ++i)
						data[i] = attributes[i][1];
					std::vector<int> sequence;
					mergeSort(data, attributes.size(), sequence, true);
					for (size_t i = 0; i < attributes.size(); ++i)
						rank[sequence[i]] = i;
					num_rank = attributes.size();
				}
				else if (m_ptnl_pred_strat == PtnlPredStrat::kPareto) {
					std::vector<std::vector<Real> *> data(attributes.size());
					for (size_t i = 0; i < attributes.size(); ++i)
						data[i] = &attributes[i];
					std::vector<OptimizeMode> opt_mode = { env->problem()->optimizeMode()[0], OptimizeMode::kMinimize };
					num_rank = nd_sort::fastSort(data, rank, opt_mode);
				}

				for (size_t i = 0; i < attributes.size(); i++)
					m_id_to_ptr[i]->potential_explore = pow((num_rank - rank[i]) / (Real)num_rank, 2);
			}
		}
		else {
			m_hills.front()->potential_explore = 1.0;
		}
		m_sum_potential_explore = 0;
		for (auto &hill : m_hills)
			m_sum_potential_explore += hill->potential_explore;
	}

	int HGHE_DE::exploreHills(Environment *env) {
		while (m_sols_explore.size() < m_num_sols_explore) {
			m_sols_explore.emplace_back(dynamic_cast<SolutionType*>(env->problem()->createSolution()));
		}
		for (size_t i = 0; i < m_num_sols_explore; ++i) {
			Hill *hill = rouletteWheelSelection();
			randomVarInHill(m_sols_explore[i]->variable(), hill);
			int rf = m_sols_explore[i]->evaluate(env);
			if (rf != kNormalEval)
				break;
			archiveSolution(*m_sols_explore[i], TaskEval::kExplore, env);
		}
		return kNormalEval;
	}

	int HGHE_DE::exploitHills(Environment *env) {
		int rf = kNormalEval;
		for (size_t k = 0; k < m_subpops_exploit.size(); ++k) {
			if (m_subpops_exploit[k].isActive()) {
				size_t best = 0, worst = 0;
				for (size_t i = 1; i < m_subpops_exploit[k].size(); ++i) {
					if (dominate(m_subpops_exploit[k][i], m_subpops_exploit[k][best], env->problem()->optimizeMode()))
						best = i;
					if (dominate(m_subpops_exploit[k][worst], m_subpops_exploit[k][i], env->problem()->optimizeMode()))
						worst = i;
				}
				if (m_subpops_exploit[k][best].objectiveDistance(m_subpops_exploit[k][worst]) < 1e-6)
					m_subpops_exploit[k].setActive(false);
			}
			if (!m_subpops_exploit[k].isActive())
				continue;
			for (size_t i = 0; i < m_subpops_exploit[k].size(); ++i) {
				do {
					m_subpops_exploit[k].mutate(i, m_random.get(), env);
					m_subpops_exploit[k].recombine(i, m_random.get(), env);
					CAST_CONOP(env->problem())->validateSolution(m_subpops_exploit[k][i].trial(), Validation::kSetToBound, m_random.get());
				} while (!isSolInHill(m_subpops_exploit[k][i].trial(), m_id_to_ptr[k]));
				rf = m_subpops_exploit[k][i].select(env);
				if (rf != kNormalEval)
					return rf;
				archiveSolution(m_subpops_exploit[k][i].trial(), TaskEval::kExploit, env);
			}
		}
		return rf;
	}

	HGHE_DE::Hill* HGHE_DE::rouletteWheelSelection() {
		if (m_hills.empty())
			return nullptr;
		Real rand_pos = m_sum_potential_explore * m_random->uniform.next();
		Real accum = 0;
		for (auto &hill : m_hills) {
			accum += hill->potential_explore;
			if (rand_pos <= accum)
				return hill.get();
		}
		return nullptr;
	}
}