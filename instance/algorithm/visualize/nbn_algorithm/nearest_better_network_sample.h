#ifndef NEAREST_BETTER_NETWORK_SAMPLE_H
#define NEAREST_BETTER_NETWORK_SAMPLE_H

#include "../../../core/problem/solution.h"
#include "../../../utility/random/newran.h"
#include "../../../utility/function/custom_function.h"
#include "../../../core/problem/problem.h"

#include<memory>


namespace ofec {
	class NearestBetterNetworkSample {
	protected:
		std::vector<std::shared_ptr<SolutionBase>> m_samples;
		std::vector<int> m_parent_ids;
		std::vector<bool> m_opt_flag;
		std::vector<double> m_min_dis;
		//std::vector<double> m_fitness;
		int m_random.get();
		int m_problem.get();
	public:

		const std::vector<std::shared_ptr<SolutionBase>>& Sample() const {
			return m_samples;
		}
		void clear() {
			m_random.get() = -1;
			m_problem.get() = -1;
			m_samples.clear();
			m_parent_ids.clear();
			//m_fitness.clear();
		}
		void initialize(Random *rnd, Problem *pro) {
			clear();
			m_random.get() = rnd;
			m_problem.get() = pro;
		}

		void getNearestBetterNetwork(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<int>& marker
			) {
			marker.clear();
			sols = m_samples;
			belong = m_parent_ids;
			fitness.resize(sols.size());
			for (int idx(0); idx < belong.size(); ++idx) {
				fitness[idx] = m_samples[idx]->fitness();
				if (m_opt_flag[idx]) {
					marker.push_back(idx);
				}
			}
			
		}



		void getNearestBetterNetwork(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong) {

			sols = m_samples;
			belong = m_parent_ids;
			fitness.resize(sols.size());
			for (int idx(0); idx < belong.size(); ++idx) {
				fitness[idx] = m_samples[idx]->fitness();
			}
		}


		void initSols(std::vector<std::shared_ptr<SolutionBase>>& sols,
			std::vector<int>& belong) {
			m_samples = sols;
			m_parent_ids = belong;
			m_opt_flag.resize(m_samples.size(), false);
			m_min_dis.resize(m_samples.size(), std::numeric_limits<double>::max());
			for (int idx(0); idx < m_samples.size(); ++idx) {
				if (m_parent_ids[idx] != idx) {
					m_min_dis[idx] = m_samples[idx]->variableDistance(*m_samples[m_parent_ids[idx]], m_problem.get());
				}
			}
		}
		void initRandomSol(int initNum, const std::function<void(SolutionBase& sol, Problem *pro)>& eval_fun) {


			if (m_problem->hasTag(ProblemTag::kCSIWDN)) {
				SolutionBase* ptr_cur_sol;
				for (int idx(0); idx < initNum; ++idx) {
					ptr_cur_sol = m_problem->createSolution();
					m_problem->initializeSolution(*ptr_cur_sol, m_random.get());
					eval_fun(*ptr_cur_sol, m_problem.get());
					std::shared_ptr<SolutionBase> cur_sol;
					cur_sol.reset(ptr_cur_sol);
					addRandomSol(cur_sol);
				}
			}
			else {
				std::vector<SolutionBase*> sols(initNum);
				UTILITY::generateRandomSolutionsMultiThreads(sols, GET_RND(m_random.get()), m_problem.get(), eval_fun);
				for (int idx(0); idx < sols.size(); ++idx) {
					std::shared_ptr<SolutionBase> cur_sol;
					cur_sol.reset(sols[idx]);
					addRandomSol(cur_sol);
				}
				sols.clear();
			}
		
			
			//sols.clear();
		}

		void addOptSols(const std::function<void(SolutionBase& sol, Problem *pro)>& eval_fun) {
			int num_vars = m_problem->numberVariables();
			int number_objectives = m_problem->numberObjectives();
			int num_cons = m_problem->numberConstraints();

			auto& optBase(m_problem->optBase());
			for (size_t idx(0); idx < optBase.numberVariables(); ++idx) {
				std::shared_ptr<SolutionBase> cur_sol = nullptr;
				cur_sol.reset(m_problem->createSolution(optBase.variable(idx)));
				//cur_sol.reset(new SolutionType(number_objectives, num_cons, num_vars));
				//cur_sol->variable() = dynamic_cast<const SolutionType&>(optBase.variable(idx)).variable();
				eval_fun(*cur_sol, m_problem.get());
				addRandomSol(cur_sol, true);
			}
		}

		void addRandomSol(std::shared_ptr<SolutionBase>& sol, bool flag_opt = false) {
			int cur_idx(m_samples.size());
			m_samples.push_back(sol);
			m_parent_ids.push_back(cur_idx);
			m_opt_flag.push_back(flag_opt);
			m_min_dis.push_back(std::numeric_limits<double>::max());
			for (int idx(0); idx < cur_idx; ++idx) {
				if (!m_opt_flag[idx]) {
					auto& other(m_samples[idx]);
					if (m_opt_flag[cur_idx] || other->fitness() < sol->fitness()) {
						double dis = other->variableDistance(*sol, m_problem.get());
						if (m_min_dis[idx] > dis) {
							m_parent_ids[idx] = cur_idx;
							m_min_dis[idx] = dis;
						}
						else if (m_min_dis[idx] == dis) {
							if (m_random->uniform.next() < 0.5) {
								m_parent_ids[idx] = cur_idx;
							}
						}
					}
					else if (other->fitness() > sol->fitness()) {
						double dis = other->variableDistance(*sol, m_problem.get());
						if (m_min_dis[cur_idx] > dis) {
							m_min_dis[cur_idx] = dis;
							m_parent_ids[cur_idx] = idx;
						}
						else if (m_min_dis[cur_idx] == dis) {
							if (m_random->uniform.next() < 0.5) {
								m_parent_ids[cur_idx] = idx;
							}
						}
					}

				}
			}
		}


	};
}


#endif 