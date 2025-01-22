#include "cdg_moea.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/metricsMOP/IGD.h"
#include "../../../../record/rcr_vec_real.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	/******************************CDGMOEA***********************************/
	void CDGMOEA::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;;
		m_pop_size = v.get<int>("population size");
		m_pop.reset();
	}

	void CDGMOEA::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		std::vector<std::vector<Real>> ref_objs;
		for (size_t i = 0; i < m_problem->optimaBase()->numberObjectives(); ++i) {
			ref_objs.push_back(m_problem->optimaBase()->objective(i));
		}
		std::vector<std::vector<Real>> pop_objs;
		for (size_t i = 0; i < m_pop->size(); ++i) {
			pop_objs.push_back(m_pop->at(i).objective());
		}
		Real temp_IGD = IGD(ref_objs, pop_objs);
		//Real IGD = m_problem->optimaBase()->invertGenDist(*m_pop);
		entry.push_back(temp_IGD);
		dynamic_cast<RecordVectorReal*>(m_record.get())->addEntry(this, entry);
	}

#ifdef OFEC_DEMO
	void CDGMOEA::updateBuffer() {
		m_solution.clear();
		m_solution.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i)
			m_solution[0].push_back(&m_pop->at(i));
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif

	void CDGMOEA::initPop() {
		m_pop.reset(new PopCDGMOEA(m_pop_size, m_problem.get()));
		m_pop->initialize(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
	}

	void CDGMOEA::run_() {
		initPop();
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	/****************************** population ***********************************/
	PopCDGMOEA::PopCDGMOEA(size_t size_pop, Problem *pro) : PopMODE(size_pop, pro) {
		//m_offspring.resize(2 * size_pop);
		int numobj = pro->numberObjectives();
		if (numobj == 2) {
			m_grid_div = 180;
			m_T = 5;
		}
		else if (numobj == 3) {
			m_grid_div = 30;
			m_T = 1;
		}
		m_ideal.resize(numobj), m_nadir.resize(numobj), m_grid_distance.resize(numobj), m_grid_min_max.resize(numobj),
			m_S.resize(numobj, std::vector<std::vector<int>>(pow(m_grid_div, numobj - 1)));// m_R.resize(2 * size_pop, std::vector<int>(M));
		for (int i = 0; i < numobj; ++i) {
			m_ideal[i] = m_nadir[i] = m_individuals[0]->objective()[i];
		}
	}

	void PopCDGMOEA::initialize(Problem *pro, Random *rnd) {
		Population<IndCDGMOEA>::initialize(pro, rnd);
		setParameter(1.0, 0.5); //set parmeter: CR=1.0, F=0.5
		updateIdealPoint(pro);
		updateNadirPoint(pro);
		gridConstruct(pro);
		assign_GN(pro);
		for (auto &i : m_individuals) {
			m_offspring.emplace_back(*i);
		}
		for (auto &i : m_individuals) {
			m_offspring.emplace_back(*i);
		}
	}

	void PopCDGMOEA::updateIdealPoint(Problem *pro) {
		int numobj = pro->numberObjectives();
		for (int i = 0; i < numobj; ++i) {
			m_ideal[i] = 1.0e14;
			for (int j = 0; j < m_individuals.size(); ++j) {
				if (m_individuals[j]->objective()[i] < m_ideal[i])
					m_ideal[i] = m_individuals[j]->objective()[i];
			}
		}
	}

	void PopCDGMOEA::updateIdealPoint(const std::vector<IndCDGMOEA> &offspring, Problem *pro) {
		int numobj = pro->numberObjectives();
		for (int i = 0; i < numobj; ++i) {
			for (int j = 0; j < offspring.size(); ++j) {
				if (offspring[j].objective()[i] < m_ideal[i])
					m_ideal[i] = offspring[j].objective()[i];
			}
		}
	}

	void PopCDGMOEA::updateNadirPoint(Problem *pro) {
		int numobj = pro->numberObjectives();
		for (int i = 0; i < numobj; ++i) {
			m_nadir[i] = -1 * 1.0e14;
			for (int j = 0; j < m_individuals.size(); ++j) {
				if (m_nadir[i] < m_individuals[j]->objective()[i])
					m_nadir[i] = m_individuals[j]->objective()[i];
			}
		}
	}

	void PopCDGMOEA::updateNadirPoint(const std::vector<IndCDGMOEA> &offspring, Problem *pro) {  //
		int numobj = pro->numberObjectives();
		std::vector<int> SP;
		for (int i = 0; i < offspring.size(); ++i) {
			int flag = 1;
			for (int j = 0; j < numobj; ++j) {
				if (offspring[i].objective()[j] < (m_ideal[j] + m_nadir[j] / 5)) {
					SP.push_back(i);
					break;
				}
			}
		}
		if (SP.size() != 0) {	//find all nondominated solution
			std::vector<IndCDGMOEA> SP1;
			for (int i = 0; i < SP.size(); ++i) {
				SP1.push_back(offspring[SP[i]]);
			}
			nondominatedSorting(SP1, pro);
			int size = SP.size(), count = 0;
			std::vector<int>().swap(SP);	//free  SP
			for (int i = 0; i < size; i++) {
				if (SP1[i].fitness() == 0) {
					count++;
					SP.push_back(i);
				}
			}
			for (int i = 0; i < numobj; ++i) {
				m_nadir[i] = -1 * 1.0e14;
				for (int j = 0; j < SP.size(); ++j) {
					if (m_nadir[i] < SP1[SP[j]].objective()[i])
						m_nadir[i] = SP1[SP[j]].objective()[i];
				}
			}
		}
	}

	void PopCDGMOEA::gridConstruct(Problem *pro) {
		int numobj = pro->numberObjectives();
		for (int i = 0; i < numobj; ++i) {
			m_grid_distance[i] = (m_nadir[i] - m_ideal[i]) * (m_grid_div + 1) / (m_grid_div * m_grid_div);
			m_grid_min_max[i].first = m_ideal[i] - (m_nadir[i] - m_ideal[i]) / (2 * m_grid_div);
			m_grid_min_max[i].second = m_nadir[i] + (m_nadir[i] - m_ideal[i]) / (2 * m_grid_div);
		}
	}

	void PopCDGMOEA::assign_GN(Problem *pro) {
		int numobj = pro->numberObjectives();
		for (int i = 0; i < m_individuals.size(); i++) {
			for (int j = 0; j < numobj; j++) {
				m_individuals[i]->Gk()[j] = (int)ceil((m_individuals[i]->objective(j) - m_grid_min_max[j].first) / m_grid_distance[j]);
			}
		}
		for (int i = 0; i < m_individuals.size(); ++i) {
			for (int j = i + 1; j < m_individuals.size(); ++j) {
				int GD_max = abs(m_individuals[i]->Gk()[0] - m_individuals[j]->Gk()[0]);
				for (int k = 0; k < numobj; ++k) {
					if (GD_max < abs(m_individuals[i]->Gk()[k] - m_individuals[j]->Gk()[k]))
						GD_max = abs(m_individuals[i]->Gk()[k] - m_individuals[j]->Gk()[k]);
				}
				if (GD_max < m_T || GD_max == m_T) {
					m_individuals[i]->GN().push_back(j);
					m_individuals[j]->GN().push_back(i);
				}
			}
		}
	}

	void PopCDGMOEA::gridConstruct_assignGN_P_reserve(std::vector<IndCDGMOEA> &offspring, std::vector<int> &P_reserve, int size, Problem *pro) {
		int numobj = pro->numberObjectives();
		for (int i = 0; i < numobj; ++i) {
			m_grid_distance[i] = (m_nadir[i] - m_ideal[i]) * (m_grid_div + 1) / (m_grid_div * m_grid_div);
			m_grid_min_max[i].first = m_ideal[i] - (m_nadir[i] - m_ideal[i]) / (2 * m_grid_div);
			m_grid_min_max[i].second = m_nadir[i] + (m_nadir[i] - m_ideal[i]) / (2 * m_grid_div);
		}
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < numobj; j++) {
				offspring[P_reserve[i]].Gk()[j] = (int)ceil((offspring[P_reserve[i]].objective()[j] - m_grid_min_max[j].first) / m_grid_distance[j]);
			}
		}
		for (int i = 0; i < size; ++i) {
			for (int j = i + 1; j < size; ++j) {
				int GD_max = abs(offspring[P_reserve[i]].Gk()[0] - offspring[P_reserve[j]].Gk()[0]);
				for (int k = 0; k < numobj; ++k) {
					if (GD_max < abs(offspring[P_reserve[i]].Gk()[k] - offspring[P_reserve[j]].Gk()[k]))
						GD_max = abs(offspring[P_reserve[i]].Gk()[k] - offspring[P_reserve[i]].Gk()[k]);
				}
				if (GD_max < m_T || GD_max == m_T) {
					offspring[P_reserve[i]].GN().push_back(P_reserve[j]);
					offspring[P_reserve[j]].GN().push_back(P_reserve[i]);
				}
			}
		}

	}

	void PopCDGMOEA::assign_S(std::vector<IndCDGMOEA> &offspring, std::vector<int> P_reserve, Problem *pro) {
		int numobj = pro->numberObjectives();
		int size = pow(m_grid_div, numobj - 1);
		if (numobj == 2) {
			for (int i = 0; i < P_reserve.size(); ++i) {
				for (int l = 0; l < numobj; ++l) {
					for (int j = 0; j < numobj; ++j) {
						if (j != l) {
							for (int k = 0; k < m_grid_div; ++k) {	//M^(K-1)							
								if (offspring[P_reserve[i]].Gk()[j] == k + 1)
									m_S[l][k].push_back(P_reserve[i]);
							}
						}
					}
				}
			}
		}
		else if (numobj == 3) {
			for (int i = 0; i < P_reserve.size(); ++i) {
				for (int l = 0; l < numobj; ++l) {
					if (l == 0) {
						for (int k = 0; k < m_grid_div; ++k) {	//M^(K-1)
							for (int t = 0; t < m_grid_div; ++t) {
								if (offspring[P_reserve[i]].Gk()[1] == t + 1 && offspring[P_reserve[i]].Gk()[2] == k + 1)
									m_S[l][m_grid_div * k + t].push_back(P_reserve[i]);
							}
						}
					}
					else if (l == 1) {
						for (int k = 0; k < m_grid_div; ++k) {	//M^(K-1)
							for (int t = 0; t < m_grid_div; ++t) {
								if (offspring[P_reserve[i]].Gk()[0] == t + 1 && offspring[P_reserve[i]].Gk()[2] == k + 1)
									m_S[l][m_grid_div * k + t].push_back(P_reserve[i]);
							}
						}
					}
					else {
						for (int k = 0; k < m_grid_div; ++k) {	//M^(K-1)
							for (int t = 0; t < m_grid_div; ++t) {
								if (offspring[P_reserve[i]].Gk()[0] == t + 1 && offspring[P_reserve[i]].Gk()[1] == k + 1)
									m_S[l][m_grid_div * k + t].push_back(P_reserve[i]);
							}
						}
					}
				}
			}
		}
		/*test m_S*/
		for (int l = 0; l < numobj; ++l) {
			int num = 0;
			for (int k = 0; k < size; ++k) {
				num += m_S[l][k].size();
			}
			if (numobj == 2 && num != P_reserve.size())
				std::cout << "assign m_S is error:num=" << num << std::endl;
		}

	}

	void PopCDGMOEA::RBS(std::vector<IndCDGMOEA> &offspring, std::vector<int> P_reserve, Problem *pro) {
		int numobj = pro->numberObjectives();
		std::vector<std::vector<int>> R_(P_reserve.size(), std::vector<int>(numobj));
		/*Decomposition_based Ranking*/
		for (int l = 0; l < numobj; ++l) {
			for (int k = 0; k < pow(m_grid_div, numobj - 1); ++k) {
				std::vector<std::pair<double, int>> v;	//first store objective value, second store index
				for (int j = 0; j < m_S[l][k].size(); ++j) {
					v.push_back(std::make_pair(offspring[m_S[l][k][j]].objective()[l], m_S[l][k][j]));
				}
				if (v.size()) {
					sortl(v);	//sort of objective value according to the ascending
					for (int i = 0; i < m_S[l][k].size(); ++i) {
						offspring[v[i].second].R()[l] = i + 1;
					}
				}
			}
		}
		/*Lexicographic Sorting*/
		for (int i = 0; i < P_reserve.size(); ++i) {
			m_R.push_back(offspring[P_reserve[i]].R());
			sortR(m_R[i], R_[i]);
		}
		std::vector<int> number;	//to get index of P_reserve in offspring
		for (int i = 0; i < P_reserve.size(); ++i) {
			number.push_back(i);
		}
		sortL(R_, number, pro);
		for (int i = 0; i < m_individuals.size(); ++i) {
			*m_individuals[i] = offspring[P_reserve[number[i]]];
		}
		/*free the memory*/
		std::vector<std::vector<int>>().swap(m_R);
		for (int l = 0; l < numobj; ++l) {
			for (int k = 0; k < pow(m_grid_div, numobj - 1); ++k)
				if (m_S[l][k].size() != 0)
					std::vector<int>().swap(m_S[l][k]);	//clear m_S	
		}

		/*for (int i = 0; i < offspring.size(); ++i) {
			std::vector<int>().swap(offspring[i].GN());
		}*/

	}

	void PopCDGMOEA::sortl(std::vector<std::pair<double, int>> &v) {	//selection sort
		for (int i = 0; i < v.size(); ++i) {
			int min = i;
			for (int j = i + 1; j < v.size(); ++j) {
				if (v[j].first < v[min].first)
					min = j;
			}
			auto temp = v[i];
			v[i] = v[min];
			v[min] = temp;
		}
	}

	void PopCDGMOEA::sortR(const std::vector<int> &R, std::vector<int> &R_) {
		for (int i = 0; i < R.size(); ++i)
			R_[i] = R[i];
		for (int i = 0; i < R_.size(); ++i) {
			int min = i;
			for (int j = i + 1; j < R_.size(); ++j) {
				if (R_[j] < R_[min])
					min = j;
			}
			auto temp = R_[i];
			R_[i] = R_[min];
			R_[min] = temp;
		}
	}

	void PopCDGMOEA::sortL(std::vector<std::vector<int>> &R_, std::vector<int> &number, Problem *pro) {
		int numobj = pro->numberObjectives();
		for (int l = 0; l < numobj; ++l) {
			for (int i = 0; i < R_.size(); ++i) {
				int min = i;
				for (int j = i + 1; j < R_.size(); ++j) {
					if (l == 0 && R_[j][l] < R_[min][l])
						min = j;
					else if (l == 1 && R_[j][0] == R_[min][0] && R_[j][l] < R_[min][l])
						min = j;
					else if (l == 2 && R_[j][0] == R_[min][0] && R_[j][1] == R_[min][1] && R_[j][l] < R_[min][l])
						min = j;
				}
				if (i != min) {
					auto temp1 = R_[i];
					R_[i] = R_[min];
					R_[min] = temp1;
					int temp2 = number[i];
					number[i] = number[min];
					number[min] = temp2;
				}
			}
		}
	}

	void PopCDGMOEA::selectRandom(int number, std::vector<int> &result, Random *rnd) {
		std::vector<int> candidate;
		for (int i = 0; i < m_individuals.size(); ++i) {
			candidate.push_back(i);
		}
		result.resize(number);
		for (int i = 0; i < number; ++i) {
			int idx = rnd->uniform.nextNonStd(0, (int)candidate.size() - i);
			result[i] = candidate[idx];
			if (idx != candidate.size() - (i + 1)) candidate[idx] = candidate[candidate.size() - (i + 1)];
		}
	}

	void PopCDGMOEA::nondominatedSorting(std::vector<IndCDGMOEA> &offspring, Problem *pro) {
		std::vector<std::vector<Real> *> objs;
		for (auto &i : offspring)
			objs.emplace_back(&i.objective());
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs, rank, pro->optimizeMode());
		for (size_t i = 0; i < offspring.size(); ++i)
			offspring[i].setFitness(rank[i]);
	}

	int PopCDGMOEA::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		int tag = kNormalEval;
		//std::cout << "the " << m_iteration << " generation" << std::endl;
		/*step2:reproduction*/
		tag = evolveMO(pro, alg, rnd);
		/*step3:update of the ideal and nadir points*/
		updateIdealPoint(m_offspring, pro);
		if ((m_iteration % 10) == 0 && m_iteration != 0)
			updateNadirPoint(m_offspring, pro);	//updated nadir point every 50 generations
		/*step4:update of the grid system*/
		evalEens(pro, rnd);
		/*record objective every generation*/
		//record_x();
		//record_x_offspring();
		//record_f();
		//record_f_offspring();

		////下面测试输出结果
		//for (const auto& i : m_individuals) {
		//	std::cout << i->objective()[0] << "    " << i->objective()[1] << std::endl;
		//}
		//std::cout << "* * * * * * * * * * * * * * *" << std::endl;
		//std::cout << "第" << "  " << m_iteration << "  代" << std::endl;
		//std::cout << "* * * * * * * * * * * * * * *" << std::endl;

		m_iteration++;
		return tag;
	}

	int PopCDGMOEA::evolveMO(Problem *pro, Algorithm *alg, Random *rnd) {		//step2:reproduction	
		int tag = kNormalEval;
		//assign NS of Solution in P
		if (m_iteration) {
			for (size_t i = 0; i < m_individuals.size(); ++i) {
				std::vector<size_t>().swap(m_individuals[i]->GN());
			}
			assign_GN(pro);
		}
		int m = 0;
		for (int i = 0; i < m_individuals.size(); ++i) {
			double rand = rnd->uniform.next();
			std::vector<size_t> result(3);
			if (rand < m_delta && m_individuals[i]->GN().size()>3)
				selectInCandidates(3, m_individuals[i]->GN(), result, rnd);
			else
				select(i, 3, result, rnd);
			/*adopt DE operator in MOEA/D*/
			m_offspring[m++] = *m_individuals[i];
			crossMutate(result, m_offspring[m], pro, rnd);
			tag = m_offspring[m].evaluate(pro, alg);
			if (tag != kNormalEval) break;
			m++;
		}
		return tag;
	}

	void PopCDGMOEA::evalEens(Problem *pro, Random *rnd) {	//step5:rank_based selection
		//delete Solution of objective>m_nadir	
		std::vector<int> P_delete, P_reserve;
		for (int i = 0; i < m_offspring.size(); ++i) {
			int flag = 0;
			for (int j = 0; j < pro->numberObjectives(); ++j) {
				if (m_nadir[j] < m_offspring[i].objective()[j])
					flag = 1;
			}
			if (flag)
				P_delete.push_back(i);
			else
				P_reserve.push_back(i);
		}
		if ((P_delete.size() + P_reserve.size()) != m_offspring.size())
			std::cout << "step3 is error" << std::endl;

		if (P_reserve.size() != 0)
			gridConstruct_assignGN_P_reserve(m_offspring, P_reserve, P_reserve.size(), pro);
		else
			//	std::cout << "the " << m_iteration << " generation : " << "reserve = 0" << std::endl;

			if (P_reserve.size() < m_individuals.size()) {	// when Preserve < N
				for (int i = 0; i < P_reserve.size(); ++i) {
					*m_individuals[i] = m_offspring[P_reserve[i]];
				}
				//select (N-P_reserve) Solution from P_delete randomly
				std::vector<int> result(m_individuals.size() - P_reserve.size());
				selectRandom(m_individuals.size() - P_reserve.size(), result, rnd);
				int m = 0;
				for (int i = P_reserve.size(); i < m_individuals.size(); ++i) {
					*m_individuals[i] = m_offspring[P_delete[result[m++]]];
					//m_individuals[i] = m_offspring[P_delete[result[m++]]];
				}
			}
			else {	// when Preserve > N , /*step5:rank_based selection*/ RBS
				assign_S(m_offspring, P_reserve, pro);
				RBS(m_offspring, P_reserve, pro);
			}
	}

	/*void PopCDGMOEA::record_x() {
		int n = global::ms_global->m_problem->variable_size();
		if (m_iteration == 1) {
			std::ofstream out;
			out.open("./result/CDG_MOEA_var.txt");
			if (out && n == 2) {
				out << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "x1" << std::setw(10) << "x2" << std::endl;
				for (int i = 0; i < m_individuals.size(); ++i) {
					for (int j = 0; j < n; ++j) {
						out << std::setw(10) << m_individuals[i]->variable()[j];
					}
					out << std::endl;
				}
			}
			else if (out && n == 3) {
				out << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "x1" << std::setw(10) << "x2" << std::setw(10) << "x3" << std::endl;
				for (int i = 0; i < m_individuals.size(); ++i) {
					for (int j = 0; j < n; ++j) {
						out << std::setw(10) << m_individuals[i]->variable()[j];
					}
					out << std::endl;
				}
			}
			out.close();
		}
		else {
			std::ofstream out;
			out.open("./result/CDG_MOEA_var.txt", std::ios::app);
			if (out && n == 2) {
				out << std::endl << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "x1" << std::setw(10) << "x2" << std::endl;
				for (int i = 0; i < m_individuals.size(); ++i) {
					for (int j = 0; j < n; ++j) {
						out << std::setw(10) << m_individuals[i]->variable()[j];
					}
					out << std::endl;
				}
			}
			else if (out && n == 3) {
				out << std::endl << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "x1" << std::setw(10) << "x2" << std::setw(10) << "x3" << std::endl;
				for (int i = 0; i < m_individuals.size(); ++i) {
					for (int j = 0; j < n; ++j) {
						out << std::setw(10) << m_individuals[i]->variable()[j];
					}
					out << std::endl;
				}
			}
			out.close();
		}

	}

	void PopCDGMOEA::record_x_offspring() {
		int n = global::ms_global->m_problem->variable_size();
		if (m_iteration == 1) {
			std::ofstream out;
			out.open("./result/CDG_MOEA_var_offspring.txt");
			if (out && n == 2) {
				out << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "x1" << std::setw(10) << "x2" << std::endl;
				for (int i = 0; i < m_offspring.size(); ++i) {
					for (int j = 0; j < n; ++j) {
						out << std::setw(10) << m_offspring[i].variable()[j];
					}
					out << std::endl;
				}
			}
			out.close();
		}
		else {
			std::ofstream out;
			out.open("./result/CDG_MOEA_var_offspring.txt", std::ios::app);
			if (out && n == 2) {
				out << std::endl << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "x1" << std::setw(10) << "x2" << std::endl;
				for (int i = 0; i < m_offspring.size(); ++i) {
					for (int j = 0; j < n; ++j) {
						out << std::setw(10) << m_offspring[i].variable()[j];
					}
					out << std::endl;
				}
			}
			out.close();
		}

	}

	void PopCDGMOEA::record_f() {
		int m = global::ms_global->m_problem->objective_size();
		if (m_iteration == 1) {
			std::ofstream out;
			out.open("./result/CDG_MOEA_obj.txt");
			if (out && m == 2) {
				out << "the  " << m_iteration << "  generation" << std::endl;
				out << "f1\t" << "f2" << std::endl;
				for (int i = 0; i < m_individuals.size(); ++i) {
					for (int j = 0; j < m; ++j) {
						out << "\t" << m_individuals[i]->objective()[j];
					}
					out << std::endl;
				}
			}
			else if (out && m == 3) {
				out << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::setw(10) << "f3" << std::endl;
				for (int i = 0; i < m_individuals.size(); ++i) {
					for (int j = 0; j < m; ++j) {
						out << std::setw(10) << m_individuals[i]->objective()[j];
					}
					out << std::endl;
				}
			}
			out.close();
		}
		else {
			std::ofstream out;
			out.open("./result/CDG_MOEA_obj.txt", std::ios::app);
			if (out && m == 2) {
				out << std::endl << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::endl;
				for (int i = 0; i < m_individuals.size(); ++i) {
					for (int j = 0; j < m; ++j) {
						out << std::setw(10) << m_individuals[i]->objective()[j];
					}
					out << std::endl;
				}
			}
			else if (out && m == 3) {
				out << std::endl << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::setw(10) << "f3" << std::endl;
				for (int i = 0; i < m_individuals.size(); ++i) {
					for (int j = 0; j < m; ++j) {
						out << std::setw(10) << m_individuals[i]->objective()[j];
					}
					out << std::endl;
				}
			}
			out.close();
		}
	}

	void PopCDGMOEA::record_f_offspring() {
		int m = global::ms_global->m_problem->objective_size();
		if (m_iteration == 1) {
			std::ofstream out;
			out.open("./result/CDG_MOEA_obj_offspring.txt");
			if (out && m == 2) {
				out << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::endl;
				for (int i = 0; i < m_offspring.size(); ++i) {
					for (int j = 0; j < m; ++j) {
						out << std::setw(10) << m_offspring[i].objective()[j];
					}
					out << std::endl;
				}
			}
			out.close();
		}
		else {
			std::ofstream out;
			out.open("./result/CDG_MOEA_obj_offspring.txt", std::ios::app);
			if (out && m == 2) {
				out << std::endl << "the  " << m_iteration << "  generation" << std::endl;
				out << std::setw(10) << "f1" << std::setw(10) << "f2" << std::endl;
				for (int i = 0; i < m_offspring.size(); ++i) {
					for (int j = 0; j < m; ++j) {
						out << std::setw(10) << m_offspring[i].objective()[j];
					}
					out << std::endl;
				}
			}
			out.close();
		}

	}*/
}
