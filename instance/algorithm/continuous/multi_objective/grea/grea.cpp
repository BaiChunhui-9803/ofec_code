#include "grea.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include <algorithm>
#include "../../../../record/rcr_vec_real.h"
#include "../../../../../utility/environment_selection/selection_methods.h"
#include "../../../../../utility/metricsMOP/IGD.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	/**********************GrEA************************/		
	void GrEA::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		m_pop.reset();
	}

	void GrEA::initPop() {
		auto size_var = CAST_CONOP(m_problem.get())->numberVariables();
		auto size_obj = CAST_CONOP(m_problem.get())->numberObjectives();
		m_pop.reset(new PopGrEA(m_pop_size, m_problem.get(), size_obj));
		m_pop->initialize(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
	}

	void GrEA::run_() {
		initPop();
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void GrEA::record() {
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
		//Real IGD = dynamic_cast<Optima<>&>(*const_cast<OptimaBase*>(m_problem->optimaBase())).invertGenDist(*m_pop);
		entry.push_back(temp_IGD);
		dynamic_cast<RecordVectorReal*>(m_record.get())->addEntry(this, entry);
	}

#ifdef OFEC_DEMO
	void GrEA::updateBuffer() {
		m_solution.clear();
		m_solution.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i)
			m_solution[0].push_back(&m_pop->at(i));
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif

	/******************************population*************************/
	PopGrEA::PopGrEA(size_t size_pop, Problem *pro, size_t size_obj) : 
		PopSBX<IndGrEA>(size_pop, pro),
		m_offspring(2 * size_pop, pro, CAST_CONOP(pro)->numberVariables()),
		m_ind_min_max(size_obj),
		m_grid_min_max(size_obj),
		m_grid_distance(size_obj) {}

	void PopGrEA::initialize(Problem *pro, Random *rnd) {
		Population<IndGrEA>::initialize(pro,rnd);
		setRate(0.9, 1.0 / CAST_CONOP(pro)->numberVariables());
		setEta(20, 20);
		gridConstruct(pro);
		assignGR_GCPD(pro);
		assignGCD(pro);
	}

	void PopGrEA::gridConstruct(Problem *pro) {
		int num_obj = pro->numberObjectives();
		for (int i = 0; i < num_obj; ++i) {
			m_ind_min_max[i].second = -1 * 1.0e14;
			m_ind_min_max[i].first = 1.0e14;
			for (int j = 0; j < m_individuals.size(); ++j) {
				if (m_individuals[j]->objective(i) > m_ind_min_max[i].second)
					m_ind_min_max[i].second = m_individuals[j]->objective(i);
				if (m_individuals[j]->objective(i) < m_ind_min_max[i].first)
					m_ind_min_max[i].first = m_individuals[j]->objective(i);
			}
			m_grid_distance[i] = (m_ind_min_max[i].second - m_ind_min_max[i].first) * (m_grid_div + 1) / (m_grid_div * m_grid_div);
			m_grid_min_max[i].second = m_ind_min_max[i].second + (m_grid_distance[i] * m_grid_div - (m_ind_min_max[i].second - m_ind_min_max[i].first)) / 2;
			m_grid_min_max[i].first = m_ind_min_max[i].first - (m_grid_distance[i] * m_grid_div - (m_ind_min_max[i].second - m_ind_min_max[i].first)) / 2;
		}
	}

	void PopGrEA::gridConstructFi(Population<IndGrEA>& offspring, std::vector<int>& Fi, int size,Problem *pro) {
		int num_obj=pro->numberObjectives();
		for (int i = 0; i < num_obj; ++i) {
			m_ind_min_max[i].second = -1 * 1.0e14;
			m_ind_min_max[i].first = 1.0e14;
			for (int j = 0; j < size; ++j) {
				if (offspring[Fi[j]].objective()[i] > m_ind_min_max[i].second)
					m_ind_min_max[i].second = offspring[Fi[j]].objective()[i];
				if (offspring[Fi[j]].objective()[i] < m_ind_min_max[i].first)
					m_ind_min_max[i].first = offspring[Fi[j]].objective()[i];
			}
			m_grid_distance[i] = (m_ind_min_max[i].second - m_ind_min_max[i].first) * (m_grid_div + 1) / (m_grid_div * m_grid_div);
			m_grid_min_max[i].second = m_ind_min_max[i].second + (m_grid_distance[i] * m_grid_div - (m_ind_min_max[i].second - m_ind_min_max[i].first)) / 2;
			m_grid_min_max[i].first = m_ind_min_max[i].first - (m_grid_distance[i] * m_grid_div - (m_ind_min_max[i].second - m_ind_min_max[i].first)) / 2;
		}
	}

	void PopGrEA::assignGR_GCPD(Problem *pro) {
		int num_obj=pro->numberObjectives();
		int flag;
		double value;
		for (int i = 0; i < m_individuals.size(); i++) {
			flag = 0;
			value = 0;
			for (int j = 0; j < num_obj; j++) {
				m_individuals[i]->Gk(j) = (int)floor((m_individuals[i]->objective(j) - m_grid_min_max[j].first) / m_grid_distance[j]);
				flag += m_individuals[i]->Gk(j);
				value += pow((m_individuals[i]->objective(j) - (m_grid_min_max[j].first + m_individuals[i]->Gk(j) * m_grid_distance[j])) / m_grid_distance[j], 2.0);
			}
			m_individuals[i]->GR() = flag;
			m_individuals[i]->GCPD() = sqrt(value);
		}
	}

	void PopGrEA::assignGR_GCPD_Fi(Population<IndGrEA>& offspring, std::vector<int>& Fi, int size,Problem *pro) {
		int num_obj = pro->numberObjectives();
		int flag;
		double value;
		for (int i = 0; i < size; i++) {
			flag = 0;
			value = 0;
			for (int j = 0; j < num_obj; j++) {
				offspring[Fi[i]].Gk(j) = (int)floor((offspring[Fi[i]].objective()[j] - m_grid_min_max[j].first) / m_grid_distance[j]);
				flag += offspring[Fi[i]].Gk(j);
				value += pow((offspring[Fi[i]].objective()[j] - (m_grid_min_max[j].first + offspring[Fi[i]].Gk(j) * m_grid_distance[j])) / m_grid_distance[j], 2.0);
			}
			offspring[Fi[i]].GR() = flag;
			offspring[Fi[i]].GCPD() = sqrt(value);
		}
	}

	void PopGrEA::assignGCD(Problem *pro) {
		int flag;
		int num_obj=pro->numberObjectives();
		for (int i = 0; i < m_individuals.size(); i++) {
			for (int j = i + 1; j < m_individuals.size(); j++) {
				flag = 0;
				for (int k = 0; k < num_obj; k++)
					flag += abs(m_individuals[i]->Gk(k) - m_individuals[j]->Gk(k));
				if (flag < num_obj) {
					m_individuals[i]->GCD() += num_obj - flag;
					m_individuals[j]->GCD() += num_obj - flag;
				}
			}
		}
	}

	int PopGrEA::checkDominanceGrid(IndGrEA & a, IndGrEA & b, Problem *pro) {
		int num_obj = pro->numberObjectives();
		int flag1 = 0, flag2 = 0;
		for (int i = 0; i < num_obj; ++i) {
			if (a.Gk(i) < b.Gk(i))
				flag1 = 1;
			else if (a.Gk(i) > b.Gk(i))
				flag2 = 1;
		}
		if (flag1 == 1 && flag2 == 0)
			return 1;
		else if (flag1 == 0 && flag2 == 1)
			return -1;
		else
			return 0;
	}

	void PopGrEA::evalEens(Problem *pro, Random *rnd) {      //environment selection
		int num_obj = pro->numberObjectives();
		int pops = 0;  //indicate parent population size be 0
		int size = m_offspring.size();
		int rank = 0;
		while (true) {
			int count = 0;
			for (size_t i = 0; i < size; i++)
				if (m_offspring[i].fitness() == rank)
					count++;
			int size2 = pops + count;
			if (size2 > m_individuals.size()) {
				break;
			}
			for (size_t i = 0; i < size; i++)
				if (m_offspring[i].fitness() == rank) {
					*(m_individuals[pops]) = m_offspring[i];
					++pops;
				}
			rank++;
			if (pops >= m_individuals.size()) break;
		}
		if (pops < m_individuals.size())
/*				return;
		else*/ {
			//get all Solution of in Fi
			std::vector<int> Fi;
			int Fi_size = 0;
			for (size_t i = 0; i < size; ++i) {
				if (m_offspring[i].fitness() == rank) {
					Fi.push_back(i);
					Fi_size++;
				}
			}
			//grid construct of Fi
			gridConstructFi(m_offspring, Fi, Fi_size,pro);
			//assign fitness of GR and GCPD
			assignGR_GCPD_Fi(m_offspring, Fi, Fi_size,pro);
			for (int j = 0; j < Fi_size; ++j) {
				m_offspring[Fi[j]].GCD() = 0;
			}

			while (pops < m_individuals.size()) {
				//******************find best solution (q) from, algorithm 7
				int q = 0;
				for (int k = 1; k < Fi_size; ++k) {
					if (m_offspring[Fi[k]].GR() < m_offspring[Fi[q]].GR())
						q = k;
					else if (m_offspring[Fi[k]].GR() == m_offspring[Fi[q]].GR()) {
						if (m_offspring[Fi[k]].GCD() < m_offspring[Fi[q]].GCD())
							q = k;
						else if (m_offspring[Fi[k]].GCD() == m_offspring[Fi[q]].GCD())
							if (m_offspring[Fi[k]].GCPD() < m_offspring[Fi[q]].GCPD())
								q = k;
					}
				}
				*(m_individuals[pops++]) = m_offspring[Fi[q]];
				int Fi_q = Fi[q];
				Fi.erase(Fi.begin() + q);
				Fi_size--;
				//******************compute GCD of Fi, algorithm 6
				int flag;
				for (int m = 0; m < Fi_size; ++m) {
					for (int n = m + 1; n < Fi_size; ++n) {
						flag = 0;
						for (int s = 0; s < num_obj; ++s)
							flag += abs(m_offspring[Fi[m]].Gk(s) - m_offspring[Fi[n]].Gk(s));
						if (flag < num_obj) {
							m_offspring[Fi[m]].GCD() = num_obj - flag;
							m_offspring[Fi[n]].GCD() = num_obj - flag;
						}
					}
				}
				//******************adjust GR of individua in Fi, aigorithm 3
				//equal, dominated, non dominate, non equal
				int GD_q_p;
				std::vector<int> PD(Fi_size, 0);
				for (int m = 0; m < Fi_size; ++m) {
					GD_q_p = 0;
					for (int n = 0; n < num_obj; ++n) {
						GD_q_p += abs(m_offspring[Fi_q].Gk(n) - m_offspring[Fi[m]].Gk(n));
					}
					if (GD_q_p == 0)
						m_offspring[Fi[m]].GR() = num_obj + 2;
					else if (checkDominanceGrid(m_offspring[Fi_q], m_offspring[Fi[m]],pro) == 1)
						m_offspring[Fi[m]].GR() = num_obj;
					else {
						if (GD_q_p < num_obj) {
							if (PD[m] < (num_obj - GD_q_p)) {
								PD[m] = num_obj - GD_q_p;
								for (int n = 0; n < Fi_size; ++n)
									if ((m != n) && (checkDominanceGrid(m_offspring[Fi[m]], m_offspring[Fi[n]], pro) == 1) && (PD[n] < PD[m]))
										PD[n] = PD[m];
							}
						}
						m_offspring[Fi[m]].GR() = m_offspring[Fi_q].GR() + PD[m];
					}
				}
			}
		}
	}

	int PopGrEA::evolveMO(Problem *pro, Algorithm *alg, Random *rnd) {
		if (m_individuals.size() % 2 != 0)
			throw "population size should be even @NSGAII_SBXRealMu::evolveMO()";
		int tag = kNormalEval;
		int m = 0;
		for (size_t n = 0; n < m_individuals.size(); n += 2) {
			std::vector<int> p(2);
			p[0] = tournamentSelection(pro, rnd);
			do { p[1] = tournamentSelection(pro, rnd); } while (p[1] == p[0]);
			crossover(p[0], p[1], m_offspring[m], m_offspring[m + 1], pro, rnd);
			mutate(m_offspring[m], pro, rnd);
			mutate(m_offspring[m + 1], pro, rnd);
			tag = m_offspring[m].evaluate(pro, alg);
			if (tag != kNormalEval) break;
			tag = m_offspring[m + 1].evaluate(pro, alg);
			if (tag != kNormalEval) break;
			m += 2;
			m_offspring[m++] = *m_individuals[n];
			m_offspring[m++] = *m_individuals[n + 1];
		}
		return tag;
	}

	int PopGrEA::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		int tag = kNormalEval;
		tag = evolveMO(pro,alg, rnd);
		//nondominatedSorting(m_offspring);
		popNondominatedSorting(m_offspring, pro->optimizeMode());
		evalEens(pro, rnd);
		m_iteration++;
		return tag;
	}
}