#include "nsgaii_de.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/environment_selection//selection_methods.h"
#include "../../../../../utility/metricsMOP/IGD.h"
#include "../../../../record/rcr_vec_real.h"
#include <algorithm>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void NSGAII_DE::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		if (m_pop_size < 5)
			throw "Population size of NSGAIIDE should not be smaller than 5.";
		m_mr = v.has("crossover rate") ? v.get<Real>("crossover rate") : 1.0 / CAST_CONOP(m_problem.get())->numberVariables();
		m_meta = v.has("crossover rate") ? v.get<Real>("crossover rate") : 20.;
		m_pop.reset();
	}

	void NSGAII_DE::initPop() {
		m_pop.reset(new PopNSGAII_DE(m_pop_size, m_problem.get()));
		m_pop->initialize(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
		m_pop->setParamMODE(m_mr, m_meta);
		updateHistorySols(*m_pop);
	}

	void NSGAII_DE::run_() {
		initPop();
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
			Population<IndDE> temp_pop;
			for (size_t i = 0; i < m_pop->size(); ++i) {
				temp_pop.append(m_pop->getCombinedPop()[i]);
			}
			updateHistorySols(temp_pop);
			updateHistoryFrontSols(temp_pop,m_problem.get());
			recordMetrics(m_problem.get(), this);

#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void NSGAII_DE::updateHistorySols(Population<IndDE>& pop) {
		for (size_t i = 0; i < pop.size(); ++i) {
			m_historical_sols.emplace_back(std::make_shared<Solution<>>(pop[i]));
		}
	}

	void NSGAII_DE::updateHistoryFrontSols(Population<IndDE>& pop, Problem *pro) {
		//更新历史前沿
		//先对pop排序
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < pop.size(); ++i)
			objs.emplace_back(&pop[i].objective());
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
		for (size_t i = 0; i < pop.size(); ++i)
			pop[i].setFitness(rank[i]);
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < pop.size(); ++i) {
			if (pop[i].fitness() == 0) {
				temp_pop.append(pop[i]);
			}
		}
		if (m_historical_front_sols.empty()) {
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				m_historical_front_sols.emplace_back(std::make_shared<Solution<>>(temp_pop[i]));
			}
		}
		else {
			//先看每个解是否被当前前沿支配
			std::vector<size_t> temp_his_sols(m_historical_front_sols.size(), 1);
			std::vector<size_t> temp_pop_sols(temp_pop.size(), 0);
			for (size_t j = 0; j < temp_pop.size(); ++j) {
				for (size_t i = 0; i < m_historical_front_sols.size(); ++i) {
					if (temp_his_sols[i] == 1) {
						Dominance dominanceship = temp_pop[j].compare(*m_historical_front_sols[i], CAST_CONOP(pro)->optimizeMode());
						if (dominanceship == Dominance::kNonDominated) {
							if (i == m_historical_front_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else if (dominanceship == Dominance::kDominant) {
							temp_his_sols[i] = 0;
							if (i == m_historical_front_sols.size() - 1) {
								temp_pop_sols[j] = 1;
							}
						}
						else {
							break;
						}
					}
				}
			}
			auto temp = m_historical_front_sols;
			m_historical_front_sols.clear();
			for (size_t i = 0; i < temp.size(); ++i) {
				if (temp_his_sols[i] == 1) {
					m_historical_front_sols.emplace_back(temp[i]);
				}
			}
			for (size_t i = 0; i < temp_pop.size(); ++i) {
				if (temp_pop_sols[i] == 1) {
					m_historical_front_sols.emplace_back(std::make_shared<Solution<>>(temp_pop[i]));
				}
			}
		}
	}

	void NSGAII_DE::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

	void NSGAII_DE::recordMetrics(Problem *pro, Algorithm *alg) {
		/************************************/
		/*            性能指标计算          */
		/************************************/
		//使用archive计算性能指标
		//排序
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < m_pop->size(); ++i)
			objs.emplace_back(&m_pop->at(i).objective());
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
		for (size_t i = 0; i < m_pop->size(); ++i)
			m_pop->at(i).setFitness(rank[i]);
		Population<IndDE> temp_pop;
		for (size_t i = 0; i < m_pop->size(); ++i) {
			if (m_pop->at(i).fitness() == 0) {
				temp_pop.append(m_pop->at(i));
			}
		}
		std::vector<std::vector<Real>> ref_objs;
		for (size_t i = 0; i < m_problem->optimaBase()->numberObjectives(); ++i) {
			ref_objs.push_back(m_problem->optimaBase()->objective(i));
		}
		std::vector<std::vector<Real>> pop_objs;
		for (size_t i = 0; i < m_pop->size(); ++i) {
			pop_objs.push_back(m_pop->at(i).objective());
		}
		Real temp_IGD = IGD(ref_objs, pop_objs);
		//Real temp_IGD = CAST_CONOP(pro)->optima()->invertGenDist(temp_pop);
		getIGD().push_back(temp_IGD);

		std::cout << "累积前沿解个数" << m_historical_front_sols.size() << std::endl;

		std::cout << alg->evaluations() << "  " << temp_IGD << std::endl;
		//record();//store metrics data
	}

#ifdef OFEC_DEMO
	void NSGAII_DE::updateBuffer() {
		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(1);
			for (size_t i = 0; i < m_pop->size(); ++i) {
				m_solution[0].push_back(&m_pop->at(i));
			}
			m_off_solution.clear();
			m_off_solution.resize(1);
			for (size_t i = 0; i < m_pop->size(); ++i) {
				m_off_solution[0].push_back(&m_pop->getCombinedPop()[i]);
			}
			m_his_solution.clear();
			auto& his_sols = getHisSols();
			for (size_t i = 0; i < his_sols.size(); ++i) {
				m_his_solution.push_back(his_sols[i].get());
			}
			ofec_demo::g_buffer->appendAlgBuffer(this);
		}
	}
#endif

	PopNSGAII_DE::PopNSGAII_DE(size_t size_pop, Problem *pro) :
		PopMODE(size_pop, pro),m_pop_combined(2 * size_pop, pro, CAST_CONOP(pro)->numberVariables()) {}

	int PopNSGAII_DE::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		int tag = kNormalEval;
		for (int i = 0; i < m_individuals.size(); i++) {
			std::vector<size_t> p(3);
			p[0] = tournamentSelection(pro, rnd);
			while (1) { p[1] = tournamentSelection(pro, rnd);  	if (p[1] != p[0]) break; }
			while (1) { p[2] = tournamentSelection(pro, rnd);  	if (p[2] != p[0] && p[2] != p[1]) break; }
			crossMutate(p, m_pop_combined[i], pro, rnd);
			tag = m_pop_combined[i].evaluate(pro, alg);
			if (tag != kNormalEval) break;
		}
		for (int i = 0; i < m_individuals.size(); i++) {
			m_pop_combined[m_individuals.size()+i] = *m_individuals[i];
		}
		//survivorSelection(*this, m_pop_combined);
		envirSelectionByNSGAII(*this,m_pop_combined,pro->optimizeMode());
		m_iteration++;
		return tag;
	}
}