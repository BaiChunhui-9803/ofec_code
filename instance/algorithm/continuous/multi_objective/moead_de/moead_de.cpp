#include "moead_de.h"
#include "../../../../record/rcr_vec_real.h"
#include "../../../../../utility/metricsMOP/IGD.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void MOEAD_DE::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		m_pop.reset();
	}

	void MOEAD_DE::record() {
		/*std::vector<Real> entry;
		entry.push_back(m_evaluations);
		Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(IGD);
		dynamic_cast<RecordVectorReal*>(m_record.get())->record(this, entry);*/

		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void MOEAD_DE::updateBuffer() {
		/*m_solution.clear();
		m_solution.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i)
			m_solution[0].push_back(&m_pop->at(i));
		ofec_demo::g_buffer->appendAlgBuffer(this);*/

		if (ofec_demo::g_buffer->algorithm().get() == this) {
			m_solution.clear();
			m_solution.resize(1);
			for (size_t i = 0; i < m_pop->size(); ++i) {
				m_solution[0].push_back(&m_pop->at(i));
			}
			m_off_solution.clear();
			m_off_solution.resize(1);
			for (size_t i = 0; i < m_pop->size(); ++i) {
				m_off_solution[0].push_back(&m_pop->getOffPop()[i]);
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

	void MOEAD_DE::initPop() {
		auto size_var = CAST_CONOP(m_problem.get())->numberVariables();
		auto size_obj = CAST_CONOP(m_problem.get())->numberObjectives();
		m_pop.reset(new PopMOEAD_DE(m_pop_size, m_problem.get()));
		m_pop->initialize_(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
		updatePopRange();
		updateHistorySols(*m_pop);
	}

	void MOEAD_DE::updatePopRange() {
		m_pop->updatePopRange();
	}

	void MOEAD_DE::run_() {
		initPop();
#ifdef OFEC_DEMO
		updateBuffer();
#endif
		while (!terminating()) {
			m_pop->evolve(m_problem.get(),this,m_random.get());

			/*Population<IndDE> new_pop;
			for (size_t i = 0; i < m_pop->size(); ++i) {
				new_pop.append(m_pop->getOffPop()[i]);
			}*/
			updateHistorySols(m_pop->getOffPop());
			updateHistoryFrontSols(m_pop->getOffPop(), m_problem.get());
			recordMetrics(m_problem.get(), this);
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void MOEAD_DE::recordMetrics(Problem *pro, Algorithm *alg) {
		/************************************/
		/*            性能指标计算          */
		/************************************/
		//使用archive计算性能指标
		Population<IndDE> temp_pop;
		//排序
		std::vector<std::vector<Real>*> objs;
		for (size_t i = 0; i < m_pop->size(); ++i)
			objs.emplace_back(&m_pop->at(i).objective());
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs, rank, CAST_CONOP(pro)->optimizeMode());
		for (size_t i = 0; i < m_pop->size(); ++i)
			m_pop->at(i).setFitness(rank[i]);
		for (size_t i = 0; i < m_pop->size(); ++i) {
			if (m_pop->at(i).fitness() == 0.) {
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

	void MOEAD_DE::updateHistorySols(Population<IndDE>& pop) {
		for (size_t i = 0; i < pop.size(); ++i) {
			m_historical_sols.emplace_back(std::make_shared<Solution<>>(pop[i]));
		}
	}

	void MOEAD_DE::updateHistoryFrontSols(Population<IndDE>& pop, Problem *pro) {
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

	std::vector<Vector> MOEAD_DE::getWeigh() {
		auto ref_vector = m_pop->getWeigh();
		return ref_vector;
	}

	PopMOEAD_DE::PopMOEAD_DE(size_t size_pop, Problem *pro) :
		PopMODE<>(size_pop,pro), MOEAD<IndDE>(CAST_CONOP(pro)->numberObjectives()), m_off_pop(size_pop, pro, CAST_CONOP(pro)->numberVariables()) { }

	void PopMOEAD_DE::initialize_(Problem *pro, Random *rnd) {
		Population<IndDE>::initialize(pro,rnd);
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		MOEAD<IndDE>::initialize(this->m_individuals, pro);
		m_pop_range.resize(num_obj);
		for (size_t i = 0; i < num_obj; ++i) {
			m_pop_range[i].first = INT_MAX;
			m_pop_range[i].second = INT_MIN;
		}
	}

	int PopMOEAD_DE::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		m_interactive_sol_pair.clear();
		int tag = kNormalEval;
		std::vector<int> perm(this->m_individuals.size());
		for (int i(0); i < perm.size(); ++i) {
			perm[i] = i;
		}
		rnd->uniform.shuffle(perm.begin(), perm.end());
		for (int i = 0; i < this->m_individuals.size(); i++) {
			int n = perm[i];
			// or int n = i;
			int type;
			double r = rnd->uniform.next();
			// mating selection based on probability
			if (r < m_realb)
				type = 1;   // neighborhood
			else
				type = 2;   // whole population

			// select the indexes of mating parents
			std::vector<int> p;
			matingSelection(p, n, 2, type, this->m_individuals.size(), rnd);  // neighborhood selection

			// produce a child Solution
			ofec::IndDE child = *this->m_individuals[0];
			std::vector<size_t> index(3);
			index[0] = n; index[1] = p[0]; index[2] = p[1];
			this->crossMutate(index, child, pro, rnd);
			tag = child.evaluate(pro, alg);
			m_off_pop[i] = child;

			//记录交互解
			Solution<> ind1(*this->m_individuals.at(index[0]));
			Solution<> ind2(*this->m_individuals.at(index[1]));
			Solution<> ind3(*this->m_individuals.at(index[2]));
			Solution<> ind4(child);
			std::vector<std::shared_ptr<Solution<>>> temp_pair;
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind1));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind2));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind3));
			temp_pair.emplace_back(std::make_shared<Solution<>>(ind4));
			m_interactive_sol_pair.emplace_back(temp_pair);

			//if the range of objs less than 0 
			for (size_t j = 0; j < child.objective().size(); ++j) {
				if (child.objective()[j] < 0) {
					size_t a = 1;
				}
			}
			if (tag != kNormalEval) break;
			// update the reference points and other TypeIndivs in the neighborhood or the whole population
			updateReference(child, pro);
			updateProblem(this->m_individuals, child, n, type, pro, rnd);

			
		}
		updatePopRange();
		m_iteration++;
		return tag;
	}

	void PopMOEAD_DE::updatePopRange() {
		for (size_t i = 0; i < m_pop_range.size(); ++i) {
			for (size_t j = 0; j < this->m_individuals.size(); ++j) {
				if (m_pop_range[i].first > this->m_individuals[j]->objective()[i]) {
					m_pop_range[i].first = this->m_individuals[j]->objective()[i];
				}
				if (m_pop_range[i].second < this->m_individuals[j]->objective()[i]) {
					m_pop_range[i].second = this->m_individuals[j]->objective()[i];
				}
			}
		}
	}
}
