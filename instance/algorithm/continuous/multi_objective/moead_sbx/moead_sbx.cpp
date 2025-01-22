#include "moead_sbx.h"
#include "../../../../record/rcr_vec_real.h"
#include "../../../../../utility/metricsMOP/IGD.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void MOEAD_SBX::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		m_cr = v.has("crossover rate") ? v.get<Real>("crossover rate") : 0.9;
		m_mr = v.has("mutation rate") ? v.get<Real>("mutation rate") : 1.0 / CAST_CONOP(m_problem.get())->numberVariables();
		m_ceta = v.has("crossover eta") ? v.get<Real>("crossover eta") : 20.;
		m_meta = v.has("mutation eta") ? v.get<Real>("mutation eta") : 20.;
	}

	void MOEAD_SBX::initPop() {
		m_pop.reset(new PopMOEAD_SBX(m_pop_size, m_problem.get()));
		m_pop->setRate(m_cr, m_mr);
		m_pop->setEta(m_ceta, m_meta);
		m_pop->initialize_(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
		updatePopRange();
		updateHistorySols(*m_pop);
	}

	void MOEAD_SBX::updatePopRange() {
		m_pop->updatePopRange();
	}

	void MOEAD_SBX::run_() {
		initPop();
#ifdef OFEC_DEMO
		updateBuffer();
#endif
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
			updateHistorySols(m_pop->getOffPop());
			updateHistoryFrontSols(m_pop->getOffPop(), m_problem.get());
			recordMetrics(m_problem.get(), this);

#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

	void MOEAD_SBX::record() {
		//std::vector<Real> entry;
		//entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		//entry.push_back(IGD);
		//dynamic_cast<RecordVectorReal*>(m_record.get())->record(this, entry);

		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

	void MOEAD_SBX::recordMetrics(Problem *pro, Algorithm *alg) {
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
		Population<Solution<>> temp_pop;
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
		Real temp_IGD = IGD(ref_objs,pop_objs);
		getIGD().push_back(temp_IGD);

		std::cout << "累积前沿解个数" << m_historical_front_sols.size() << std::endl;


		std::cout << alg->evaluations() << "  " << temp_IGD << std::endl;
		//record();//store metrics data
	}

	void MOEAD_SBX::updateHistorySols(Population<Solution<>>& pop) {
		for (size_t i = 0; i < pop.size(); ++i) {
			m_historical_sols.emplace_back(std::make_shared<Solution<>>(pop[i]));
		}
	}

	void MOEAD_SBX::updateHistoryFrontSols(Population<Solution<>>& pop, Problem *pro) {
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

#ifdef OFEC_DEMO
	void MOEAD_SBX::updateBuffer() {
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

	std::vector<Vector> MOEAD_SBX::getWeigh() {
		auto ref_vector = m_pop->getWeigh();
		return ref_vector;
	}

	PopMOEAD_SBX::PopMOEAD_SBX(size_t size_pop, Problem *pro) :
		PopSBX<>(size_pop, pro), MOEAD<Solution<>>(CAST_CONOP(pro)->numberObjectives()), m_off_pop(size_pop, pro, CAST_CONOP(pro)->numberVariables()) {}

	void PopMOEAD_SBX::initialize_(Problem *pro, Random *rnd) {
		Population<Solution<>>::initialize(pro, rnd);
		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		MOEAD<Solution<>>::initialize(this->m_individuals, pro);
		m_pop_range.resize(num_obj);
		for (size_t i = 0; i < num_obj; ++i) {
			m_pop_range[i].first = INT_MAX;
			m_pop_range[i].second = INT_MIN;
		}
	}

	int PopMOEAD_SBX::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		int tag = kNormalEval;
		std::vector<int> perm(this->m_individuals.size());
		for (int i(0); i < perm.size(); ++i) {
			perm[i] = i;
		}
		rnd->uniform.shuffle(perm.begin(), perm.end());
		for (int i = 0; i < this->m_individuals.size(); i += 2) {
			int n = perm[i];
			// or int n = i;
			int type = 1; //only select parent Solution from neighborhood
			double rand = rnd->uniform.next();

			// select the indexes of mating parents
			std::vector<int> p;
			matingSelection(p, n, 1, type, this->m_individuals.size(), rnd);  // neighborhood selection

			//produce a child Solution
			Solution<> child1 = *this->m_individuals[0];
			Solution<> child2 = *this->m_individuals[0];
			std::vector<size_t> index(2);
			index[0] = n; index[1] = p[0];
			this->crossover(index[0], index[1], child1, child2, pro, rnd);
			mutate(child1, pro, rnd);
			mutate(child2, pro, rnd);
			tag = child1.evaluate(pro, alg);
			if (tag != kNormalEval) break;
			tag = child2.evaluate(pro, alg);
			if (tag != kNormalEval) break;
			m_off_pop[i] = child1;
			m_off_pop[i+1] = child2;
			// update the reference points and other TypeIndivs in the neighborhood or the whole population
			updateReference(child1, pro);
			updateReference(child2, pro);
			updateProblem(this->m_individuals, child1, n, type, pro, rnd);
			updateProblem(this->m_individuals, child2, n, type, pro, rnd);
		}
		updatePopRange();
		m_iteration++;
		return tag;
	}

	void PopMOEAD_SBX::updatePopRange() {
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
