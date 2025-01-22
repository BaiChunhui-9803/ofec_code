#include "nsgaii_sbx.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../record/rcr_vec_real.h"
#include "../../../../../utility/environment_selection/selection_methods.h"
#include "../../../../../utility/metricsMOP/IGD.h"
#include <fstream>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void NSGAII_SBX::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		if (m_pop_size % 2)
			throw MyExcept("Population size of NSGAII should be even.");

		m_cr = v.has("crossover rate") ? v.get<Real>("crossover rate") : 0.9;
		m_mr = v.has("mutation rate") ? v.get<Real>("mutation rate") : 1.0 / CAST_CONOP(m_problem.get())->numberVariables();
		m_ceta = v.has("crossover eta") ? v.get<Real>("crossover eta") : 20.;
		m_meta = v.has("mutation eta") ? v.get<Real>("mutation eta") : 20.;
		

		//std::string pro_name = std::get<std::string>(v.at("problem name"));
		//initPop();
		////为种群赋新值
		////读取目标点，构造解，为种群赋值
		//std::vector<std::vector<Real>> obj_sols;
		//std::ifstream infile;
		//std::stringstream os;
		//os << g_working_dir << "/instance/problem/continuous/multi_objective/glt/data/" << pro_name << ".dat";
		//infile.open(os.str());
		//if (!infile) {
		//	throw MyExcept("open PF file of GLT problem is fail");
		//}
		//std::string str;
		//size_t line = 0;
		//while (getline(infile, str))
		//	++line;
		//obj_sols.resize(line);
		//infile.close();
		//infile.clear();
		//infile.open(os.str());
		//size_t num_obj = std::get<int>(v.at("number of objectives"));
		//for (size_t i = 0; i < line; i++) {
		//	std::vector<Real> temp_obj(num_obj);
		//	for (size_t j = 0; j < num_obj; j++) {
		//		infile >> temp_obj[j];
		//	}
		//	obj_sols[i]=temp_obj;
		//}
		//infile.close();

		////为组合的个体赋采样点
		//for (size_t i = 0; i < m_pop->size(); ++i) {
		//	m_pop->getCombinPop()[i].objective() = obj_sols[i];
		//}

	}

	void NSGAII_SBX::initPop() {
		m_pop.reset(new PopNSGAII_SBX(m_pop_size, m_problem.get()));
		m_pop->setRate(m_cr, m_mr);
		m_pop->setEta(m_ceta, m_meta);
		m_pop->initialize(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
		updateHistorySols(*m_pop);
	}

 	void NSGAII_SBX::run_() {
		initPop();
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
			Population<Solution<>> temp_pop;
			for (size_t i = 0; i < m_pop->size(); ++i) {
				temp_pop.append(m_pop->getCombinedPop()[i]);
			}
			updateHistorySols(temp_pop);
			//updateHistoryFrontSols(temp_pop, m_problem.get());
			recordMetrics(m_problem.get(), this);
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}

#ifdef OFEC_DEMO
	void NSGAII_SBX::updateBuffer() {
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

	void NSGAII_SBX::updateHistorySols(Population<Solution<>>& pop) {
		for (size_t i = 0; i < pop.size(); ++i) {
			m_historical_sols.emplace_back(std::make_shared<Solution<>>(pop[i]));
		}
	}

	void NSGAII_SBX::updateHistoryFrontSols(Population<Solution<>>& pop, Problem *pro) {
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

	void NSGAII_SBX::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		//Real IGD = m_problem->optima().invertGenDist(*m_pop);
		entry.push_back(getIGD().back());
		dynamic_cast<RecordVecRealMOEA*>(m_record.get())->record(this, entry);
	}

	void NSGAII_SBX::recordMetrics(Problem *pro, Algorithm *alg) {
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
		Real temp_IGD = IGD(ref_objs, pop_objs);
		//Real temp_IGD = CAST_CONOP(pro)->optima()->invertGenDist(temp_pop);
		getIGD().push_back(temp_IGD);

		/*std::cout << "累积前沿解个数" << m_historical_front_sols.size() << std::endl;

		std::cout << alg->evaluations() << "  " << temp_IGD << std::endl;*/
		//record();//store metrics data
	}



	/* ============== population ============== */

	PopNSGAII_SBX::PopNSGAII_SBX(size_t size_pop, Problem *pro) :
		PopSBX<>(size_pop, pro),m_pop_combined(size_pop * 2, pro, CAST_CONOP(pro)->numberVariables()) {}

	int PopNSGAII_SBX::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		reproduction(m_pop_combined, pro, rnd);
		int tag = kNormalEval;
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			tag = m_pop_combined[i].evaluate(pro, alg);
			if (tag != kNormalEval) return tag;
		}
		//survivorSelection(getPop(), m_pop_combined);
		envirSelectionByNSGAII(*this, m_pop_combined, pro->optimizeMode());
		m_iteration++;
		return tag;
	}
}