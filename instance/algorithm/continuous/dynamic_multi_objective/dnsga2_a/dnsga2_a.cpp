#include "dnsga2_a.h"
#include "../../../../problem/continuous/multi_objective/dynamic/DMOPs.h"
#include "../../../../../utility/metricsMOP/IGD.h"
#include "../../../../../utility/environment_selection/selection_methods.h"
#include "../../../../record/rcr_vec_real.h"
#include <algorithm>
#include<fstream>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void DNSGAII_A::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		m_cr = v.has("crossover rate") ? v.get<Real>("crossover rate") : 0.9;
		m_mr = v.has("mutation rate") ? v.get<Real>("mutation rate") : 1.0 / CAST_CONOP(m_problem.get())->numberVariables();
		m_ceta = v.has("crossover eta") ? v.get<Real>("crossover eta") : 20.;
		m_meta = v.has("mutation eta") ? v.get<Real>("mutation eta") : 20.;
		
		m_pop.reset();
	}

	void DNSGAII_A::initPop() {
		auto size_var = CAST_CONOP(m_problem.get())->numberVariables();
		auto size_obj = CAST_CONOP(m_problem.get())->numberObjectives();
		m_pop.reset(new DNSGAII_A_pop(m_pop_size, m_problem.get(), size_var, size_obj, CAST_CONOP(m_problem.get())->optimizeMode()));
		m_pop->initialize_(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
	}

	void DNSGAII_A::run_() {
		initPop();
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			m_solution.clear();
			m_solution.resize(1);
			for (size_t i = 0; i < m_pop->size(); ++i)
				m_solution[0].emplace_back(&m_pop->at(i));
			ofec_demo::g_buffer->appendAlgBuffer(this);
#endif
		}
	}

	void DNSGAII_A::record() {
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

	DNSGAII_A_pop::DNSGAII_A_pop(size_t size_pop, Problem* pro, size_t size_var, size_t size_obj, const std::vector<OptimizeMode>& opt_mode) :
		PopSBX<>(size_pop, pro){}

	void DNSGAII_A_pop::initialize_(Problem* pro, Random* rnd) {
		Population<Solution<>>::initialize(pro, rnd);
		setDetectRate(0.1);
		setReplaceRate(0.2);
		for (auto& i : this->m_individuals) {
			m_offspring.emplace_back(*i);
		}
		for (auto& i : this->m_individuals) {
			m_offspring.emplace_back(*i);
		}
		
	}

	int DNSGAII_A_pop::evolve(Problem* pro, Algorithm* alg, Random* rnd) {
		if (this->m_individuals.size() % 2 != 0)
			throw "population size should be even @NSGAII_SBXRealMu::evolve()";
		int tag = kNormalEval;
		for (size_t i = 0; i < this->m_individuals.size(); i += 2) {
			std::vector<size_t> p(2);
			p[0] = tournamentSelection(pro, rnd);
			do { p[1] = tournamentSelection(pro, rnd); } while (p[1] == p[0]);
			crossover(p[0], p[1], m_offspring[i], m_offspring[i + 1], pro, rnd);
			mutate(m_offspring[i], pro, rnd);
			mutate(m_offspring[i + 1], pro, rnd);
		}

		for (size_t i = 0; i < this->m_individuals.size(); i++) {
			//store location of population before change
			size_t eval = alg->evaluations();
			if ((eval + 1) % dynamic_cast<DMOPs*>(pro)->get_change_fre() == 0) {
				std::vector<Solution<>> temp_pop0;
				for (size_t j = 0; j < this->m_individuals.size() + i + 1; ++j) {
					if (j < i + 1)
						temp_pop0.emplace_back(m_offspring[j]);
					else
						temp_pop0.emplace_back(*this->m_individuals[j - i - 1]);
				}
				////output performance metrics
				//nondominatedSorting(temp_pop0);
				//std::vector<Solution<>> temp_pop;
				//for (size_t i = 0; i < temp_pop0.size(); ++i) {
				//	if (temp_pop0[i].fitness() == 0)
				//		temp_pop.emplace_back(temp_pop0[i]);
				//}
				//CONTINUOUS_CAST->get_optima().location_before_change(temp_pop);//for time-linkage
				//Real IGD = CONTINUOUS_CAST->get_optima().IGD_to_PF(temp_pop);
				//std::string file_name = static_cast<std::string>(g_working_dir) + "instance/algorithm/dynamic/dmoa/data/IGD4.txt";
				//std::ofstream out(file_name, std::ios::app);
				//if (out) {
				//	out << IGD << std::endl;
				//	out.close();
				//}
			}
			tag = m_offspring[i].evaluate(pro, alg);
			//if structure of the problems changed
			if (tag & kChangeVariableMemory || tag & kChangeObjectiveMemory) {
				handleEvaluationTag(tag);
			}
		}

		//detect changes or status
		if (1) {
			auto changed = problemChanged(m_detect_rate, pro, alg, rnd);
			if (changed) {
				responseChange(pro, alg, rnd);
			}
			else {
				//environment selection
				for (size_t i = 0; i < this->m_individuals.size(); ++i) {
					m_offspring[i + this->m_individuals.size()] = *this->m_individuals[i];
				}
				envirSelectionByNSGAII(*this,m_offspring,pro->optimizeMode());
				//survivorSelection(*this, m_offspring);

				//#ifdef OFEC_DEMO
				//				m_solution.clear();
				//				m_solution.resize(1);
				//				std::vector<std::vector<Solution<>*>> pops(1);
				//				for (size_t i = 0; i < this->m_individuals.size(); ++i)
				//					pops[0].emplace_back(&(this->m_individuals[i]));
				//				ofec_demo::g_buffer->appendAlgBuffer(alg);
				//#endif
			}
		}
		else {
			//environment selection
			for (size_t i = 0; i < this->m_individuals.size(); ++i) {
				m_offspring[i + this->m_individuals.size()] = *this->m_individuals[i];
			}
			envirSelectionByNSGAII(*this, m_offspring, pro->optimizeMode());
			//survivorSelection(*this, m_offspring);
			//#ifdef OFEC_DEMO
			//			std::vector<std::vector<Solution<>*>> pops(1);
			//			for (size_t i = 0; i < this->m_individuals.size(); ++i)
			//				pops[0].emplace_back(&(this->m_individuals[i]));
			//			ofec_demo::g_buffer->appendAlgBuffer(alg);
			//#endif
			if (populationConverged()) {
				addDiversity(pro, alg, rnd);
			}
		}
		//output objectives
		for (int i = 0; i < this->m_individuals.size(); i++) {
			std::cout << this->m_individuals[i]->objective()[0] << " " << this->m_individuals[i]->objective()[1] << std::endl;
		}

		m_iteration++;
		return tag;
	}

	bool DNSGAII_A_pop::problemChanged(Real r, Problem* pro, Algorithm* alg, Random* rnd) {
		Population<Solution<>> temp_pop;
		for (size_t i = 0; i < this->m_individuals.size(); ++i) {
			temp_pop.append(*(this->m_individuals[i]));
		}
		return ifProblemChanged(temp_pop, m_detect_rate, pro, alg, rnd);
	}

	bool DNSGAII_A_pop::populationConverged() {
		//CONTINUOUS_CAST->get_location_before_change();
		return ifPopConverged(this->m_individuals);
	}

	void DNSGAII_A_pop::addDiversity(Problem* pro, Algorithm* alg, Random* rnd) {

	}

	void DNSGAII_A_pop::responseChange(Problem* pro, Algorithm* alg, Random* rnd) {
		updatePop(pro, rnd);
		for (auto& i : this->m_individuals) {
			i->evaluate(pro, alg);
			//if (tag != evaluation_tag::Normal) break;
		}
		std::vector<std::vector<Real>*> objs2;
		for (auto& i : this->m_individuals)
			objs2.emplace_back(&i->objective());
		std::vector<int> rank;
		nd_sort::fastSort<Real>(objs2, rank, pro->optimizeMode());
		for (size_t i = 0; i < this->m_individuals.size(); ++i) {
			this->m_individuals[i]->setFitness(rank[i]);
		}
	}

	void DNSGAII_A_pop::updatePop(Problem* pro, Random* rnd) {
		size_t num = std::floor(this->m_individuals.size() * m_replace_rate);//get the number of replace Solutions
		std::vector<size_t> index;//store the index of the replace Solutions
		int tag = kNormalEval;
		for (size_t i = 0; i < num;) {
			size_t ind = std::floor(this->m_individuals.size() * rnd->uniform.next());
			if (index.empty() || *(index.end() - 1) != ind) {
				index.push_back(ind);
				this->m_individuals[ind]->initialize(pro, rnd);//re-initilize the replaced Solutions
				i++;
			}
		}
	}

}