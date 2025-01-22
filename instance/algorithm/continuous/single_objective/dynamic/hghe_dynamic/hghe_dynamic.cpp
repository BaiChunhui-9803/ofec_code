#include "hghe_dynamic.h"
#include "../../../../core/problem/continuous/continuous.h"
#include "../../../../core/problem/solution.h"
#include "../../../record/rcr_vec_real.h"

namespace ofec {
	void HGHE_DYNAMIC::initialize_() {
		HGHE::initialize_();
		MetricsDynamicConOEA::initialize_();
		auto& v = *m_param;;
		m_dynamic_strategy = v.has("dynamic strategy") ? v.get<int>("dynamic strategy") : 1;
		auto adaptor = new AdaptorDynamicHGHEC(m_problem.get(), m_random.get(), this, 1);
		int tns = v.has("threshold number of solutions") ? v.get<int>("threshold number of solutions") : 100;
		int ins = v.has("interval number of solutions") ? v.get<int>("interval number of solutions") : 100;
		Real phi = v.has("phi") ? v.get<Real>("phi") : 6;
		auto pps = v.has("potential prediction strategy") ?
			static_cast<AdaptorDynamicHGHEC::PtnlPredStrat>(v.get<int>("potential prediction strategy")) :
			AdaptorDynamicHGHEC::PtnlPredStrat::kMean;
		adaptor->setPhi(phi);
		adaptor->setNumSolsThold(tns);
		adaptor->setNumSolsIntvl(ins);
		adaptor->setPtnlPredStrat(pps);
		m_adaptor.reset(adaptor);

		m_scaling_factor = v.has("scaling factor") ? v.get<Real>("scaling factor") : 0.5;
		m_crossover_rate = v.has("crossover rate") ? v.get<Real>("crossover rate") : 0.6;
		m_pop.resize(m_num_sols_exploit, m_problem.get());
		m_pop.setParameter(m_crossover_rate, m_scaling_factor);
		m_last_ner = m_last_nei = 0;
		m_keep_candidates_updated = true;
	}

	void HGHE_DYNAMIC::exploreHills(){
		if (m_pre_best_vars.empty()) {
			m_adaptor->calculateExplorationPotential();
			auto& ps = m_adaptor->getExplorationPotential();
			int hill;
			for (size_t i = 0; i < m_num_sols_explore; ++i) {
				hill = m_random->uniform.spinWheel(ps.begin(), ps.end()) - ps.begin();
				m_adaptor->exploreInHill(*m_sols_explore[i], hill);
			}
		}
		else{
			m_pre_best_vars.clear();
		}
	}

	void HGHE_DYNAMIC::updatePreBestVars(const std::vector<const SolutionBase*>& sols) {
		for (auto sol : sols) {
			m_pre_best_vars.push_back(adaptorDynamicHGHEC()->getSol(sol)->variable());
		}
	}

	void HGHE_DYNAMIC::reproduceInHill(size_t id_hill) {
		const auto& ids_pop = m_pop_each_hill[id_hill];
		const auto& ids_off = m_off_each_hill[id_hill];
		for (size_t i = 0; i < m_num_sols_exploit; ++i)
			m_pop[i] = dynamic_cast<Solution<>&>(*m_population[ids_pop[i]]);
		for (size_t i = 0; i < m_num_sols_exploit; ++i) {
			do {
				m_pop.mutate(i, m_random.get(), m_problem.get());
				m_pop.recombine(i, m_random.get(), m_problem.get());
				CAST_CONOP(m_problem.get())->validateSolution(m_pop[i].trial(), Validation::kSetToBound, m_random.get());
			} while (adaptorDynamicHGHEC()->theHillsLocated(m_pop[i].trial()).count(id_hill) == 0);
			dynamic_cast<Solution<>&>(*m_offspring[ids_off[i]]) = m_pop[i].trial();
		}
	}

	int HGHE_DYNAMIC::classifyPopulation(){
		int rf = kNormalEval;
		std::vector<decltype(m_population)> temp_inds(m_adaptor->numHills());
		if (!m_pre_best_vars.empty() && m_dynamic_strategy == 2) {
			/* Environment Change, use pre sols for explore */
			std::vector<int> inds_index(m_adaptor->numHills());
			for (size_t id_hill = 0; id_hill < m_adaptor->numHills(); ++id_hill) {
				inds_index[id_hill] = -1;
				temp_inds[id_hill].resize(m_num_sols_exploit);
				for (auto& sol_ptr : temp_inds[id_hill])
					m_adaptor->createSolution(sol_ptr);
			}
			for (const auto& pre_best_var : m_pre_best_vars) {
				const int id_hill = adaptorDynamicHGHEC()->theHillLocated(pre_best_var);
				auto id_hills = adaptorDynamicHGHEC()->theHillsLocated(pre_best_var);
				if (id_hill != -1) {
					adaptorDynamicHGHEC()->getSol(temp_inds[id_hill][++inds_index.at(id_hill)].get())->variable() = pre_best_var;
					rf = temp_inds[id_hill][inds_index.at(id_hill)]->evaluate(m_problem.get(), this);
					if ((rf & kTerminate) || (rf & kChangeNextEval))
						return rf;
					m_adaptor->archiveSolution(*temp_inds[id_hill][inds_index.at(id_hill)], TaskHGHE::kExplore);
				}
			}
			for (size_t id_hill = 0; id_hill < m_adaptor->numHills(); ++id_hill) {
				while(inds_index[id_hill] < m_num_sols_exploit - 1){
					m_adaptor->exploreInHill(*temp_inds[id_hill][++inds_index.at(id_hill)], id_hill);
					rf = temp_inds[id_hill][inds_index.at(id_hill)]->evaluate(m_problem.get(), this);
					if ((rf & kTerminate) || (rf & kChangeNextEval))
						return rf;
					m_adaptor->archiveSolution(*temp_inds[id_hill][inds_index.at(id_hill)], TaskHGHE::kExplore);
				}
			}
		}
		else {
			/* Update affiliation */
			// TODO 这里把所有的公共区域的解给丢弃了
			m_off_each_hill.clear();
			m_off_each_hill.resize(m_adaptor->numHills());
			for (size_t i = 0; i < m_offspring.size(); ++i) {
				int id_hill = m_adaptor->theHillLocated(*m_offspring[i]);
				if (id_hill > -1)
					m_off_each_hill[id_hill].push_back(i);
			}
			m_pop_each_hill.clear();
			m_pop_each_hill.resize(m_adaptor->numHills());
			for (size_t i = 0; i < m_population.size(); ++i) {
				int id_hill = m_adaptor->theHillLocated(*m_population[i]);
				if (id_hill > -1)
					m_pop_each_hill[id_hill].push_back(i);
			}
			/* Selection in each hill */
			for (size_t id_hill = 0; id_hill < m_adaptor->numHills(); ++id_hill) {
				std::vector<const SolutionBase*> cand;
				for (size_t id_ind : m_pop_each_hill[id_hill])
					cand.push_back(m_population[id_ind].get());
				for (size_t id_ind : m_off_each_hill[id_hill])
					cand.push_back(m_offspring[id_ind].get());
				temp_inds[id_hill].resize(m_num_sols_exploit);
				for (auto& sol_ptr : temp_inds[id_hill])
					m_adaptor->createSolution(sol_ptr);
				if (cand.size() > m_num_sols_exploit) {
					selectSurvivors(cand, temp_inds[id_hill]);
				}
				else {
					for (size_t i = 0; i < cand.size(); ++i)
						m_adaptor->replaceSolution(*temp_inds[id_hill][i], *cand[i]);
					for (size_t i = cand.size(); i < m_num_sols_exploit; ++i) {
						m_adaptor->exploreInHill(*temp_inds[id_hill][i], id_hill);
						rf = temp_inds[id_hill][i]->evaluate(m_problem.get(), this);
						if ((rf & kTerminate) || (rf & kChangeNextEval))
							return rf;
						// m_adaptor->archiveSolution(*temp_inds[id_hill][i], TaskHGHE::kExplore);
						m_adaptor->archiveSolution(*temp_inds[id_hill][i], TaskHGHE::kExploit);
					}
				}
			}
		}
		m_population.resize(m_adaptor->numHills() * m_num_sols_exploit);
		for (auto& sol_ptr : m_population)
			m_adaptor->createSolution(sol_ptr);
		for (size_t id_hill = 0; id_hill < m_adaptor->numHills(); ++id_hill) {
			m_pop_each_hill[id_hill].resize(m_num_sols_exploit);
			for (size_t i = 0; i < m_num_sols_exploit; ++i) {
				m_adaptor->replaceSolution(*m_population[id_hill * m_num_sols_exploit + i], *temp_inds[id_hill][i]);
				m_pop_each_hill[id_hill][i] = id_hill * m_num_sols_exploit + i;
			}
		}
		return rf;
	}

	AdaptorDynamicHGHEC* HGHE_DYNAMIC::adaptorDynamicHGHEC() const {
		return dynamic_cast<AdaptorDynamicHGHEC*>(m_adaptor.get());
	}

	int HGHE_DYNAMIC::evaluation(){
		int rf = kNormalEval;
		for (size_t i = 0; i < m_offspring.size(); ++i) {
			rf = m_offspring[i]->evaluate(m_problem.get(), this);
			if ((rf & kChangeNextEval) || (rf & kTerminate))
				return rf;
			m_adaptor->archiveSolution(*m_offspring[i], TaskHGHE::kExploit);
		}
		for (size_t i = 0; i < m_num_sols_explore; ++i) {
			rf = m_sols_explore[i]->evaluate(m_problem.get(), this);
			if ((rf & kChangeNextEval) || (rf & kTerminate))
				return rf;
			m_adaptor->archiveSolution(*m_sols_explore[i], TaskHGHE::kExplore);
		}
		return rf;
	}

	void HGHE_DYNAMIC::handleDynamicChange(int rf, int strategy){
		if(rf & kChangeNextEval){
			if (strategy == 1) {
				m_population.clear();
				m_offspring.clear();
				auto& v = *m_param;
				auto adaptor = new AdaptorDynamicHGHEC(m_problem.get(), m_random.get(), this, 1);
				int tns = v.has("threshold number of solutions") ? v.get<int>("threshold number of solutions") : 100;
				int ins = v.has("interval number of solutions") ? v.get<int>("interval number of solutions") : 100;
				Real phi = v.has("phi") ? v.get<Real>("phi") : 6;
				auto pps = v.has("potential prediction strategy") ?
					static_cast<AdaptorDynamicHGHEC::PtnlPredStrat>(v.get<int>("potential prediction strategy")) :
					AdaptorDynamicHGHEC::PtnlPredStrat::kMean;
				adaptor->setPhi(phi);
				adaptor->setNumSolsThold(tns);
				adaptor->setNumSolsIntvl(ins);
				adaptor->setPtnlPredStrat(pps);
				m_adaptor.reset(adaptor);
			}
			else if(strategy == 2){
				std::vector<const SolutionBase*> best_sols;
				best_sols = adaptorDynamicHGHEC()->getBestSolsEachHill();
				updatePreBestVars(best_sols);
				adaptorDynamicHGHEC()->clearPreHis();
				m_population.clear();
				m_offspring.clear();
			}
		}
	}

	void HGHE_DYNAMIC::selectSurvivors(
		std::vector<const SolutionBase *> &candidates, 
		std::vector<std::unique_ptr<SolutionBase>> &population)
	{
		for (size_t i = 0; i < m_num_sols_exploit; ++i) {
			auto pu = i + m_num_sols_exploit;
			if (pu < candidates.size() && candidates[pu]->dominate(*candidates[i], m_problem.get()))
				m_adaptor->replaceSolution(*population[i], *candidates[pu]);
			else
				m_adaptor->replaceSolution(*population[i], *candidates[i]);
		}
	}

	void HGHE_DYNAMIC::run_() {
		m_sols_explore.resize(m_num_sols_explore);
		for (auto& sol_ptr : m_sols_explore)
			m_adaptor->createSolution(sol_ptr);
		while (!terminating()) {
			m_adaptor->updateHills();
			int rf = exploitHills();
			if (rf & kChangeNextEval) {
				handleDynamicChange(rf, m_dynamic_strategy);
				continue;
			}
			if (rf & kTerminate) break;
			exploreHills();
			rf = evaluation();
			if (rf & kChangeNextEval) {
				handleDynamicChange(rf, m_dynamic_strategy);
				continue;
			}
			if (rf & kTerminate) break;
#ifdef OFEC_DEMO
			updateBuffer();
#endif 
		}
	}

	void HGHE_DYNAMIC::record() {

	}
}