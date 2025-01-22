#include "cpsor.h"
#include "../../../../../record/rcr_vec_real.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void CPSOR::initialize_() {
		MetricsDynamicConOEA::initialize_();
		auto& v = *m_param;
		//int numvar = GET_DYNCONOP(m_problem.get())->numberVariables();
		//int numpeak = GET_DYNCONOP(m_problem.get())->getNumPeak();
		m_init_popsize = v.has("population size") ? v.get<int>("population size") : 70;
		m_overlapDegree = v.has("overlap degree") ? v.get<Real>("overlap degree") : 0.1;
		m_rconv = v.has("converge radius") ? v.get<Real>("converge radius") : 10e-3;
		m_alpha = v.has("diversity degree") ? v.get<Real>("diversity degree") : 0.3;
		m_max_subpopSize = v.has("max subpop size") ? v.get<int>("max subpop size") : 3;
	}

	void CPSOR::measureMultiPop(bool flag){
		for (int i(0); i < m_sub_pop.size(); ++i) {
			for (int j(0); j < m_sub_pop[i].size(); j++) {
				m_sub_pop[i][j].pbest().evaluate(m_problem.get(), flag ? this : nullptr);
			}
		}
		for (int i(0); i < m_all_indis.size(); ++i) {
			m_all_indis[i]->pbest().evaluate(m_problem.get(), flag ? this : nullptr);
		}
	}

	int CPSOR::initializeOriginalPop(int num) {
		int rf = kNormalEval;
		const size_t num_obj = CAST_CONOP(m_problem.get())->numberObjectives();
		const size_t num_con = CAST_CONOP(m_problem.get())->numberConstraints();
		const size_t num_var = CAST_CONOP(m_problem.get())->numberVariables();
		int origin_size = m_all_indis.size();
		m_all_indis.resize(num + origin_size);

		for (int i(origin_size); i < m_all_indis.size(); ++i) {
			m_all_indis[i].reset(new templateParticle(num_obj, num_con, num_var));
			m_all_indis[i]->initialize(i, m_problem.get(), m_random.get());
			m_all_indis[i]->initializeVmax(m_problem.get(), m_random.get());
			m_all_indis[i]->initVelocity(m_problem.get(), m_random.get());
			rf = m_all_indis[i]->evaluate(m_problem.get(), this);
			if (rf & kChangeCurEval) return rf;
			m_all_indis[i]->pbest() = m_all_indis[i];
		}

		return rf;
	}

	int CPSOR::createSubswarms(int num)
	{
		int rf = kNormalEval;
		rf = initializeOriginalPop(num);
		if (rf & kChangeCurEval) {
			m_all_indis.clear();
			return rf;
		}
		HSLH<templateParticle> cluster(m_all_indis,m_problem.get());
		cluster.clusteringRough(m_max_subpopSize,m_min_subpopSize,m_problem.get());
		for(int i(0);i<cluster.size();i++){
			std::unique_ptr<CPSOSwarm> new_pop = std::make_unique<CPSOSwarm>(cluster[i].size(), m_problem.get());
			auto iter = cluster[i].begin();
			for (int j(0); j < cluster[i].size(); ++j, ++iter) {
				(*new_pop)[j] = *(iter->second);
			}
			new_pop->initializeParameters(m_problem.get(), this, m_random.get());
			new_pop->updateCurRadius(m_problem.get());
			new_pop->initializeRadius();
			m_sub_pop.append(new_pop);
		}
#ifdef OFEC_DEMO
				for (size_t i(0); i < m_sub_pop.size(); i++) m_gbestslocation.push_back(m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest());
				updateBuffer();
#endif
		return rf;
	}

	int CPSOR::checkOverlapping() {
		for (size_t i = 0; i < m_sub_pop.size(); i++) {
			if (this->m_sub_pop[i].getConverge()) continue;
			if (this->m_sub_pop[i].size() == 0) continue;
			for (size_t j = i + 1; j < m_sub_pop.size(); j++) {
				if (this->m_sub_pop[j].getConverge()) continue;
				if (this->m_sub_pop[j].size() == 0) continue;
				const int best_i_inx = m_sub_pop[i].best_idx(m_problem.get());
				const int best_j_inx = m_sub_pop[j].best_idx(m_problem.get());
				Real dist = m_problem->variableDistance(m_sub_pop[i][best_i_inx], m_sub_pop[j][best_j_inx]);
				if (dist < m_sub_pop[i].getInitCurRadius() && dist < m_sub_pop[j].getInitCurRadius()) {
					int c1(0);
					int c2(0);
					for (size_t k = 0; k < m_sub_pop[j].size(); k++) {
						dist = m_problem->variableDistance(*m_sub_pop[i].center(), m_sub_pop[j][k]);
						if (dist < m_sub_pop[i].getInitCurRadius()) c1++;
					}

					for (size_t k = 0; k < m_sub_pop[i].size(); k++) {
						dist = m_problem->variableDistance(*m_sub_pop[j].center(), m_sub_pop[i][k]);
						if (dist < m_sub_pop[i].getInitCurRadius()) c2++;
					}
					if(c1 > m_overlapDegree * m_sub_pop[i].size() && c2 > m_overlapDegree * m_sub_pop[j].size()){
						int idx = -1;
						if (m_sub_pop[i][best_i_inx].pbest().dominate(m_sub_pop[j][best_j_inx].pbest(), m_problem.get())) {
							m_sub_pop[i].merge(m_sub_pop[j],m_problem.get(), m_random.get());
							m_sub_pop.remove(m_sub_pop.begin() + j);
							idx = j;
						}
						else{
							m_sub_pop[j].merge(m_sub_pop[i], m_problem.get(), m_random.get());
							m_sub_pop.remove(m_sub_pop.begin() + i);
							idx = i;
						}
						return idx;
					}
				}
			}
		}
		return -1;
	}

	void CPSOR::checkConverging() {
		for (size_t i(0); i < m_sub_pop.size(); i++){
			if (m_sub_pop[i].getConverge()) continue;
			if(m_sub_pop[i].getCurRadius() <m_rconv){
				m_sub_pop[i].setConverge(true);
			}
		}
	}

	void CPSOR::run_() {
		createSubswarms(m_init_popsize);
		while (!terminating()) {
#ifdef OFEC_DEMO
			m_gbestslocation.clear();
#endif
			int rf = kNormalEval;

			int total = 0;
			for (size_t i(0); i < m_sub_pop.size(); i++) {
				if(m_sub_pop[i].getConverge()) continue;
				rf = m_sub_pop[i].localSearch(m_problem.get(), this, m_random.get(), m_init_popsize);	//localSearch
				total += m_sub_pop[i].size();
				if (rf != kNormalEval) break;
			}

			if (rf != kNormalEval) {
				informChange(rf);
				continue;
			}

			while (checkOverlapping() != -1);
			for (int i(0); i < m_sub_pop.size(); i++) m_sub_pop[i].checkOverCrowd(m_max_subpopSize, m_problem.get(), m_random.get());

			checkConverging();

			int num_survived_inidividuals(0);
			for (int i(0); i < m_sub_pop.size(); i++) {
				if (!m_sub_pop[i].getConverge()) num_survived_inidividuals += m_sub_pop[i].size();
			}

			if(num_survived_inidividuals < m_init_popsize * m_alpha){
				m_all_indis.clear();
				for (size_t i(0); i < m_sub_pop.size(); i++) {
					if (m_sub_pop[i].getConverge()) {
						m_sub_pop.remove(m_sub_pop.begin() + i);
						i--;
					}
				}
				int size = m_init_popsize - num_survived_inidividuals;
				rf = createSubswarms(size);
				if (rf != kNormalEval) {
					informChange(rf);
					continue;
				}
			}
#ifdef OFEC_DEMO
			for (size_t i(0); i < m_sub_pop.size(); i++) m_gbestslocation.push_back(m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest());
			updateBuffer();
#endif
		}
	}

	void CPSOR::informChange(int rf){
		if (rf & kChangeCurEval) {
			measureMultiPop(true);
		}
	}

	void CPSOR::record(){
		//std::vector<Real> entry;
		//entry.push_back(m_offline_error);
		//entry.push_back(m_best_before_change_error);
		////entry.push_back(m_evaluations);
		////entry.push_back(GET_DYNCONOP(m_problem.get())->getCount());
		////entry.push_back(GET_DYNCONOP(m_problem.get())->optima().objective(0)[0]);
		////entry.push_back(m_candidates.front()->objective(0));
		//dynamic_cast<RecordVecRealDynamic*>(m_record.get())->record(this, entry);

		std::vector<Real> entry(7, 0);
		entry[0] = m_evaluations;
		for (auto &sol : m_candidates) {
			auto &var = dynamic_cast<VariableVector<Real>&>(sol->variableBase());
			if (var[0] < exp(-5 * OFEC_PI / 20))
				entry[1]++;
			else if (var[0] < exp(-1 * OFEC_PI / 20))
				entry[2]++;
			else if (var[0] < exp(3 * OFEC_PI / 20))
				entry[3]++;
			else if (var[0] < exp(7 * OFEC_PI / 20))
				entry[4]++;
			else if (var[0] < exp(11 * OFEC_PI / 20))
				entry[5]++;
			else
				entry[6]++;
		}
		dynamic_cast<RecordVectorReal*>(m_record.get())->addEntry(this, entry);
	}

#ifdef OFEC_DEMO
	void CPSOR::updateBuffer() {
		m_solution.clear();
		m_solution.resize(m_sub_pop.size());
		for (size_t k = 0; k < m_sub_pop.size(); ++k) {
			for (size_t i = 0; i < m_sub_pop[k].size(); ++i) {
				m_solution[k].push_back(&m_sub_pop[k][i]);
			}
		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}

	std::vector<bool> CPSOR::getPopHiberState() const {
		std::vector<bool> hiber_state;
		for (size_t i(0); i < m_sub_pop.size(); i++) {
			hiber_state.push_back(m_sub_pop[i].getConverge());
		}
		return hiber_state;
	}
#endif
}
