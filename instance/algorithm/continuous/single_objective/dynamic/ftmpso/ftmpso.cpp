#include "ftmpso.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void FTMPSO::initialize_() {
		MetricsDynamicConOEA::initialize_();
		m_xr.reset(new Particle(m_problem->numberObjectives(), m_problem->numberConstraints(), CAST_CONOP(m_problem.get())->numberVariables()));
		auto& v = *m_param;
		int numvar = GET_DYNCONOP(m_problem.get())->numberVariables();
		int numpeak = GET_DYNCONOP(m_problem.get())->getNumPeak();
		m_etry = 20;
		m_rcloud = 0.2;
		m_cfmin = 0.8;
		m_rconv = v.has("Converge Radius") ? v.get<Real>("Converge Radius") : 1.;
		m_rexcl = v.has("Exclusion Radius") ? v.get<Real>("Exclusion Radius") : 0.5*(100./pow(numpeak, (1./numvar)));
		m_rs = v.has("Search Radius") ? v.get<Real>("Search Radius") : 0.5;
		m_qs = v.has("Velocity Search") ? v.get<Real>("Velocity Search") : 0.5;
		m_limit = v.has("Hibernate Limit") ? v.get<Real>("Hibernate Limit") : 0.4;
	}

	void FTMPSO::initializeFinder(int num){
		int rf = kNormalEval;
		auto newPop = std::make_unique<FTMPSOSwarm>(num, m_problem.get(), FTMPSOSwarm::SwarmType::FINDER);
		newPop->initialize(m_problem.get(), m_random.get());
		rf = newPop->evaluate(m_problem.get(), this);
		newPop->initVelocity(m_problem.get(), m_random.get());
		newPop->initPbest(m_problem.get());
		informChange(rf);
		m_sub_pop.append(newPop);
	}

	int FTMPSO::initializeTracker(FTMPSOSwarm& parent, int num){
		int rf = kNormalEval;
		sortFinder(parent);
		std::unique_ptr<FTMPSOSwarm> newPop(new FTMPSOSwarm(num, m_problem.get(), FTMPSOSwarm::SwarmType::TRACKER));
		newPop->weight() = 0.729844;
		newPop->accelerator1() = 1.49618;
		newPop->accelerator2() = 1.49618;
		newPop->setNeighborhood(m_random.get());
		for (size_t i(0); i < newPop->size(); i++)	(*newPop)[i] = parent[parent.size() - i - 1];
		m_sub_pop.append(newPop);
		parent.initialize(m_problem.get(), m_random.get());
		rf = parent.evaluate(m_problem.get(), this);
		parent.initPbest(m_problem.get());
		informChange(rf);
		return rf;
	}

	void FTMPSO::sortFinder(FTMPSOSwarm& parent){
		Particle temp(m_problem->numberObjectives(), m_problem->numberConstraints(), CAST_CONOP(m_problem.get())->numberVariables());
		//std::cout << "Before sort:" << std::endl;
		//for (int i = 0; i < parent.size(); i++) {
		//	for (int j = 0; j < parent[i].pbest().variable().size(); j++) {
		//		std::cout << parent[i].pbest().variable()[j] << '\t';
		//	}
		//	std::cout << '\n' << parent[i].pbest().objective()[0] << '\n';
		//}

		for (int i = parent.size() - 1; i > 0; i--) {
			for (int j = 0; j < i; j++) {
				if (parent[j].pbest().dominate(parent[j + 1].pbest(),m_problem.get())) {
					temp = m_sub_pop[0][j + 1];
					parent[j + 1] = m_sub_pop[0][j];
					parent[j] = temp;
				}
			}
		}

		//std::cout << "After sort:" << std::endl;
		//for (int i = 0; i < parent.size(); i++) {
		//	for (int j = 0; j < parent[i].pbest().variable().size(); j++) {
		//		std::cout << parent[i].pbest().variable()[j] << '\t';
		//	}
		//	std::cout << '\n' << parent[i].pbest().objective()[0] << '\n';
		//}
	}

	int FTMPSO::getBestTracker(){
		int bestindex = -1;
		for (size_t i(1); i < m_sub_pop.size(); i++) {
			if (bestindex == -1) bestindex = i;
			if (m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest().dominate(m_sub_pop[bestindex][m_sub_pop[bestindex].best_idx(m_problem.get())].pbest(), m_problem.get())) {
				bestindex = i;
			}
		}
		return bestindex;
	}

	void FTMPSO::updateFinderGBestMemory(){
		if (m_finderGbests.size() < 3) m_finderGbests.push_back(m_sub_pop[0][m_sub_pop[0].best_idx(m_problem.get())].pbest());
		else {
			for(size_t i(0);i< m_finderGbests.size()-1;i++)
				m_finderGbests[i] = m_finderGbests[i+1];
			m_finderGbests[m_finderGbests.size() - 1] = m_sub_pop[0][m_sub_pop[0].best_idx(m_problem.get())].pbest();
		}
	}

	int FTMPSO::finderExclusion(){
		int rf = kNormalEval;
		for (size_t i = 1; i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i].getexcelFlag()) continue;
			if(m_sub_pop[0][m_sub_pop[0].best_idx(m_problem.get())].pbest().variableDistance(m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest(), m_problem.get()) < m_rexcl) {
				m_sub_pop[0].initialize(m_problem.get(), m_random.get());
				rf = m_sub_pop[0].evaluate(m_problem.get(), this);
				m_sub_pop[0].initPbest(m_problem.get());
				informChange(rf);									//inform change if there are evaluations.
				if (rf != kNormalEval) return rf;
			}
		}
		return rf;
	}

	void FTMPSO::trackerExclusion(){
		for (size_t i = 1; i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i].getexcelFlag()) continue;
			for (size_t j = i + 1; j < m_sub_pop.size(); j++) {
				if (m_sub_pop[j].getexcelFlag()) continue;
				if (m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest().variableDistance(m_sub_pop[j][m_sub_pop[j].best_idx(m_problem.get())].pbest(), m_problem.get()) < m_rexcl) {
					if (m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest().dominate(m_sub_pop[j][m_sub_pop[j].best_idx(m_problem.get())].pbest(), m_problem.get())) {
						m_sub_pop[j].setexcelFlag(true);
					}
					else m_sub_pop[i].setexcelFlag(true);
				}
			}
		}
	}

	void FTMPSO::removeExclusion() {
		if (m_sub_pop.size() < 2) return;
		for (size_t i = 1; i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i].getexcelFlag()) {
				m_sub_pop.remove(m_sub_pop.begin() + i);
				i--;
			}
		}
	}

	int FTMPSO::isFinderInTracker(){
		for (decltype(m_sub_pop.size()) i = 1; i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i].getexcelFlag()) continue;
			if (m_sub_pop[0][m_sub_pop[0].best_idx(m_problem.get())].pbest().variableDistance(m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest(), m_problem.get()) < m_rexcl) {
				return i;
			}
		}
		return -1;
	}

	int FTMPSO::activation(){
		int rf = kNormalEval;
		updateFinderGBestMemory();
		if (m_finderGbests.size() < 3) return rf;
		if (m_sub_pop[0].iteration() % 2 == 0) {
			if (m_finderGbests[2].variableDistance(m_finderGbests[0], m_problem.get()) < m_rconv) {
				int idx = isFinderInTracker();
				if (idx == -1) {
					rf = initializeTracker(m_sub_pop[0], 5);
					if (rf != kNormalEval) return rf;
				}
				else if (idx != -1) {
					m_sub_pop[idx].setexcelFlag(true);
				}
			}
		}
		return rf;
	}

	int FTMPSO::trackerEvolve() {
		int rf = kNormalEval;
		if (m_sub_pop.size() > 1) {
			int bestIdx = getBestTracker();
			for (size_t i(1); i < m_sub_pop.size(); i++) {
				if (m_sub_pop[i].gethiberFlag() || m_sub_pop[i].getexcelFlag()) continue;
				rf = m_sub_pop[i].evolve(m_problem.get(), this, m_random.get());
				informChange(rf); if (rf != kNormalEval) return rf;
				if (i != bestIdx && judgeHiberPop(m_sub_pop[i])) m_sub_pop[i].sethiberFlag(true);
			}
		}
		return rf;
	}

	int FTMPSO::besttrackerTry(){
		int rf = kNormalEval;
		int best_popindex = getBestTracker();
		if (best_popindex != -1) {
			int effectivecount(0);
			int best_indivindex = m_sub_pop[best_popindex].best_idx(m_problem.get());
			for (size_t i(0); i < m_etry; i++) {
				for (size_t j(0); j < CAST_CONOP(m_problem.get())->numberVariables(); j++) {
					(*m_xr).variable()[j] = m_sub_pop[best_popindex][best_indivindex].pbest().variable()[j] +
						m_random->uniform.nextNonStd<Real>(-1, 1) * m_rcloud;
				}
				rf = (*m_xr).evaluate(m_problem.get(), this,true);
				if ((*m_xr).dominate(m_sub_pop[best_popindex][best_indivindex].pbest(),m_problem.get())) {
					m_sub_pop[best_popindex][best_indivindex] = (*m_xr);
					m_sub_pop[best_popindex][best_indivindex].pbest() = m_sub_pop[best_popindex][best_indivindex];
					effectivecount++;
				}
				informChange(rf);
				if (rf != kNormalEval)  break;
			}
			m_cf = m_cfmin * m_random->uniform.next() * (1 - m_cfmin);
			m_rcloud *= m_cf;
			//cout << effectivecount << endl;
		}
		return rf;
	}

	bool FTMPSO::judgeHiberPop(FTMPSOSwarm& pop){
		bool hiberFlag = true;
		for (size_t i(0); i < pop.size(); i++) {
			for (size_t j(0); j < pop[i].velocity().size(); j++) {
				if (fabs(pop[i].velocity()[j]) > m_limit) {
					hiberFlag = false;
					break;
				}
			}
		}
		return hiberFlag;
	}

	void FTMPSO::wakeupHiberPop(){
		for (size_t i(1); i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i].gethiberFlag()) {
				m_sub_pop[i].sethiberFlag(false);
			}
			relaunchTrackerSwarm(m_sub_pop[i]);
		}
	}

	void FTMPSO::relaunchTrackerSwarm(FTMPSOSwarm& tracker){
		int rf = kNormalEval;
		int best_idx = tracker.best_idx(m_problem.get());
		for (size_t i(0); i < tracker.size(); i++) {
			if (i != best_idx) {
				for (size_t j(0); j < tracker[i].variable().size(); j++) {
					tracker[i].variable()[j] = tracker[best_idx].pbest().variable()[j] + m_random->uniform.nextNonStd<Real>(-1, 1) * m_rs;
					tracker[i].velocity()[j] = m_random->uniform.nextNonStd<Real>(-1, 1) * m_qs;
				}
				rf = tracker[i].evaluate(m_problem.get(), this);
				tracker[i].pbest() = tracker[i];
				informChange(rf);
			}
		}
	}

	void FTMPSO::informChange(int rf){
		if (rf & kChangeCurEval) {
			m_finderGbests.clear();
			measureMultiPop(true);
			wakeupHiberPop();
		}
	}

	void FTMPSO::measureMultiPop(bool flag){
		m_rcloud = 0.2;
		for (int i(0); i < m_sub_pop.size(); ++i) {
			for (int j(0); j < m_sub_pop[i].size(); j++) {
				m_sub_pop[i][j].pbest().evaluate(m_problem.get(), this, flag);
			}
		}
	}

	void FTMPSO::run_() {
		initializeFinder(10);
		while (!terminating()) {
#ifdef OFEC_DEMO
			m_gbestslocation.clear();
#endif
			int rf = kNormalEval;
			rf = m_sub_pop[0].evolve(m_problem.get(), this, m_random.get());	//finder evolve
			informChange(rf);		//inform environment change
			finderExclusion();		//remove exclusion finder swarm
			activation();		//activate tracker swarm
			trackerEvolve();							//tracker evolve
			besttrackerTry();
			trackerExclusion();
			removeExclusion();
#ifdef OFEC_DEMO
			for (size_t i(0); i < m_sub_pop.size(); i++) m_gbestslocation.push_back(m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest());
			updateBuffer();
#endif
		}
	}

	void FTMPSO::record(){
		std::vector<Real> entry;
		entry.push_back(m_offline_error);
		//entry.push_back(m_evaluations);
		//entry.push_back(GET_DYNCONOP(m_problem.get())->getCount());
		//entry.push_back(GET_DYNCONOP(m_problem.get())->optima().objective(0)[0]);
		//entry.push_back(m_candidates.front()->objective(0));
		dynamic_cast<RecordVecRealDynamic*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void FTMPSO::updateBuffer(){
		m_solution.clear();
		m_solution.resize(m_sub_pop.size());
		for (size_t k = 0; k < m_sub_pop.size(); ++k) {
			for (size_t i = 0; i < m_sub_pop[k].size(); ++i) {
				m_solution[k].push_back(&m_sub_pop[k][i]);
			}
		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}

	std::vector<bool> FTMPSO::getPopHiberState() const {
		std::vector<bool> hiber_state;
		for (size_t i(0); i < m_sub_pop.size(); i++) {
			hiber_state.push_back(m_sub_pop[i].gethiberFlag());
		}
		return hiber_state;
	}
#endif
}
