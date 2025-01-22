#include "mqso.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void mQSO::initialize_()
	{
		MetricsDynamicConOEA::initialize_();
		auto& v = *m_param;
		mc_init_indiSize = v.has("population size") ? (v.get<int>("population size") / 10) * 10 : 100;
		double u, l;
		int d = 0;
		int peaks = 0;
		l = CAST_CONOP(m_problem.get())->range(0).first;
		u = CAST_CONOP(m_problem.get())->range(0).second;
		d = CAST_CONOP(m_problem.get())->numberVariables();
		if (CAST_CONOP(m_problem.get())->hasTag(ProblemTag::kDOP)) peaks = GET_DYNCONOP(m_problem.get())->getNumPeak();
		else peaks = 1;
		m_Rexcl = 0.5 * (u - l) / pow((double)peaks, 1. / d);		//paper function(14) r_excl= 0.5d_boa, we need to know how many peaks the problem have, it's not very well for Algorithm.
		m_Rconv = 1 * m_Rexcl;

		m_exclusion = true;
		m_M = mc_init_indiSize / 10;
	}

	void mQSO::initializeSolution(int num) {
		int rf = kNormalEval;
		std::unique_ptr<MQSOSwarm> newPop(new MQSOSwarm(num, m_problem.get()));
		newPop->initialize(this, m_problem.get(), m_random.get());
		rf = newPop->evaluate(m_problem.get(), this);
		if (rf & kChangeCurEval) {
			measureMultiPop(true);
			newPop->evaluate(m_problem.get(), this);
			newPop->initVelocity(m_problem.get(), m_random.get());
			newPop->initPbest(m_problem.get());
		}
		newPop->initVelocity(m_problem.get(), m_random.get());
		newPop->initPbest(m_problem.get());
		newPop->updateCurRadius(m_problem.get());
		newPop->getworstFlag() = false;
		m_sub_pop.append(newPop);
		//cout << "New Pop Initilized!" << endl;
	}

	bool mQSO::checkConvergenceAll()
	{
		for (int i = 0; i < m_M; i++) {
			if (m_sub_pop[i].getCurRadius() > m_Rconv) return false;
		}
		return true;
	}

	void mQSO::exclude()
	{
		if (m_exclusion) {
			for (int i = 0; i < m_M; i++) {
				for (int j = i + 1; j < m_M; j++) {
					auto best_i = m_sub_pop[i].bestIdx(m_problem.get());
					auto best_j = m_sub_pop[j].bestIdx(m_problem.get());
					if (m_sub_pop[i].getworstFlag() == false && m_sub_pop[j].getworstFlag() == false && m_sub_pop[i][best_i].variableDistance(m_sub_pop[j][best_j], m_problem.get()) < m_Rexcl) {
						count_exclusion++;
						if (m_sub_pop[i][best_i].pbest().dominate(m_sub_pop[j][best_j].pbest(), m_problem.get())) {
							m_sub_pop[j].getworstFlag() = true;
						}
						else m_sub_pop[i].getworstFlag() = true;
					}
				}
			}
		}
	}


	void mQSO::evolve(){
		int rf = kNormalEval;
		for (int i(0); i < m_M; i++) {
			if (!m_sub_pop[i].getworstFlag()) {
				rf = m_sub_pop[i].evolve(m_problem.get(), this, m_random.get());
				m_sub_pop[i].updateCurRadius(m_problem.get());
			}
			if (rf & kChangeCurEval) {
				measureMultiPop(true); break;
			}
		}
	}

	void mQSO::measureMultiPop(bool flag)
	{
		for (int i(0); i < m_sub_pop.size(); i++) {
			m_sub_pop[i].getworstFlag() = false;
			for (int j(0); j < m_sub_pop[i].size(); j++) {
				m_sub_pop[i][j].pbest().evaluate(m_problem.get(), flag ? this : nullptr);
			}
		}
	}

	void mQSO::run_() {
		for (int i(0); i < m_M; i++) initializeSolution(10);
		int iter = 0;
		while (!terminating()) {
			//std::cout << "iter#####" << iter++ << std::endl;
			evolve();
			updateSubPopBestIdx();
			m_bestsubPopIdx = updateBestSubPopFlag();
			exclude();
			autiConvergence();
			if (m_sub_pop[m_bestsubPopIdx].getworstFlag()) {
				std::cout << std::endl;
			}
			removeRedundancePop();
			//std::cout << std::endl;
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
		std::cout << "converge: " << count_converge << std::endl;
		std::cout << "exclude: " << count_exclusion << std::endl;
	}

	void mQSO::updateSubPopBestIdx() {
		for (int i(0); i < m_sub_pop.size(); i++) {
			m_sub_pop[i].updateBestIdx(m_sub_pop[i].bestIdx(m_problem.get()));
			//std::cout << i << ": " << "best Idx: " << m_sub_pop[i].bestIdx(m_problem.get()) << std::endl;
			//for (int j(0); j < m_sub_pop[i].size(); j++) {
			//	std::cout << m_sub_pop[i][j].objective()[0] << '\t' << std::endl;
			//}
			//std::cout << std::endl;
		}
	}

	int mQSO::updateBestSubPopFlag(){
		if (m_sub_pop.size() == 0) return -1;
		int idx = 0;
		for (int i(1); i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i][m_sub_pop[i].bestIdx(m_problem.get())].pbest().dominate(m_sub_pop[idx][m_sub_pop[idx].bestIdx(m_problem.get())].pbest(), m_problem.get())) {
				idx = i;
			}
		}
		for (int i(0); i < m_sub_pop.size(); i++) {
			if (i == idx) m_sub_pop[i].setBestSubPop(true);
			else m_sub_pop[i].setBestSubPop(false);
		}
		//std::cout << "best subpop: " << idx << std::endl;
		return idx;
	}

	int mQSO::findWorstPop()
	{
		if (m_sub_pop.size() == 0) return -1;
		int idx = 0;
		for (int i(1); i < m_sub_pop.size(); i++) {
			if (m_sub_pop[idx][m_sub_pop[idx].bestIdx(m_problem.get())].pbest().dominate(m_sub_pop[i][m_sub_pop[i].bestIdx(m_problem.get())].pbest(), m_problem.get())) {
				idx = i;
			}
		}
		return idx;
	}

	void mQSO::autiConvergence(){
		if (checkConvergenceAll()) {
			count_converge++;
			m_sub_pop[findWorstPop()].getworstFlag() = true;
		}
	}

	void mQSO::removeRedundancePop(){
		int rf = kNormalEval;
		for (int i(0); i < m_M; i++) {
			if (m_sub_pop[i].getworstFlag()) {
				rf = m_sub_pop[i].reInitialize(this, m_problem.get(), m_random.get());
				if (rf & kChangeCurEval) {
					measureMultiPop(true); break;
				}
			}
		}
	}

	void mQSO::record() {
		std::vector<Real> entry;
		entry.push_back(m_offline_error);
		//entry.push_back(m_evaluations);
		//entry.push_back(GET_DYNCONOP(m_problem.get())->getCount());
		//entry.push_back(GET_DYNCONOP(m_problem.get())->optima().objective(0)[0]);
		//entry.push_back(m_candidates.front()->objective(0));
		dynamic_cast<RecordVecRealDynamic*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void mQSO::updateBuffer() {
		m_solution.clear();
		m_solution.resize(m_sub_pop.size());
		for (size_t k = 0; k < m_sub_pop.size(); ++k) {
			for (size_t i = 0; i < m_sub_pop[k].size(); ++i) {
				m_solution[k].push_back(&m_sub_pop[k][i]);
			}
		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif
}