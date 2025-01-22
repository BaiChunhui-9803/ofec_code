#include "samo.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void SAMO::initialize_() {
		MetricsDynamicConOEA::initialize_();
		auto& v = *m_param;
		m_M = 0;
		m_Mfree = 0;
		m_Necxess = v.has("number of free populations") ? v.get<int>("number of free populations") : 1;

		m_exclusion = true;
	}

	void SAMO::initializeSolution(int num) {
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
		m_M++;
		m_Mfree++;
	}

	void SAMO::exclude() {
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
						else{
							m_sub_pop[i].getworstFlag() = true;
						}
					}
				}
			}
		}
	}

	void SAMO::measureMultiPop(bool flag) {
		for (int i(0); i < m_M; i++) {
			m_sub_pop[i].getworstFlag() = false;
			for (int j(0); j < m_sub_pop[i].size(); j++) {
				m_sub_pop[i][j].pbest().evaluate(m_problem.get(), flag ? this : nullptr);
			}
		}
	}

	void SAMO::updateRexcludeRaduis() {
		double u, l;
		int d = 0;
		int peaks = 0;
		l = CAST_CONOP(m_problem.get())->range(0).first;
		u = CAST_CONOP(m_problem.get())->range(0).second;
		d = CAST_CONOP(m_problem.get())->numberVariables();
		m_Rexcl = 0.5 * (u - l)/(pow(m_M, 1./d));
		m_Rconv = 2 * m_Rexcl;
	}

	void SAMO::evolve() {
		int rf = kNormalEval;
		for (int i(0); i < m_M; i++) {
			if (!m_sub_pop[i].getworstFlag()) {
				rf = m_sub_pop[i].evolve(m_problem.get(), this, m_random.get());
				m_sub_pop[i].updateCurRadius(m_problem.get());
			}
			else if (m_sub_pop[i].getworstFlag()) {
				rf = reInitializePop(i);
				if (m_sub_pop[i].getswarmType() == MQSOSwarm::SwarmType::SWARM_CONVERGING) {
					m_Mfree++;
				}
				m_sub_pop[i].setswarmType(MQSOSwarm::SwarmType::SWARM_FREE);
			}
			if (rf & kChangeCurEval) {
				measureMultiPop(true); break;
			}
		}
	}

	void SAMO::run_() {
		initializeSolution(10);
		while (!terminating()) {
			updateRexcludeRaduis();
			checkSwarmConverge();
			exclude();
			evolve();
			//removeRedundancePop();
			checkSwarmState();
#ifdef OFEC_DEMO
 			updateBuffer();
#endif
			}
		}

	int SAMO::findWorstPop() {
		if (m_M == 0) return -1;
		int idx = 0;
		while (m_sub_pop[idx].getswarmType() != MQSOSwarm::SwarmType::SWARM_FREE) idx++;
		for (int i(idx+1); i < m_M; i++) {
			if (m_sub_pop[i].getswarmType() != MQSOSwarm::SwarmType::SWARM_FREE) continue;
			if (m_sub_pop[idx][m_sub_pop[idx].bestIdx(m_problem.get())].pbest().dominate(m_sub_pop[i][m_sub_pop[i].bestIdx(m_problem.get())].pbest(), m_problem.get())) {
				idx = i;
			}
		}
		return idx;
	}

	void SAMO::checkSwarmConverge() {
		int num_converging = 0;
		bool converging;
		for (int i = 0; i < m_M; i++) {
			converging = true;
			for (int j = 0; j < m_sub_pop[i].size(); j++) {
				for (int k = j + 1; k < m_sub_pop[i].size(); k++) {
					if (m_sub_pop[i][j].getType() != MQSOParticle::ParticleType::PARTICLE_NEUTRAL || m_sub_pop[i][k].getType() != MQSOParticle::ParticleType::PARTICLE_NEUTRAL) continue;

					if (fabs(m_sub_pop[i][j].variableDistance(m_sub_pop[i][k],m_problem.get())) > m_Rconv) {
						// a free swarm becomes a converging swarm when the distance of two farthest particles is less than ms_Rconv
						converging = false;
						break;
					}

				}
				if (!converging) break;
			}
			if (converging && m_sub_pop[i].getswarmType()== MQSOSwarm::SwarmType::SWARM_FREE) {
				m_sub_pop[i].setswarmType(MQSOSwarm::SwarmType::SWARM_CONVERGING);
				m_Mfree--;
			}
		}
	}

	void SAMO::checkSwarmConverge2() {
		int count = 0;
		//updateRexcludeRaduis();
		for (size_t i(0); i < m_M; i++) {
			if (m_sub_pop[i].getCurRadius() < m_Rconv) {
				if (m_sub_pop[i].getswarmType() == MQSOSwarm::SwarmType::SWARM_FREE) {
					m_sub_pop[i].setswarmType(MQSOSwarm::SwarmType::SWARM_CONVERGING);
					m_Mfree--;
				}
			}
		}

		int num_converging = 0;
		bool converging;
		for (int i = 0; i < m_M; i++) {
			converging = true;
			for (int j = 0; j < m_sub_pop[i].size(); j++) {
				for (int k = j + 1; k < m_sub_pop[i].size(); k++) {
					if (m_sub_pop[i][j].getType() != MQSOParticle::ParticleType::PARTICLE_NEUTRAL || m_sub_pop[i][k].getType() != MQSOParticle::ParticleType::PARTICLE_NEUTRAL) continue;

					if (fabs(m_sub_pop[i][j].variableDistance(m_sub_pop[i][k], m_problem.get())) > m_Rconv) {
						// a free swarm becomes a converging swarm when the distance of two farthest particles is less than ms_Rconv
						converging = false;
						break;
					}

				}
				if (!converging) break;
			}
			if (converging && m_sub_pop[i].getswarmType() == MQSOSwarm::SwarmType::SWARM_FREE) {
				m_sub_pop[i].setswarmType(MQSOSwarm::SwarmType::SWARM_CONVERGING);
				m_Mfree--;
			}
		}
	}

	void SAMO::checkSwarmState() {
		if(m_Mfree<=0) {
			initializeSolution(10);
		}
		while(m_Mfree > m_Necxess) {
			int idx = findWorstPop();
			m_sub_pop.remove(m_sub_pop.begin() + idx);
			m_Mfree--;
			m_M--;
		}
	}

	int SAMO::reInitializePop(int idx) {
		int rf = kNormalEval;
		rf = m_sub_pop[idx].reInitialize(this, m_problem.get(), m_random.get());
		return rf;
	}

	void SAMO::removeRedundancePop(){
		int rf = kNormalEval;
		for (int i = 0; i < m_M; i++) {
			if (m_sub_pop[i].getworstFlag()) {
				rf = reInitializePop(i);
				if (m_sub_pop[i].getswarmType() == MQSOSwarm::SwarmType::SWARM_CONVERGING) {
					m_Mfree++;
				}
				m_sub_pop[i].setswarmType(MQSOSwarm::SwarmType::SWARM_FREE);
				if (rf & kChangeCurEval) {
					measureMultiPop(true);
				}
			}
		}
	}

	void SAMO::record() {
		std::vector<Real> entry;
		entry.push_back(m_offline_error);
		//entry.push_back(m_evaluations);
		//entry.push_back(GET_DYNCONOP(m_problem.get())->getCount());
		//entry.push_back(GET_DYNCONOP(m_problem.get())->optima().objective(0)[0]);
		//entry.push_back(m_candidates.front()->objective(0));
		dynamic_cast<RecordVecRealDynamic*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void SAMO::updateBuffer() {
		m_solution.clear();
		m_solution.resize(m_M);
		for (size_t k = 0; k < m_M; ++k) {
			for (size_t i = 0; i < m_sub_pop[k].size(); ++i) {
				m_solution[k].push_back(&m_sub_pop[k][i]);
			}
		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif
}