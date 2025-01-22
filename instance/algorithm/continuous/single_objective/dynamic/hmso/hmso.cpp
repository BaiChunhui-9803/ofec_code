#include "hmso.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void HmSO::initialize_() {
		MetricsDynamicConOEA::initialize_();
		auto& v = *m_param;
		m_preiter.reset(new Solution<>(m_problem->numberObjectives(), m_problem->numberConstraints(), CAST_CONOP(m_problem.get())->numberVariables()));
		m_curiter.reset(new Solution<>(m_problem->numberObjectives(), m_problem->numberConstraints(), CAST_CONOP(m_problem.get())->numberVariables()));
		m_rpc = v.has("Radius between Parent and Child") ? v.get<Real>("Radius between Parent and Child") : 30.;
		m_rconv = v.has("Converge Radius") ? v.get<Real>("Converge Radius") : 1.0;
		m_rexcl = v.has("Exclusion Radius") ? v.get<Real>("Exclusion Radius") : 30.0;
		m_rs = v.has("Search Radius") ? v.get<Real>("Search Radius") : 0.5;
		m_margin = 5.;
	}

	void HmSO::initializeParent(int num)
	{
		int rf = kNormalEval;
		std::unique_ptr<HmSOSwarm> newPop(new HmSOSwarm(num, m_problem.get(), HmSOSwarm::SwarmType::Parent));
		for (size_t i(0); i < newPop->size(); i++) (*newPop)[i].initialize(i, m_problem.get(), m_random.get());
		newPop->weight() = 0.729844;
		newPop->accelerator1() = 1.49618;
		newPop->accelerator2() = 1.49618;
		newPop->setNeighborhood(m_random.get());
		rf = newPop->evaluate(m_problem.get(), this);
		newPop->initPbest(m_problem.get());
		informChange(rf);
		//cout << "Parent: " << newPop->neighborhoodBest(0,m_problem.get()).objective()[0] << endl;
		m_sub_pop.append(newPop);
		*m_preiter = m_sub_pop[0][m_sub_pop[0].best_idx(m_problem.get())].pbest();
	}

	void HmSO::initializeChild(HmSOSwarm& parent, int num) {
		int rf;
		int count = 0;
		std::unique_ptr<HmSOSwarm> newPop(new HmSOSwarm(num, m_problem.get(), HmSOSwarm::SwarmType::Child));
		newPop->weight() = 0.729844;
		newPop->accelerator1() = 1.49618;
		newPop->accelerator2() = 1.49618;
		newPop->setNeighborhood(m_random.get());
		(*newPop)[0] = parent[parent.best_idx(m_problem.get())];
		count++;
		for (size_t i(0); i < parent.size(); i++) {
			if (parent[i].pbest().variableDistance(parent[parent.best_idx(m_problem.get())].pbest(), m_problem.get()) < m_rpc/* && i != parent.best_idx(m_problem.get())*/) {
				if (count < num) {
					(*newPop)[count++] = parent[i];
				}
				parent[i].initialize(i, m_problem.get(), m_random.get());
				rf = parent[i].evaluate(m_problem.get(), this);
				informChange(rf);
				if (rf != kNormalEval) break;
				parent[i].pbest() = parent[i];
			}
		}
		while (count < num) {
			Particle add(m_problem->numberObjectives(), m_problem->numberConstraints(), CAST_CONOP(m_problem.get())->numberVariables());
			//add.initialize(count, m_problem.get(), m_random.get());
			for (size_t j = 0; j < add.variable().size(); j++) {
				Real temp = m_random->uniform.nextNonStd<Real>(-1, 1);
				add.variable()[j] = newPop->neighborhoodBest(0, m_problem.get()).variable()[j] + (m_rpc / 3.) * temp;
				add.velocity()[j] = m_random->uniform.nextNonStd<Real>(-1,1) * (m_rpc / 3.);
			}
			rf = add.evaluate(m_problem.get(), this); informChange(rf);
			add.pbest() = add;
			add.pbest() = add;
			(*newPop)[count++] = add;
		}
		//rf = newPop->evaluate(m_problem.get(), this);
		newPop->initPbest(m_problem.get());
		m_sub_pop.append(newPop);
		//informChange(rf);
	}

	void HmSO::checkOverlapping() {
		for (size_t i = 0; i < m_sub_pop.size(); i++) {
			for (size_t j = i + 1; j < m_sub_pop.size(); j++) {
				if (m_sub_pop[i].getswarmType() == HmSOSwarm::SwarmType::Child && m_sub_pop[j].getswarmType() == HmSOSwarm::SwarmType::Child
					&& m_sub_pop[i].neighborhoodBest(0, m_problem.get()).variableDistance(m_sub_pop[j].neighborhoodBest(0, m_problem.get()), m_problem.get()) < m_rexcl) {
					if (m_sub_pop[i].neighborhoodBest(0, m_problem.get()).dominate(m_sub_pop[j].neighborhoodBest(0, m_problem.get()), m_problem.get())) {
						m_sub_pop[j].setexcelFlag(true);
					}
					else m_sub_pop[i].setexcelFlag(true);
				}
			}
		}
	}

	void HmSO::removeOverlapping() {
		for (size_t i = 0; i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i].getswarmType() == HmSOSwarm::SwarmType::Child && m_sub_pop[i].getexcelFlag()) {
				m_sub_pop.remove(m_sub_pop.begin() + i);
				i--;
			}
		}
	}

	void HmSO::createChildSwarm(){
		if ((*m_curiter).dominate((*m_preiter), m_problem.get())) 
			initializeChild(m_sub_pop[0], 10);
		*m_preiter = *m_curiter;

	}

	int HmSO::updateBestSubPopIdx() {
		if (m_sub_pop.size() == 0) return -1;
		int idx = 0;
		for (int i(1); i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i].neighborhoodBest(0,m_problem.get()).dominate(m_sub_pop[idx].neighborhoodBest(0, m_problem.get()), m_problem.get())) {
				idx = i;
			}
		}
		//std::cout << "best subpop: " << idx << std::endl;
		return idx;
	}

	void HmSO::updateChildSwarm() {
		const int best_subpop_idx = updateBestSubPopIdx();
		for (auto& i : m_sub_pop)
		{
			if (i->getswarmType() == HmSOSwarm::SwarmType::Child && !i->gethiberFlag()) {
				int rf = kNormalEval;
				rf = i->evolve(m_problem.get(), this, m_random.get());
				informChange(rf);
				i->updateCurRadius(m_problem.get());
				//cout << m_sub_pop[i].getCurRadius() << '\t';
				if (i->getCurRadius() < m_rconv
					&& (*i)[i->best_idx(m_problem.get())].pbest().objective()[0] < (m_sub_pop[best_subpop_idx][m_sub_pop[best_subpop_idx].best_idx(m_problem.get())].pbest().objective()[0] - m_margin))
					i->sethiberFlag(true);
			}
			//cout << endl;
		}
		checkOverlapping();
		removeOverlapping();
	}

	void HmSO::updateChildAttractor(){
		int rf = kNormalEval;
		rf = m_sub_pop[0].evolve(m_problem.get(), this, m_random.get());
		informChange(rf);
		*m_curiter = m_sub_pop[0][m_sub_pop[0].best_idx(m_problem.get())].pbest();
		for (size_t i = 0; i < m_sub_pop[0].size(); i++) {
			for (size_t j = 1; j < m_sub_pop.size(); j++) {
				if (m_sub_pop[0][i].pbest().variableDistance(m_sub_pop[j][m_sub_pop[j].best_idx(m_problem.get())].pbest(), m_problem.get()) < m_rpc) {
					if (m_sub_pop[0][i].pbest().dominate(m_sub_pop[j][m_sub_pop[j].best_idx(m_problem.get())].pbest(), m_problem.get())) {
						m_sub_pop[j][m_sub_pop[j].best_idx(m_problem.get())] = m_sub_pop[0][i];
					}
					m_sub_pop[0][i].initialize(i, m_problem.get(), m_random.get());
					rf = m_sub_pop[0][i].evaluate(m_problem.get(), this);
					informChange(rf);
					m_sub_pop[0][i].pbest() = m_sub_pop[0][i];
				}
			}
		}
		createChildSwarm();
	}

	void HmSO::wakeupHiberPop() {
		for (size_t i(0); i < m_sub_pop.size(); i++) {
			if (m_sub_pop[i].getswarmType() == HmSOSwarm::SwarmType::Child && m_sub_pop[i].gethiberFlag()) {
				m_sub_pop[i].sethiberFlag(false);
			}
			if (m_sub_pop[i].getswarmType() == HmSOSwarm::SwarmType::Child) {
				relaunchChildSwarm(m_sub_pop[i]);
			}
		}
	}

	void HmSO::relaunchChildSwarm(HmSOSwarm& child) {
		const int best_idx = child.best_idx(m_problem.get());
		for (size_t i(0); i < child.size(); i++) {
			if (i != best_idx) {
				for (size_t j(0); j < child[i].variable().size(); j++) {
					child[i].variable()[j] = child[best_idx].pbest().variable()[j] + m_random->uniform.nextNonStd<Real>(-1, 1) * m_rs;
					child[i].velocity()[j] = m_random->uniform.nextNonStd<Real>(-1, 1) * m_rs;

				}
				child[i].evaluate(m_problem.get(), this);
				//child[i].initVelocity(m_problem.get(), m_random.get());
				child[i].pbest() = child[i];
			}
		}
	}

	void HmSO::informChange(int rf)
	{
		if (rf & kChangeCurEval) {
			measureMultiPop(true);
			wakeupHiberPop();
		}
	}

	bool HmSO::checkParentConvergence()
	{
		bool flag = true;
		for (int j = 0; j < m_sub_pop[0].size() && flag; j++) {
			for (int k = j + 1; k < m_sub_pop[0].size() && flag; k++) {
				for (int d = 0; d < CAST_CONOP(m_problem.get())->numberVariables(); d++) {
					if (fabs(static_cast<double>(m_sub_pop[0][j].variable()[d] - m_sub_pop[0][k].variable()[d])) > m_rconv) {
						flag = false;
						break;
					}
				}
			}
		}
		return flag;
	}

	void HmSO::preventParentConvergence(){
		int rf = kNormalEval;
		if (checkParentConvergence()) {
			for (size_t i = 0; i < m_sub_pop[0].size(); i++) m_sub_pop[0][i].initialize(i, m_problem.get(), m_random.get());
			rf = m_sub_pop[0].evaluate(m_problem.get(), this);
			m_sub_pop[0].initPbest(m_problem.get());
			informChange(rf);
		}
	}

	void HmSO::run_() {
		initializeParent(5);
		while (!terminating()) {
#ifdef OFEC_DEMO
			m_gbestslocation.clear();
#endif
			updateChildAttractor();
			createChildSwarm();
			//preventParentConvergence();
			updateChildSwarm();
#ifdef OFEC_DEMO
			for (size_t i(0); i < m_sub_pop.size(); i++) m_gbestslocation.push_back(m_sub_pop[i][m_sub_pop[i].best_idx(m_problem.get())].pbest());
			updateBuffer();
#endif
		}
	} 

	void HmSO::measureMultiPop(bool flag)
	{
		for (int i(0); i < m_sub_pop.size(); ++i) {
			for (int j(0); j < m_sub_pop[i].size(); j++) {
				m_sub_pop[i][j].pbest().evaluate(m_problem.get(), this);
			}
		}
	}

	void HmSO::record() {
		std::vector<Real> entry;
		entry.push_back(m_offline_error);
		entry.push_back(m_best_before_change_error);
		//entry.push_back(m_evaluations);
		//entry.push_back(GET_DYNCONOP(m_problem.get())->getCount());
		//entry.push_back(GET_DYNCONOP(m_problem.get())->optima().objective(0)[0]);
		//entry.push_back(m_candidates.front()->objective(0));
		dynamic_cast<RecordVecRealDynamic*>(m_record.get())->record(this, entry);
	}

#ifdef OFEC_DEMO
	void HmSO::updateBuffer() {
		m_solution.clear();
		m_solution.resize(m_sub_pop.size());
		for (size_t k = 0; k < m_sub_pop.size(); ++k) {
			for (size_t i = 0; i < m_sub_pop[k].size(); ++i) {
				m_solution[k].push_back(&m_sub_pop[k][i]);
			}
		}
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}

	std::vector<bool> HmSO::getPopHiberState() const {
		std::vector<bool> hiber_state;
		for (size_t i(0); i < m_sub_pop.size(); i++) {
			hiber_state.push_back(m_sub_pop[i].gethiberFlag());
		}
		return hiber_state;
	}
#endif

}