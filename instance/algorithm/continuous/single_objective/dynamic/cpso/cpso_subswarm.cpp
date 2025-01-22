#include "cpso_subswarm.h"

namespace ofec {

	CPSOSwarm::CPSOSwarm(size_t pop, Problem *pro) :Swarm(pop, pro){
		m_center.reset(new templateParticle(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables()));
	}

	void CPSOSwarm::initializeParameters(Problem *pro, Algorithm *alg, Random *rnd) {
		if(alg->name() == "CPSO"){
			m_maxW = 0.6;
			m_minW = 0.3;
		}
		else if(alg->name() == "CPSOR"){
			m_maxW = 0.6;
			m_minW = 0.6;
		}
		m_accelerator1 = 1.7;
		m_accelerator2 = 1.7;

		m_weight = m_maxW;
		if (pro->hasTag(ProblemTag::kDOP))
			m_u = GET_DOP(pro)->getFrequency();
		else
			m_u = alg->maximumEvalutions();
		m_converge_flag = false;
		setNeighborhood(rnd);
	}

	//void CPSOSwarm::initialize(Problem *pro, Random *rnd) {
	//	m_weight = 0.6;
	//	m_accelerator1 = 1.7;
	//	m_accelerator2 = 1.7;

	//	m_maxW = 0.6;
	//	m_minW = 0.3;
	//		m_u = GET_DYNCONOP(pro)->getFrequency();	
	//	setNeighborhood(rnd);
	//	m_converge_flag = false;
	//	for (size_t i(0); i < this->size(); i++) (*this)[i].initialize(i, pro, rnd);
	//	updateCurRadius(pro);
	//	initializeRadius();
	//}

	void CPSOSwarm::initializeRadius() {
		m_initradius = m_radius;
	}

	int CPSOSwarm::localSearch(Problem *pro, Algorithm *alg, Random *rnd, int pop_size) {
		int rf = kNormalEval;
		//generate a permutation of templateParticle index
		std::vector<int> rindex(this->m_individuals.size());
		std::iota(rindex.begin(), rindex.end(), 0);
		rnd->uniform.shuffle(rindex.begin(), rindex.end());
		//this->setNeighborhood(rnd);

		bool flag = false;
		for (int i = 0; i < this->m_individuals.size(); i++) {
			auto& x = neighborhoodBest(rindex[i], pro);

			this->m_individuals[rindex[i]]->nextVelocity(&x, m_weight, m_accelerator1, m_accelerator2, rnd);
			this->m_individuals[rindex[i]]->move();
			this->m_individuals[rindex[i]]->clampVelocity(pro,rnd);

			rf = this->m_individuals[rindex[i]]->evaluate(pro, alg);
			if (rf != kNormalEval) return rf;
			if (this->m_individuals[rindex[i]]->dominate(this->m_individuals[rindex[i]]->pbest(), pro->optimizeMode())) {
				this->m_individuals[rindex[i]]->pbest() = this->m_individuals[rindex[i]];
				rf = learnGbest(this->m_individuals[rindex[i]], this->m_individuals[best_idx(pro)], pro, alg);
				if (rf != kNormalEval) return rf;
				if (this->updateBest(this->m_individuals[rindex[i]], pro))
					flag = true;
			}
		}
		if (rf == kNormalEval) {
			m_flag_best_impr = flag;
			this->m_iteration++;
			Real ratio = 0.;
			m_evals = alg->evaluations() % m_u;
			ratio = (m_iteration) / ((m_u - m_evals) / static_cast<Real> (pop_size));
			if (ratio > 1.0) ratio = 1.0;
			m_weight = m_maxW - (m_maxW - m_minW) * ratio;
		}
		updateCurRadius(pro);
		return rf;
	}

	void CPSOSwarm::setNeighborhood(Random *rnd) {
		m_link.assign(this->m_individuals.size(), std::vector<bool>(this->m_individuals.size(), true));
	}

	int CPSOSwarm::learnGbest(const std::unique_ptr<templateParticle>& pi, std::unique_ptr<templateParticle>& gbest,
		Problem *pro, Algorithm *alg) const{
		int rf = kNormalEval;

		size_t num_obj = CAST_CONOP(pro)->numberObjectives();
		size_t num_con = CAST_CONOP(pro)->numberConstraints();
		size_t num_var = CAST_CONOP(pro)->numberVariables();
		for (size_t i(0); i < num_var; i++) {
			std::unique_ptr<templateParticle> temp_templateParticle = std::make_unique<templateParticle>(num_obj, num_con, num_var);
			*temp_templateParticle = *gbest;

			temp_templateParticle->pbest().variable()[i] = pi->variable()[i];
			rf = temp_templateParticle->pbest().evaluate(pro, alg);
			if (temp_templateParticle->pbest().dominate(gbest->pbest(), pro)) {
				gbest->pbest() = temp_templateParticle->pbest();
//				gbest = gbest->pbest();
			}
			if (rf != kNormalEval) return rf;
		}
		return rf;
	}

	int CPSOSwarm::best_idx(Problem *pro) {
		int l = 0;
		for (int i = 0; i < this->m_individuals.size(); ++i) {
			if (this->m_individuals[i]->pbest().dominate(this->m_individuals[l]->pbest(), pro)) {
				l = i;
			}
		}
		return l;
	}

	void CPSOSwarm::computeCenter(Problem *pro) {
		if (this->size() < 1) return;
		if (pro->hasTag(ProblemTag::kConOP)) {
			for (size_t i = 0; i < CAST_CONOP(pro)->numberVariables(); i++) {
				double x = 0.;
				for (int j = 0; j < this->size(); j++) {
					auto chr = this->m_individuals[j]->pbest().variable();
					x = chr[i] + x;
				}
				m_center->variable()[i] = x / this->size();
			}
			//m_center.evaluate(true, caller::Algorithm);
		}
	}

	void CPSOSwarm::updateCurRadius(Problem *pro) {
		m_radius = 0;
		computeCenter(pro);
		if (this->size() < 2) return;
		for (int j = 0; j < this->size(); j++) m_radius += this->m_individuals[j]->pbest().variableDistance(*m_center, pro);
		m_radius = m_radius / this->size();
	}

	void CPSOSwarm::checkOverCrowd(size_t max_sub_size, Problem *pro, Random *rnd){
		if (m_converge_flag) return;
		if (size() <= max_sub_size) return;
		this->sortPbest(pro);
		for(size_t i(max_sub_size); i < size();i++){
			m_individuals.erase(m_individuals.begin() + i);
			i--;
		}
		setNeighborhood(rnd);
	}

	void CPSOSwarm::sortPbest(Problem *pro)
	{
		templateParticle temp(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables());
		for (int i = 0; i < this->size() - 1; i++) {
			for (int j = this->size() - 1; j > 0; j--) {
				if (m_individuals[j]->pbest().dominate(m_individuals[j - 1]->pbest(), pro)) {
					temp = (*m_individuals[j]);
					(*m_individuals[j]) = (*m_individuals[j - 1]);
					(*m_individuals[j - 1]) = temp;
				}
			}
		}
	}

	void CPSOSwarm::merge(CPSOSwarm& pop, Problem *pro, Random *rnd)
	{
		const size_t pre_size = this->size();
		this->m_individuals.resize(pre_size + pop.size());
		setNeighborhood(rnd);
		for (size_t i = 0; i < pop.size(); i++) {
			std::swap(this->m_individuals[pre_size + i], pop.m_individuals[i]);
		}
		m_initradius = (m_initradius + pop.m_initradius) / 2;
		updateCurRadius(pro);
	}
}
