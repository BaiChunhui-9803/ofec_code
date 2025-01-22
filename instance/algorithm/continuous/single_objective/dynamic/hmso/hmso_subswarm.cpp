#include "hmso_subswarm.h"

namespace ofec {
	HmSOSwarm::HmSOSwarm(size_t pop, Problem *pro, SwarmType t) :Swarm(pop, pro), m_swarmtype(t){
		m_center.reset(new Solution<>(pro->numberObjectives(), pro->numberConstraints(), CAST_CONOP(pro)->numberVariables()));
	}

	//void HmSOSwarm::initialize()
	//{
	//	//set_parameters(0.729844, 1.49618, 1.49618);
	//	/*swarm<particle>::initialize();*/
	//	m_center.reset(new Solution<>(pro->numberObjectives(), pro->numberConstraints(), CAST_CONOP(pro)->numberVariables()));

	//	m_pop.reset(new Swarm<Particle07>(m_pop_size, m_problem.get()));
	//	Swarm::initialize(m_problem.get(), m_random.get());
	//	m_pop->initVelocity(m_problem.get(), m_random.get());
	//	m_pop->evaluate(m_problem.get(), this);
	//	m_pop->initPbest(m_problem.get());
	//	m_pop->W() = m_w;
	//	m_pop->C1() = m_c1;
	//	m_pop->C1() = m_c2;
	//	setNeighborhood();
	//}

	int HmSOSwarm::evolve(Problem *pro, Algorithm *alg, Random *rnd)
	{
		int rf = kNormalEval;
		//generate a permutation of particle index
		std::vector<int> rindex(this->m_individuals.size());
		std::iota(rindex.begin(), rindex.end(), 0);
		rnd->uniform.shuffle(rindex.begin(), rindex.end());

		//this->setNeighborhood(rnd);

		bool flag = false;
		for (int i = 0; i < this->m_individuals.size(); i++) {
			auto& x = neighborhoodBest(rindex[i], pro);

			this->m_individuals[rindex[i]]->nextVelocity(&x, m_weight, m_accelerator1, m_accelerator2, rnd);
			this->m_individuals[rindex[i]]->move();
			this->m_individuals[rindex[i]]->clampVelocity(pro, rnd);

			rf = this->m_individuals[rindex[i]]->evaluate(pro, alg);

			if (this->m_individuals[rindex[i]]->dominate(this->m_individuals[rindex[i]]->pbest(), pro->optimizeMode())) {
				this->m_individuals[rindex[i]]->pbest() = this->m_individuals[rindex[i]];
				if (this->updateBest(this->m_individuals[rindex[i]], pro))
					flag = true;
			}

			if (rf != kNormalEval) break;
		}
		m_flag_best_impr = flag;
		this->m_iteration++;
		return rf;
	}
	void HmSOSwarm::setNeighborhood(Random *rnd) {
		m_link.assign(this->m_individuals.size(), std::vector<bool>(this->m_individuals.size(), true));
	}

	void HmSOSwarm::computeCenter(Problem *pro) {
		if (this->size() < 1) return;
		if (pro->hasTag(ProblemTag::kConOP)) {
			for (size_t i = 0; i < CAST_CONOP(pro)->numberVariables(); i++) {
				double x = 0.;
				for (int j = 0; j < this->size(); j++) {
					auto chr = this->m_individuals[j]->variable();
					x = chr[i] + x;
				}
				m_center->variable()[i] = x / this->size();
			}
			//m_center.evaluate(true, caller::Algorithm);
		}
	}
	void HmSOSwarm::updateCurRadius(Problem *pro) {
		m_radius = 0;
		computeCenter(pro);
		if (this->size() < 2) return;
		for (int j = 0; j < this->size(); j++) m_radius += this->m_individuals[j]->variableDistance(*m_center, pro);
		m_radius = m_radius / this->size();
	}
}
