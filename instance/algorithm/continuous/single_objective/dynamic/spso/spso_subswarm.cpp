#include "spso_subswarm.h"

namespace ofec {

	SPSOSwarm::SPSOSwarm(size_t pop, Problem *pro) :Swarm(pop, pro){
		
	}

	void SPSOSwarm::initialize(Problem *pro, Random *rnd){
		m_weight = 0.729844;
		m_accelerator1 = 1.49618;
		m_accelerator2 = 1.49618;
		setNeighborhood(rnd);
		Swarm<SPSOParticle>::initialize(pro,rnd);
		//this->evaluate(pro, alg);
		//this->initPbest();
	}

	void SPSOSwarm::findpopSeed(){}

	int SPSOSwarm::evolve(Problem *pro, Algorithm *alg, Random *rnd){
		int rf = kNormalEval;
		//generate a permutation of particle index
		std::vector<int> rindex(this->m_individuals.size());
		std::iota(rindex.begin(), rindex.end(), 0);
		rnd->uniform.shuffle(rindex.begin(), rindex.end());
		bool flag = false;
		for (int i = 0; i < this->m_individuals.size(); i++) {
			auto& x = neighborhoodBest(rindex[i], pro);
			this->m_individuals[rindex[i]]->nextVelocity(&x, m_weight, m_accelerator1, m_accelerator2, rnd);
			this->m_individuals[rindex[i]]->move();
			this->m_individuals[rindex[i]]->clampVelocity(pro,rnd);
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

	int SPSOSwarm::best_idx(Problem *pro) {
		int l = 0;
		for (int i = 0; i < this->m_individuals.size(); ++i) {
			if (this->m_individuals[i]->pbest().dominate(this->m_individuals[l]->pbest(), pro)) {
				l = i;
			}
		}
		return l;
	}

	void SPSOSwarm::setNeighborhood(Random *rnd){
		m_link.assign(this->m_individuals.size(), std::vector<bool>(this->m_individuals.size(), true));
	}
}