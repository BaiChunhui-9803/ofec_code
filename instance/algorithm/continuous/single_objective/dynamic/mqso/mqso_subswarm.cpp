#include "mqso_subswarm.h"

#ifdef __GNUC__
#include <float.h>
#endif

namespace ofec {
	MQSOSwarm::MQSOSwarm(size_t pop, Problem *pro) :Swarm(pop,pro){}

	void MQSOSwarm::initialize(Algorithm *alg, Problem *pro, Random *rnd)
	{
		m_weight = 0.729844;
		m_accelerator1 = 1.49618;
		m_accelerator2 = 1.49618;
		m_Nplus = 5, m_Nq = 5, m_N = 5, m_Q = (pow(1. / 4.9, 1. / 0.6)), m_Rcloud = (0.5);	//parameters in swarm.
		m_center.reset(new Solution<>(CAST_CONOP(pro)->numberObjectives(), CAST_CONOP(pro)->numberConstraints(), CAST_CONOP(pro)->numberVariables()));
		assignType(alg);
		m_swarmType = SwarmType::SWARM_FREE;
		for (auto& i : this->m_individuals) i->initializeVmax(pro);
		setNeighborhood(rnd);
		for (size_t i(0); i < m_individuals.size(); i++) 
			(*this)[i].initialize(i, pro,rnd);
	}

	int MQSOSwarm::reInitialize(Algorithm *alg, Problem *pro, Random *rnd){
		//m_swarmType = SwarmType::SWARM_FREE;
		int rf = kNormalEval;
		for (size_t i(0); i < m_individuals.size(); i++)
			m_individuals[i]->initialize(i, pro, rnd);
		rf = evaluate(pro, alg);
		if (rf != kNormalEval) return rf;
		initVelocity(pro, rnd);
		initPbest(pro);
		updateCurRadius(pro);
		m_worstFlag = false;
		return rf;
	}

	void MQSOSwarm::assignType(Algorithm *alg) {
		if (alg->name() == "mCPSO") {
			for (int i = 0; i < size(); i++) {
				if (i < m_N) m_individuals[i]->getType() = MQSOParticle::ParticleType::PARTICLE_NEUTRAL;
				else  m_individuals[i]->getType() = MQSOParticle::ParticleType::PARTICLE_CHARGED;
			}
		}
		else if (alg->name() == "mQSO" || alg->name() == "SAMO") {
			for (int i = 0; i < size(); i++) {
				if (i < m_N) m_individuals[i]->getType() = MQSOParticle::ParticleType::PARTICLE_NEUTRAL;
				else  m_individuals[i]->getType() = MQSOParticle::ParticleType::PARTICLE_QUANTUM;
			}
		}
	}

	void MQSOSwarm::setswarmType(SwarmType r)
	{
		m_swarmType = r;
	}

	int MQSOSwarm::bestIdx(Problem *pro)
	{
		int l = 0;
		for (int i = 0; i < this->m_individuals.size(); ++i) {
			if (this->m_individuals[i]->pbest().dominate(this->m_individuals[l]->pbest(),pro)) {
				l = i;
			}
		}
		return l;
	}

	void MQSOSwarm::computeRepulsion(const int idx,Problem *pro) {
		int numVar = CAST_CONOP(pro)->numberVariables();
		for (int d = 0; d < numVar; d++) {
			m_individuals[idx]->getRepulse()[d] = 0;
		}
		for (int i = 0; i < size(); i++) {
			if (i != idx && m_individuals[i]->getType() == MQSOParticle::ParticleType::PARTICLE_CHARGED) {
				double m = 0, x, y;
				for (int d = 0; d < numVar; d++) {
					x = m_individuals[idx]->variable()[d];
					y = m_individuals[i]->variable()[d];
					m += (x - y) * (x - y);
				}
				m = sqrt(m);
				m = pow(m, 3.);
				for (int d = 0; d < numVar; d++) {
					x = m_individuals[idx]->variable()[d];
					y = m_individuals[i]->variable()[d];
					m_individuals[idx]->getRepulse()[d] += (x - y) * m_Q * m_Q / m;
				}
			}
		}

		// clamp to max<double>
		long double m = 0;
		for (int d = 0; d < numVar; d++) {
			m += m_individuals[idx]->getRepulse()[d] * m_individuals[idx]->getRepulse()[d];
		}
		m = sqrt(m);
#ifdef __GNUC__
		if (m > DBL_MAX) {
#elif defined _MSC_VER
		if (m > DBL_MAX) {
#endif
			for (int d = 0; d < numVar; d++) {
				m_individuals[idx]->getRepulse()[d] = m_individuals[idx]->getRepulse()[d] * DBL_MAX / m;
			}
		}
		for (int d = 0; d < numVar; d++) {
			if (std::isnan(m_individuals[idx]->getRepulse()[d])) m_individuals[idx]->getRepulse()[d] = 1;	//if calculate repulse too big(to NAN), let it be 1;
		}
	}

	void MQSOSwarm::computeCenter(Problem *pro) {
		if (this->size() < 1) return;
		int neutral_count = 0;
		if (pro->hasTag(ProblemTag::kConOP)) {
			for (size_t i = 0; i < CAST_CONOP(pro)->numberVariables(); i++) {
				double x = 0.;
				for (int j = 0; j < this->size(); j++) {
					if (this->m_individuals[j]->getType() == MQSOParticle::ParticleType::PARTICLE_NEUTRAL) {
						neutral_count++;
						auto chr = this->m_individuals[j]->pbest().variable();
						x = chr[i] + x;
					}
				}
				m_center->variable()[i] = x / neutral_count;
				neutral_count = 0;
			}
			//m_center.evaluate(true, caller::Algorithm);
		}
	}

	void MQSOSwarm::updateCurRadius(Problem *pro)
	{
		m_radius = 0;
		int neutral_count = 0;
		computeCenter(pro);
		if (this->size() < 2) return;
		for (int j = 0; j < this->size(); j++) {
			if (this->m_individuals[j]->getType() == MQSOParticle::ParticleType::PARTICLE_NEUTRAL) {
				m_radius += this->m_individuals[j]->pbest().variableDistance(*m_center, pro);
				neutral_count++;
			}
		}
		m_radius = m_radius / neutral_count;
		neutral_count = 0;
	}

	int MQSOSwarm::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		int rf = kNormalEval;
		if (this->size() < 1) return rf = kNormalEval;
		//generate a permutation of particle index
		std::vector<int> rindex(this->m_individuals.size());
		std::iota(rindex.begin(), rindex.end(), 0);
		rnd->uniform.shuffle(rindex.begin(), rindex.end());
		//this->setNeighborhood(rnd);

		bool flag = false;
		for (int i = 0; i < this->m_individuals.size(); i++) {
			auto& x = neighborhoodBest(rindex[i], pro);
			if (m_individuals[rindex[i]]->getType() == MQSOParticle::ParticleType::PARTICLE_CHARGED) computeRepulsion(rindex[i], pro);
			if (m_individuals[rindex[i]]->getType() != MQSOParticle::ParticleType::PARTICLE_QUANTUM) {
				this->m_individuals[rindex[i]]->nextVelocity(&x, m_weight, m_accelerator1, m_accelerator2, rnd);
				this->m_individuals[rindex[i]]->move();
			}
			if (m_individuals[rindex[i]]->getType() == MQSOParticle::ParticleType::PARTICLE_QUANTUM) {
				m_individuals[rindex[i]]->quantumMove(&x, m_Rcloud, rnd);
			}
			this->m_individuals[rindex[i]]->clampVelocity(pro, rnd);
			rf = this->m_individuals[rindex[i]]->evaluate(pro, alg);

			//rf = m_individuals[i]->moveto(x, m_weight, m_C1, m_C2, m_Rcloud);
			//if (rf != kNormalEval) break;
			if (this->m_individuals[rindex[i]]->dominate(this->m_individuals[rindex[i]]->pbest(),pro)) {
				this->m_individuals[rindex[i]]->pbest() = this->m_individuals[rindex[i]];
				if (this->updateBest(this->m_individuals[rindex[i]],pro))
					flag = true;
			}
			if (rf != kNormalEval) break;
			//cout << m_individuals[i]->variable()[0] << " & "<< m_individuals[i]->variable()[1] << endl;
		}
		m_flag_best_impr = flag;
		this->m_iteration++;
		return rf;
	}

	void MQSOSwarm::setNeighborhood(Random *rnd) {
		m_link.assign(this->m_individuals.size(), std::vector<bool>(this->m_individuals.size(), true));
	}
}