#include "csiwdn_individual.h"
#include "../../../problem/realworld/csiwdn/csiwdn.h"
#include "../../../../utility/functional.h"

namespace ofec {
	IndCSIWDN& IndCSIWDN::operator=(IndCSIWDN& indi) {
		Solution::operator=(indi);
		m_pv = indi.m_pv;
		m_pu = indi.m_pu;

		return *this;
	}

	IndCSIWDN& IndCSIWDN::operator=(IndCSIWDN&& indi) {
		Solution::operator=(std::move(indi));
		m_pv = std::move(indi.m_pv);
		m_pu = std::move(indi.m_pu);

		return *this;
	}

	void IndCSIWDN::mutateFirstPart(const std::vector<std::vector<Real>> &prob, Environment *env, Random *rnd, const std::pair<int, int>& source_index) {
		m_pu.variable() = variable();
		m_pu.variable().isDetected() = false;

		mutateLocation(env,rnd, prob, source_index);
	}

	void  IndCSIWDN::mutateSecondPart(Environment *env, Real F, int jdx, const std::vector<std::vector<Real>> &prob,
		Solution<VarCSIWDN> *r1,
		Solution<VarCSIWDN> *r2,
		Solution<VarCSIWDN> *r3,
		Solution<VarCSIWDN> *r4,
		Solution<VarCSIWDN> *r5) {
	
		m_pv.variable() = variable();
		m_pv.variable().isDetected() = false;

		mutateMass(env,F, jdx, r1, r2, r3, r4, r5);
	}

	void  IndCSIWDN::mutateMass(Environment *env, Real F, int jdx, Solution<VarCSIWDN> *r1,
		Solution<VarCSIWDN> *r2,
		Solution<VarCSIWDN> *r3,
		Solution<VarCSIWDN> *r4,
		Solution<VarCSIWDN> *r5) {

		auto pro = env->problem();

		/// mutate multiplier
		float u_multi =  CAST_CSIWDN(pro)->maxMultiplier();
		float l_multi =  CAST_CSIWDN(pro)->minMultiplier();

		
		size_t pv_size = m_pv.variable().duration(jdx) /  CAST_CSIWDN(pro)->patternStep();
		if (m_pv.variable().duration(jdx) %  CAST_CSIWDN(pro)->patternStep() != 0)
			++pv_size;
		long time_step =  CAST_CSIWDN(pro)->patternStep();
		std::vector<Real> vec_r1(r1->variable().multiplier(jdx));
		std::vector<Real> vec_r2(r2->variable().multiplier(jdx));
		std::vector<Real> vec_r3(r3->variable().multiplier(jdx));
		vec_r1.resize(pv_size);
		vec_r2.resize(pv_size);
		vec_r3.resize(pv_size);
		

		std::vector<Real> vec_r4;
		std::vector<Real> vec_r5;
		if (r4&&r5) {
			vec_r4 = r4->variable().multiplier(jdx);
			vec_r5 = r5->variable().multiplier(jdx);
			vec_r4.resize(pv_size);
			vec_r5.resize(pv_size);
		}
		float temp;
		for (size_t i = 0; i < pv_size; ++i) {
			m_pv.variable().multiplier(jdx)[i] = (vec_r1[i]) + F * ((vec_r2[i]) - (vec_r3[i]));
			if (r4 && r5) m_pv.variable().multiplier(jdx)[i] += F * ((vec_r4[i]) - (vec_r5[i]));
			if ((m_pv.variable().multiplier(jdx)[i]) > u_multi) {
				m_pv.variable().multiplier(jdx)[i] = ((vec_r1[i]) + u_multi) / 2;
			}
			else if ((m_pv.variable().multiplier(jdx)[i]) < l_multi) {
				m_pv.variable().multiplier(jdx)[i] = ((vec_r1[i]) + l_multi) / 2;
			}
		}

	}

	void  IndCSIWDN::mutateLocation(Environment *env, Random *rnd,const std::vector<std::vector<Real>> &prob, const std::pair<int, int>& source_index) {
		auto pro = env->problem();
		// roulette wheel selection for node location
		std::vector<std::vector<Real>> roulette_node(CAST_CSIWDN(pro)->numberSource(), std::vector<Real>(CAST_CSIWDN(pro)->numberNode() + 1, 0.));
		int q = source_index.first, z = source_index.second;
		for (int j = q; j <= z; ++j) {
			for (size_t i = 0; i < prob[j].size(); ++i) {
				roulette_node[j][i + 1] = roulette_node[j][i] + prob[j][i];
			}
		}
		Real p;
		for (int j = q; j <= z; ++j) {
			for (size_t i = 0; i < roulette_node[j].size() - 1; ++i) {
				p = rnd->uniform.nextNonStd<Real>(0, roulette_node[j][roulette_node[j].size() - 1]);
				if (p >= roulette_node[j][i] && p < roulette_node[j][i + 1])
					m_pu.variable().index(j) = i + 1;
			}
		}
	}

	void  IndCSIWDN::recombine(Environment *env, Random *rnd, size_t jdx, Real CR) {
		m_pu = m_pv;
		size_t pv_size_multiplier = m_pv.variable().multiplier(jdx).size();
		std::vector<Real> mult_temp(variable().multiplier(jdx));
		mult_temp.push_back(variable().startTime(jdx));
		++pv_size_multiplier;
		mult_temp.resize(pv_size_multiplier);

		size_t I = rnd->uniform.nextNonStd<size_t>(0, pv_size_multiplier);
		for (size_t i = 0; i < pv_size_multiplier; ++i) {
			Real p = rnd->uniform.next();
			if (i == pv_size_multiplier - 1) {
				if (p <= CR || i == I) {
					m_pu.variable().startTime(jdx) = m_pv.variable().startTime(jdx);
				}
				else {
					m_pu.variable().startTime(jdx) = mult_temp[i];
				}
			}
			else {
				if (p <= CR || i == I) {
					m_pu.variable().multiplier(jdx)[i] = m_pv.variable().multiplier(jdx)[i];
				}
				else {
					m_pu.variable().multiplier(jdx)[i] = mult_temp[i];
				}
			}
		}

	}

	int  IndCSIWDN::select(Environment *env,  bool is_stable, const std::pair<int, int>& source_index) {

		auto pro = env->problem();
		
		int tag = m_pu.evaluate(env, true);

		//if ((is_stable && (m_pu_distance_fitness > m_distance_fitness)) || (!is_stable && m_pu.dominate(m_objectives, env->optimizeMode()))) {
		if ((is_stable && (m_pu_distance_fitness > m_distance_fitness)) || (!is_stable && Dominance::kDominant == objectiveCompare(m_pu.objective(), m_objectives, pro->optimizeMode()))) {
			int z = source_index.second;
			SolutionType temp(1, 1);
			if (z > 0) {
				temp.variable() = variable();
				temp.variable().getEpa(z) = m_pu.variable().getEpa(z);
				tag = temp.evaluate(env, true);

			}
			else
				temp = m_pu;

			//auto& best = Dominance::kDominant == objectiveCompare(a.objective(), b.objective(), pro->optimizeMode()) ? temp : m_pu;
			auto &best = dominate(temp, m_pu, pro->optimizeMode()) ? temp : m_pu;
			variable() = best.variable();
			m_objectives = best.objective();
			m_constraints = best.constraint();
			m_improved = true;
		}
		else {
			m_improved = false;
		}
		return tag;
	}

	bool IndCSIWDN::sameLocation(Environment *env, IndCSIWDN & indi, const std::pair<int, int>& source_index) {
		int z = source_index.second;
		return variable().index(z) == indi.variable().index(z);
	}

	void IndCSIWDN::initialize(Environment *env, Random *rnd ) {
		Solution::initialize(env, rnd);
		m_pu.variable() = variable();
		m_pu.objective() = objective();
		m_pu.constraint() = constraint();

		m_pv.variable() = variable();
		m_pv.objective() = objective();
		m_pv.constraint() = constraint();
	}

	Solution< VarCSIWDN> &  IndCSIWDN::trial() {
		return m_pu;
	}

	void  IndCSIWDN::coverFirstPart(VarCSIWDN& indi, Environment *env, const std::pair<int, int>& source_index) {
		int q = source_index.first, z = source_index.second;
		variable().isDetected() = false;
		for (size_t i = q; i <= z; i++) {
			std::strcpy(variable().location(i), indi.location(i));
			variable().slocation(i) = indi.slocation(i);
			variable().index(i) = indi.index(i);
		}
	}

	void  IndCSIWDN::coverSecondPart(VarCSIWDN& indi, Environment *env, const std::pair<int, int>& source_index) {
		int q = source_index.first, z = source_index.second;
		variable().isDetected() = false;
		for (size_t i = q; i <= z; i++) {
			variable().startTime(i) = indi.startTime(i);
			variable().duration(i) = indi.duration(i);
			variable().multiplier(i) = indi.multiplier(i);
		}
	}
}