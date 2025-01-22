#include "peak_domain.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	PeakDomain::PeakDomain(int num_vars, int number_objectives) :
		Solution(num_vars, number_objectives),
		m_radius(0),
		m_start_radius(0),
		m_vbest(num_vars),
		m_step(0),
		m_ready(false),
		m_order(1),
		m_derating_obj(0) {}

	PeakDomain::PeakDomain(const Solution& s, Problem *pro, Algorithm *alg, Random *rnd, Real startRadi, Real step) :
		Solution<>(s),
		m_radius(0),
		m_start_radius(startRadi),
		m_vbest(s.variable().vect()),
		m_step(step),
		m_ready(false),
		m_order(1)
	{
		m_derating_obj = pro->optimizeMode(0) == OptimizeMode::kMinimize ?
			alg->maxObjFound(0) :
			alg->minObjFound(0);
		if (m_start_radius == 0) {
			m_start_radius = CAST_CONOP(pro)->variableAccuracy() / 3;
		}
		if (m_step == 0) {
			m_step = m_start_radius / 3;
		}
		generateSample(pro, alg, rnd);
	}

	void PeakDomain::setRadius(Real r) {
		m_radius = r;
	}

	Real PeakDomain::getRadius() {
		return m_radius;
	}

	void PeakDomain::setDeratingObj(int num) {
		m_order = num;
	}

	Real PeakDomain::getDeratingObj(Problem *pro, Algorithm *alg) {
		Real w = pow(1.001, -sqrt(m_order));
		if (pro->optimizeMode(0) == OptimizeMode::kMinimize)
			m_derating_obj = alg->minObjFound(0) * (1 - w) + alg->maxObjFound(0) * w;
		else
			m_derating_obj = alg->maxObjFound(0) * (1 - w) + alg->minObjFound(0) * w;
		return m_derating_obj;
	}

	Real PeakDomain::getRefRadi(const Solution<>& s, Problem *pro, Algorithm *alg, Random *rnd) {
		if (!m_ready) return 0;
		Vector v(s.variable().vect());
		v -= m_vbest;
		if (v.length() == 0) return m_radius;
		std::vector<Real> angle;
		for (auto& i : m_sample) {
			angle.push_back(v.angle(i));
		}
		int idx = min_element(angle.begin(), angle.end()) - angle.begin();
		Real mang = angle.empty() ? OFEC_PI : angle[idx];
		if (mang > 3 * OFEC_PI / 180) {
			m_ready = false;
			creatOneSample(v, pro, alg, rnd);
			m_ready = true;
			return m_sample.back().length();
		}
		return m_sample[idx].length();
	}

	void PeakDomain::generateSample(Problem *pro, Algorithm *alg, Random *rnd) {
		int dim = variable().size();
		for (int tr = 0; tr < dim * 2; ++tr) {
			Vector vnor(dim);
			vnor.randomize(&rnd->uniform, -1, 1);
			vnor.normalize();
			creatOneSample(vnor, pro, alg, rnd);
		}
		m_ready = true;
	}

	void PeakDomain::creatOneSample(const Vector& vnor, Problem *pro, Algorithm *alg, Random *rnd) {
		Vector vtr(vnor);
		Solution<> s0(*this), s1(*this);
		Real r = m_start_radius;
		bool flag = true, first = true;
		Real dr0 = 1, dr1 = 0;
		do {
			s0 = s1;
			vtr = vnor * r;
			vtr += m_vbest;
			copy(vtr.begin(), vtr.end(), s1.variable().begin());
			if (!s1.boundaryViolated(pro)) {
				s1.evaluate(pro, alg);
			}
			else {
				pro->validateSolution(s1, Validation::kSetToBound, rnd);
				//	SolutionValidation mode = VALIDATION_SETTOBOUND;
				//	s1.validate(&mode);
				s0 = s1;
				flag = false;
				break;
			}
			//June 19, 2016 TO DO: need to test
			/*if (dynamic_cast<FreePeak*>(Global::msp_global->mp_problem.get())) {
				dr1 = s0.getObjDistance_(s1.data().m_objectives);
				if (!first) {
					if (dr1 / dr0 > 10) {
						flag = false;
						break;
					}
				}
				dr0 = dr1;
				if (first) first = false;
			}*/
			//***********
			r += m_step;
		} while (s1.dominate(s0, pro->optimizeMode()));
		if (flag) {
			binarySearch(s0, s1, pro, alg);
		}
		Vector vr(s0.variable().vect());
		vr -= m_vbest;
		m_sample.push_back(std::move(vr));
		if (m_sample.size() == 1 || m_radius > r) m_radius = r;
	}

	PeakDomain& PeakDomain::operator=(const PeakDomain& rhs) {
		Solution::operator=(rhs);
		m_radius = rhs.m_radius;
		m_sample = rhs.m_sample;
		m_start_radius = rhs.m_start_radius;
		m_vbest = rhs.m_vbest;
		m_derating_obj = rhs.m_derating_obj;
		m_step = rhs.m_step;
		m_ready = rhs.m_ready;
		m_order = rhs.m_order;
		return *this;
	}

	PeakDomain::PeakDomain(const PeakDomain& rhs) :Solution(rhs) {
		m_radius = rhs.m_radius;
		m_sample = rhs.m_sample;
		m_start_radius = rhs.m_start_radius;
		m_vbest = rhs.m_vbest;
		m_derating_obj = rhs.m_derating_obj;
		m_step = rhs.m_step;
		m_ready = rhs.m_ready;
		m_order = rhs.m_order;
	}

	Vector& PeakDomain::operator[](int i) {
		return m_sample[i];
	}

	int PeakDomain::size() {
		return m_sample.size();
	}

	Vector& PeakDomain::getBest() {
		return m_vbest;
	}

	void PeakDomain::binarySearch(Solution<>& s0, Solution<>& s1, Problem *pro, Algorithm *alg) {
		Solution<> s(s0);
		for (int i = 0; i < s.variable().size(); ++i) {
			s.variable()[i] = (s0.variable()[i] + s1.variable()[i]) / 2.;
		}
		s.evaluate(pro, alg);
		if (s0.variableDistance(s, pro) <= 1.e-3) {
			s0 = s;
			return;
		}
		if (s0.dominate(s, pro->optimizeMode()))
			binarySearch(s, s1, pro, alg);
		if (s.dominate(s0, pro->optimizeMode()))
			binarySearch(s0, s, pro, alg);
	}
}
