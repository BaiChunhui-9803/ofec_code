#ifndef PEAK_DOMAIN_H
#define PEAK_DOMAIN_H

#include "../../../../../../core/problem/solution.h"
#include "../../../../../../utility/linear_algebra/vector.h"

namespace ofec {
	class PeakDomain : public Solution<> {
	protected:
		using Solution = Solution<>;
		Real m_radius;
		std::vector<Vector> m_sample;
		Real m_start_radius;
		Vector m_vbest;
		Real m_derating_obj;
		Real m_step;
		bool m_ready;
		int m_order;

	public:
		PeakDomain(int dim, int numObj);
		PeakDomain(const Solution<>& s, Problem *pro, Algorithm *alg, Random *rnd, Real r = 0, Real step = 0);

		void setRadius(Real r);
		Real getRefRadi(const Solution<>& s, Problem *pro, Algorithm *alg, Random *rnd);
		Real getRadius();
		void setDeratingObj(int num);
		Real getDeratingObj(Problem *pro, Algorithm *alg);
		PeakDomain& operator=(const PeakDomain& rhs);
		PeakDomain(const PeakDomain& rhs);
		Vector& operator[](int i);
		int size();
		Vector& getBest();

	protected:
		void generateSample(Problem *pro, Algorithm *alg, Random *rnd);
		void creatOneSample(const Vector& vnor, Problem *pro, Algorithm *alg, Random *rnd);
		void binarySearch(Solution<>& s0, Solution<>& s1, Problem *pro, Algorithm *alg);

	};
}


#endif