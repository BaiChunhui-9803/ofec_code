#ifndef HIBERNATED_ZONE_H
#define HIBERNATED_ZONE_H

#include "peak_domain.h"
#include <climits>

namespace ofec {
	class HibernatedZone {
	protected:
		std::vector<PeakDomain> m_found_optimums;
	public:
		HibernatedZone() = default;

		void add_optimum(const Solution<>& sol, Problem *pro, Algorithm *alg, Random *rnd) {
			m_found_optimums.emplace_back(sol, pro, alg, rnd);
		}

		PeakDomain& operator[](int idx) {
			return m_found_optimums[idx];
		}

		void get_nearest_optimum(const Solution<>& sol, int& idx, Real& minDis, Problem *pro) {
			idx = -1;
			minDis = std::numeric_limits<Real>::max();
			Real cur_dis(0);
			for (int i(0); i < m_found_optimums.size(); ++i) {
				cur_dis = m_found_optimums[i].variableDistance(sol, pro);
				if (minDis > cur_dis) {
					minDis = cur_dis;
					idx = i;
				}
			}
		}

		Real getDeratingFactor(Real dis, Real ref) {
			if (ref <= 0) return 0;
			double f = dis * dis * dis / (ref * ref * ref);
			f = exp(-f) * (1 - f);
			return f > 0 ? f : 0;;
		}


		bool derateFitness(Solution<>& x, Problem *pro, Algorithm *alg, Random *rnd) {
			if (m_found_optimums.empty()) return false;
			int idx(0);
			Real dis(0);
			get_nearest_optimum(x, idx, dis, pro);
			if (idx == -1)return false;
			double refdis = m_found_optimums[idx].getRefRadi(x, pro, alg, rnd);
			if (dis < refdis) {
				x.objective(0) = m_found_optimums[idx].getDeratingObj(pro, alg);
				return true;
			}
			return false;
		}
	};



}


#endif