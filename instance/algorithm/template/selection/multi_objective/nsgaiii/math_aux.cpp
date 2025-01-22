#include "math_aux.h"
#include <math.h>
#include <limits>
#include <algorithm>

namespace ofec {
	namespace math_aux {
		Real ASF(const std::vector<Real>& objs, const std::vector<Real>& weight) {
			Real max_ratio = -std::numeric_limits<Real>::max();
			for (size_t f = 0; f < objs.size(); f += 1) {
				Real w = weight[f] ? weight[f] : 0.00001;
				max_ratio = std::max(max_ratio, objs[f] / w);
			}
			return max_ratio;
		}

		void guassianElimination(std::vector<Real>* px, std::vector<std::vector<Real> > A, const std::vector<Real>& b) {
			std::vector<Real>& x = *px;
			const size_t N = A.size();
			for (size_t i = 0; i < N; i += 1)
				A[i].push_back(b[i]);
			for (size_t base = 0; base < N - 1; base += 1) {
				for (size_t target = base + 1; target < N; target += 1) {
					Real ratio = A[target][base] / A[base][base];
					for (size_t term = 0; term < A[base].size(); term += 1) {
						A[target][term] -= A[base][term] * ratio;
					}
				}
			}
			x.resize(N);
			for (int i = N - 1; i >= 0; i -= 1) {
				for (size_t known = i + 1; known < N; known += 1)
					A[i][N] -= A[i][known] * x[known];
				x[i] = A[i][N] / A[i][i];
			}
		}


		// ---------------------------------------------------------------------
		// PerpendicularDistance:
		//
		// Given a direction vector (w1, w2) and a point P(x1, y1),
		// we want to find a point Q(x2, y2) on the line connecting (0, 0)-(w1, w2)
		// such that (x1-x2, y1-y2) is perpendicular to (w1, w2).
		//
		// Since Q is on the line (0, 0)-(w1, w2), it should be (w1*k, w2*k).
		// (x1-w1*k, y1-w2*k).(w1, w2) = 0. (inner product)
		// => k(w1^2 + w2^2) = w1x1 + w2x2
		// => k = (w1x1 + w2x2)/(w1^2 +w2^2).
		//
		// After obtaining k, we have Q = (w1*k, w2*k) and the distance between P and Q.
		//
		// Code example:
		//    vector<Real> dir{1, 3}, point{5.5, 1.5};
		//    cout << PerpendicularDistance(dir, point) << endl;
		// ---------------------------------------------------------------------
		Real perpendicularDistance(const std::vector<Real>& direction, const std::vector<Real>& point) {
			Real numerator = 0, denominator = 0;
			for (size_t i = 0; i < direction.size(); i += 1) {
				numerator += direction[i] * point[i];
				denominator += direction[i] * direction[i];
			}
			Real k = numerator / denominator;
			Real d = 0;
			for (size_t i = 0; i < direction.size(); i += 1)
				d += (k * direction[i] - point[i]) * (k * direction[i] - point[i]);
			return sqrt(d);
		}
	}
}