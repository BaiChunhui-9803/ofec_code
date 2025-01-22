#include "peak2peak.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	double peak2peak(const Ref<const VectorXd> &x) {
		return x.maxCoeff() - x.minCoeff();
	}
}