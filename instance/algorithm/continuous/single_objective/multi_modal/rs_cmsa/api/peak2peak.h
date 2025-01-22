#ifndef OFEC_RS_CMA_ES_PEAK2PEAK_H
#define OFEC_RS_CMA_ES_PEAK2PEAK_H

#include "../../../../../../../utility/linear_algebra/eigen.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	double peak2peak(const Ref<const VectorXd> &x);
}

#endif //!OFEC_RS_CMA_ES_PEAK2PEAK_H