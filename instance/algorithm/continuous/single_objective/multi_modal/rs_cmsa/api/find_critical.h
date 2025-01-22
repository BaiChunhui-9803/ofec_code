#ifndef OFEC_RS_CMSA_ES_FIND_CRITICAL_H
#define OFEC_RS_CMSA_ES_FIND_CRITICAL_H

#include "../../../../../../../utility/linear_algebra/eigen.h"
#include "struct.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	VectorXi findCritical(const MatrixXd &all_better, const VectorXd &all_better_NormTabDis,
		const Ref<const RowVectorXd, 0, InnerStride<>> &Xmean,
		const VectorXd &Str, const MatrixXd &Cinv, double Smean, double criticality_treshhold, TypeDisMetric distance_metric);
}

#endif //!OFEC_RS_CMSA_ES_FIND_CRITICAL_H