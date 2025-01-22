#ifndef OFEC_RS_CMSA_ES_CALC_DIS_H
#define OFEC_RS_CMSA_ES_CALC_DIS_H

#include "../../../../../../../utility/linear_algebra/eigen.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	VectorXd calcDis(const Ref<const MatrixXd> &X, const Ref<const RowVectorXd, 0, InnerStride<>> &x);
}

#endif //!OFEC_RS_CMSA_ES_CALC_DIS_H