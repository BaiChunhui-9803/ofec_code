#ifndef OFEC_RS_CMSA_ES_CALC_NORMALIZED_DIS_H
#define OFEC_RS_CMSA_ES_CALC_NORMALIZED_DIS_H

#include "../../../../../../../utility/linear_algebra/eigen.h"
#include "struct.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	VectorXd calcNormalizedDis(const Ref<const MatrixXd> &X, const Ref<const RowVectorXd, 0, InnerStride<>> &x,
		const RowVectorXd &str, const MatrixXd &Cinv, double S, const Opt &opt);
}
#endif //!OFEC_RS_CMSA_ES_CALC_NORMALIZED_DIS_H