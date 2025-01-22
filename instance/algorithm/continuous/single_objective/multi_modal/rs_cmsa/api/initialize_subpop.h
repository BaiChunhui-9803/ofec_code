#ifndef OFEC_RS_CMSA_ES_INITIALIZE_SUBPOP_H
#define OFEC_RS_CMSA_ES_INITIALIZE_SUBPOP_H

#include "../../../../../../../utility/linear_algebra/eigen.h"
#include "struct.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	void initializeSubpop(const MatrixXd &archived_solutions, const VectorXd &NormTabDis, size_t NsubP,
		const RowVectorXd &D_X, const RowVectorXd &U_X, const Opt &opt, Random *rnd, MatrixXd &X,
		VectorXd &Smean, MatrixXd &C, VectorXd &default_NormTabDis);
}

#endif //!OFEC_RS_CMSA_ES_INITIALIZE_SUBPOP_H