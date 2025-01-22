#ifndef OFEC_RS_CMSA_ES_EVOLVE_SUBPOP_H
#define OFEC_RS_CMSA_ES_EVOLVE_SUBPOP_H

#include "../../../../../../../utility/linear_algebra/eigen.h"
#include "struct.h"
#include "../../../../../../../core/algorithm/population.h"
#include "../../../../../../../core/problem/solution.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	void evolveSubpop(const Ref<const MatrixXd> &sup_Xmean, const Ref<const VectorXd> &sup_default_NormTabDis,
		Ref<RowVectorXd, 0, InnerStride<>> Xmean, double &Smean, size_t lambda, int mu, MatrixXd &eltX, VectorXd &eltf,
		VectorXd &eltS, MatrixXd &eltZ, const RowVectorXd &D_X, const RowVectorXd &U_X, const MatrixXd &Rot,
		const VectorXd &Str, MatrixXd &C, VectorXd &bestf_NE, VectorXd &medf_NE, Environment *env, Random *rnd,
		const MatrixXd &Xarchive, const VectorXd &Farchive, const VectorXd &NormTabDis, const RowVectorXd &W,
		const Opt &opt, Ref<RowVectorXd, 0, InnerStride<>> bestXf, size_t &used_evals,
		Population<Solution<>> &subpop);
}

#endif //!OFEC_RS_CMSA_ES_EVOLVE_SUBPOP_H