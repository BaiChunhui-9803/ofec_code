#ifndef OFEC_RS_CMSA_ES_RUN_RS_CMSA_H
#define OFEC_RS_CMSA_ES_RUN_RS_CMSA_H

#include "../../../../../../../utility/linear_algebra/eigen.h"
#include "struct.h"
#include "../../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../../core/algorithm/population.h"
#include "../../../../../../../core/problem/solution.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	void runRSCMSA(size_t lambda, size_t NsubP, size_t D, const RowVectorXd &D_X, const RowVectorXd &U_X,
		size_t &remained_eval, MatrixXd &all_converged, MatrixXd &archiveXf, VectorXd &all_N_rep,
		VectorXd &NormTabDis, const Opt &opt, Environment *env, Random *rnd, double &average_used_iter,
		MultiPopulation<Population<Solution<>>> &subpops);
}

#endif // !OFEC_RS_CMSA_ES_RUN_RS_CMSA_H
