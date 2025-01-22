#include "initialize_subpop.h"
#include "calc_normalized_dis.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	void initializeSubpop(const MatrixXd &archived_solutions, const VectorXd &NormTabDis, size_t NsubP,
		const RowVectorXd &D_X, const RowVectorXd &U_X, const Opt &opt, Random *rnd, MatrixXd &X,
		VectorXd &Smean, MatrixXd &C, VectorXd &default_norm_tab_dis)
	{
		auto d0 = opt.d0hat;					// Initial default value of the normalized taboo distance
		size_t D = D_X.size();
		RowVectorXd str = U_X - D_X;							// Initial scaling factors
		C = str.array().square().matrix().asDiagonal();	// Initial covariance matrix of subpopulations
		MatrixXd Cinv = str.array().pow(-2).matrix().asDiagonal();
		if (NormTabDis.size() > 0)
			d0 = prctile(NormTabDis, opt.InidhatPrctile);	// Set the default value of the normalized taboo distance based on the distribution of the current values
		auto Smean_ = opt.MaxIniStepSize;
		int m = 0, k = 0;
		X.resize(NsubP, D);
		auto unif = &rnd->uniform;
		while (m < NsubP) {
			if (k > 100) {	// Too many sample points were rejected
				k = 0;
				Smean_ /= opt.RedRatio;	// reduce the stepsize
			}
			X.row(m) = D_X + (U_X - D_X).cwiseProduct(randVecXd(D, unif).transpose()); // Sample random point as a candidate for a subpopulation center
			if (m > 0 || NormTabDis.size() > 0.5) {
				VectorXd dis;
				if (m > 0)
					dis = calcNormalizedDis(joinRow(archived_solutions, X.topRows(m)), X.row(m), str, Cinv, Smean_, opt);
				else
					dis = calcNormalizedDis(archived_solutions, X.row(m), str, Cinv, Smean_, opt);
				ArrayXd tmp(NormTabDis.size() + m, 1);
				tmp << NormTabDis.array() + d0 * opt.InitializeDisMult, (1 + opt.InitializeDisMult) *d0 *ArrayXd::Ones(m, 1);
				ArrayXb accept_it = dis.array() >= tmp;	// Acceptance criterion
				if (accept_it.minCoeff() == true) {	// Accept the sampled point as the center of the m-th subpopulation 
					m++;	// Try finding a center for the next subpopulation
					k = 0;
				}
				else	// Reject the sampled point 
					k++;	// Count the number of consecutive rejections
			}
			else	// there is no taboo point, accept this point as the center of the first subpopulation
				m++;
		}
		if (opt.iniC == TypeIniC::kEye) {	// Use the Identity matrix for the initial covariance matrix (Not the default choice)
			auto mean_str = str.array().pow(1.0 / D).prod();
			C = pow(mean_str, 2) * MatrixXd::Identity(D, D);
		}
		default_norm_tab_dis = d0 * VectorXd::Ones(NsubP);
		Smean = Smean_ * VectorXd::Ones(NsubP);
	}
}