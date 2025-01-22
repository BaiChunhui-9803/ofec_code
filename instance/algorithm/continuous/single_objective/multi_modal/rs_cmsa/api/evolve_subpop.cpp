#include "evolve_subpop.h"
#include "find_critical.h"
#include "calc_normalized_dis.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	void evolveSubpop(const Ref<const MatrixXd> &sup_Xmean, const Ref<const VectorXd> &sup_default_NormTabDis,
		Ref<RowVectorXd, 0, InnerStride<>> Xmean, double &Smean, size_t lambda, int mu, MatrixXd &eltX, VectorXd &eltf,
		VectorXd &eltS, MatrixXd &eltZ, const RowVectorXd &D_X, const RowVectorXd &U_X, const MatrixXd &Rot,
		const VectorXd &Str, MatrixXd &C, VectorXd &bestf_NE, VectorXd &medf_NE, Environment *env, Random *rnd,
		const MatrixXd &Xarchive, const VectorXd &Farchive, const VectorXd &NormTabDis, const RowVectorXd &W,
		const Opt &opt, Ref<RowVectorXd, 0, InnerStride<>> bestXf, size_t &used_evals,
		Population<Solution<>> &subpop)
	{
		auto D = D_X.size();	// Problem dimension
		auto eltNo = eltf.size();	// Number of elite Solutions
		MatrixXd Cinv = Rot * Str.array().square().matrix().asDiagonal() * Rot.transpose();	// Inverse of the covariance matrix
		used_evals = 0;
		double zaribR = 1;	// For temporary reduction of the taboo region
		MatrixXd all_better = sup_Xmean;	// Center of fitter subpopualations
		VectorXd all_better_NormTabDis = sup_default_NormTabDis;	// The normalized taboo distance of fitter subpopulations
		for (size_t o1 = 0; o1 < NormTabDis.size(); ++o1) {	// Consider fitter archived points as taboo points
			if ((Farchive(o1) < eltf.array() + opt.ArchiveConsiderTabooTol).minCoeff()) {	// Margin of superiority 
				all_better = joinRow(all_better, Xarchive.row(o1));
				all_better_NormTabDis = join(all_better_NormTabDis, NormTabDis(o1));
			}
		}
		// Now determine which taboo points are critical
		VectorXi sorted_critical_taboo_ind = findCritical(all_better, all_better_NormTabDis, Xmean, Str, Cinv,
			Smean, opt.CriticThresh, opt.DisMetric);
		auto No_critical_taboo = sorted_critical_taboo_ind.size();
		MatrixXd X(lambda, D); MatrixXd Z = X; VectorXd S(lambda); VectorXd f = S.array() + 1e100;
		//Solution<> tmp_sol(1, 0, D);
		OptimizeMode opt_mode = env->problem()->optimizeMode(0);
		auto norm = &rnd->normal;
		for (size_t k = 0; k < lambda; ++k) {
			bool accept_it = false;
			while (!accept_it) {
				S(k) = Smean * exp(norm->next() * (sqrt(0.5 * opt.TSigmaCoeff / D)));	// Sampling
				Z.row(k) = Rot * Str.cwiseProduct(randnVecXd(D, norm));
				X.row(k) = Xmean + S(k) * Z.row(k);
				if (No_critical_taboo == 0)	// There is no critical taboo point
					accept_it = true;
				else {
					for (int TabNo : sorted_critical_taboo_ind) {	// Check against all critical taboo points
						auto dis_ratio = calcNormalizedDis(all_better.row(TabNo), X.row(k), Str, Cinv, Smean, opt)(0)
							/ all_better_NormTabDis(TabNo);
						accept_it = dis_ratio > zaribR;
						if (accept_it == false) {	// If the sample was rejected
							zaribR = zaribR / (opt.RedRatio);	// Temporary shrink the size of the taboo regions
							break;	// Do not check wrt other critical taboo points
						}
					}
				}
			}	// A taboo acceptable solution was generated
			if ((X.row(k).array() > U_X.array()).maxCoeff() || (X.row(k).array() < D_X.array()).maxCoeff()) {	// Out of the search range
				auto pen1 = -(X.row(k).cwiseAbs() - D_X);
				auto pen2 = -(U_X - X.row(k).cwiseAbs());
				f(k) = 1e100 * (1 + ((pen1 + pen1.cwiseAbs()).array().pow(2) + (pen2 + pen2.cwiseAbs()).array().pow(2)).sum());	// Death penalty : No function
			}
			else {
				for (size_t j = 0; j < D; ++j) subpop[k].variable()[j] = X(k, j);
				subpop[k].evaluate(env);
				f(k) = opt_mode == OptimizeMode::kMinimize ? subpop[k].objective(0) : -subpop[k].objective(0);
				used_evals++;
			}
		}	// lambda taboo acceptable solutions were generated and evaluated
		// update meassures for termination
		bestf_NE = join(bestf_NE, f.minCoeff());
		medf_NE = join(medf_NE, median(f));
		for (size_t k0 = 0; k0 < eltNo; ++k0) {	// Add elite solutions to the recently generated ones
			if (eltf(k0) < 1e100) {	// Only if they are feasible and are not in the taboo regions
				if (No_critical_taboo == 0) {	// There is no taboo point
					S = join(S, eltS(k0));
					Z = joinRow(Z, eltZ.row(k0));
					X = joinRow(X, eltX.row(k0));
					f = join(f, eltf(k0));
				}
				else {	// There are taboo points. Check the elite solution wrt them
					bool accept_it;
					for (auto TabNo : sorted_critical_taboo_ind) {	// Should not be in the taboo regions
						auto dis_ratio = calcNormalizedDis(all_better.row(TabNo), eltX.row(k0), Str, Cinv, Smean, opt)(0)
							/ all_better_NormTabDis(TabNo);
						accept_it = (dis_ratio > zaribR);
						if (accept_it == false)
							break;
					}
					if (accept_it == true) {	// The elite was not in the a taboo region
						S = join(S, eltS(k0));
						Z = joinRow(Z, eltZ.row(k0));
						X = joinRow(X, eltX.row(k0));
						f = join(f, eltf(k0));
					}
				}
			}
		}
		// Recombine and update parameters of the subpopulation
		VectorXi ind = sort(f);
		Xmean = W * X(ind.head(mu), indexing::all);	// Recombination of design variables
		Smean = exp(W * S(ind.head(mu)).array().log().matrix()) / exp(S.array().log().mean()) * Smean;	// Bias compensation for recombination of the global step size
		double Tc = 1 + D * (D + 1) / mu * opt.TCovCoeff;
		MatrixXd suggC = MatrixXd::Zero(D, D);
		for (size_t k0 = 0; k0 < mu; ++k0)
			suggC += W(k0) * Z.row(ind(k0)).transpose() * Z.row(ind(k0));
		C = (1 - 1 / Tc) * C + suggC / Tc;
		C = (C + C.transpose()) / 2;
		// Update the elite solutions
		eltX = X(ind.head(eltNo), indexing::all);
		eltf = f(ind.head(eltNo));
		eltS = S(ind.head(eltNo));
		eltZ = Z(ind.head(eltNo), indexing::all);
		Xmean = U_X.cwiseMin(D_X.cwiseMax(Xmean));
		bestXf = join(X.row(ind(0)), f(ind(0)));
	}
}