#include "analyze_converged_sol.h"
#include "calc_dis.h"
#include "is_new_basin.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	void analyzeConvergedSol(const Ref<const MatrixXd> &Y, Environment *env, MatrixXd &archiveXf, VectorXd &all_N_rep,
		VectorXd &NormTabDis, const Opt &opt, size_t remained_eval, size_t &used_eval)
	{
		int D = Y.cols() - 1;
		VectorXd old_NormTabDis = NormTabDis;
		VectorXd N_rep = NormTabDis * 0;	// Consider only repitition in the recent restart
		used_eval = 0;
		double d0hat = opt.d0hat;
		VectorXi globY_ind;
		if (archiveXf.rows() > 0) {
			if (Y.col(D).minCoeff() + opt.TargetTolFun < archiveXf.col(D).minCoeff()) {	// The archived solutions were not global minimum, discard them
				archiveXf.resize(0, NoChange); NormTabDis.resize(0); N_rep.resize(0); all_N_rep.resize(0);
				globY_ind = find(Y.col(D).array() <= opt.TargetTolFun + Y.col(D).minCoeff());
			}
			else {
				d0hat = prctile(NormTabDis, opt.InidhatPrctile);
				globY_ind = find(Y.col(D).array() <= opt.TargetTolFun + join(archiveXf.col(D), Y.col(D)).minCoeff());	// Consider only global minima, among the recently generated solutions
			}
		}
		else
			globY_ind = find(Y.col(D).array() <= opt.TargetTolFun + Y.col(D).minCoeff());	// Consider only global minima, among the recently generated solutions
		for (int k : globY_ind) {	// Check the converged solutions one after another. This must be horizontal vector
			VectorXd dis = calcDis(archiveXf.leftCols(D), Y.row(k).head(D));	// Euclidean distance of the k-th Y to all the archived solutions
			VectorXi check_archive_ind = sort(dis, true);	// Check whether the k - th solution share the same basin with one of the archived solutions 
			bool isnew = true;	// By default, this solution is a new niche unless concluded otherwise
			for (int m : check_archive_ind) {	// Apply DetectMultimodal (hill-valley) function for the candidate solution and the archived solutionss. The closer ones are checked earlier
				if (isnew == true) {	// No sharing archived point is found yet
					bool isbas;
					int usedval;
					if (remained_eval - used_eval >= opt.DetectMultBudget)	// If we are not out of the evaluation budget
						isNewBasin(Y.row(k).head(D), archiveXf.row(m).head(D), Y(k, D), archiveXf(m, D), env,
							opt, isbas, usedval);
					else {
						isbas = true;
						usedval = 0;
					}
					used_eval += usedval;	// Add the used evalaution budget to the evalaution count
					if (isbas == false) {	// If it shares the same basin
						isnew = false;	// It is not a new basin
						N_rep(m)++;	// Update the number of  solutions that have converged to the basin of the m-th archived solution
						if (archiveXf(m, D) > Y(k, D)) {	// If the newly converged solution is fitter than the current archived soltuion, replace the current with the new one
							archiveXf.row(m) = Y.row(k);
							break;	// No need to check wrt other archived solutions
						}
					}	// The corresponding archived point was identified
				}
			}	// It is determined whether the k-th converged solution is a new basin or not
			if (isnew == true) {	// If the converged solution is a new basin
				archiveXf = joinRow(archiveXf, Y.row(k));	// Add it to the archive
				all_N_rep = join(all_N_rep, 0);	// No repitition yet (All restarts so far)
				N_rep = join(N_rep, 1e-14);	// No repition yet (Current restart)
				NormTabDis = join(NormTabDis, d0hat);	// Assign the default value of the normalized taboo distance 
			}
		}	// All of the (converged) subpopulations were analyzed. Now update the normalized taboo distance of the archived solutions
		size_t archive_size = NormTabDis.size();	// The archive size
		all_N_rep += N_rep;	// The overal number of the converged subpopulations to each archived point in all restarts
		double mean_rep = (1 - opt.TargetNewNiche) * globY_ind.size() / archive_size;	// Allowable fraction of subpopulations that may converge to the same basin without increasing NormTabDis of that basin
		for (size_t l = 1; l < archive_size; ++l) {	// Update the normalized taboo distance
			if (N_rep(l) >= mean_rep)	// Increase the normalized taboo distance
				NormTabDis(l) *= exp(log(1 + N_rep(l) - mean_rep) * opt.tau_d);
			else	// Decrease the normalized taboo distance
				NormTabDis(l) *= exp(-log(1 + mean_rep - N_rep(l)) * opt.tau_d);
		}
	}
}