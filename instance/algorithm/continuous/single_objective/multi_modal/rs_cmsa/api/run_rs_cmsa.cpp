#include "run_rs_cmsa.h"
#include "initialize_subpop.h"
#include "peak2peak.h"
#include "evolve_subpop.h"
#include "analyze_converged_sol.h"
#include "../../../../../../../core/problem/continuous/continuous.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../../datum/datum_inclusion.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	void runRSCMSA(size_t lambda, size_t NsubP, size_t D, const RowVectorXd &D_X, const RowVectorXd &U_X,
		size_t &remained_eval, MatrixXd &all_converged, MatrixXd &archiveXf, VectorXd &all_N_rep,
		VectorXd &NormTabDis, const Opt &opt, Environment *env, Random *rnd, double &average_used_iter,
		MultiPopulation<Population<Solution<>>> &subpops)
	{
		MatrixXd archiveX;	VectorXd archivef;
		if (NormTabDis.size() == 0)
			archiveX.resize(0, D);
		else {
			archiveX = archiveXf.leftCols(D); archivef = archiveXf.col(D);
		}
		int eltNo = ceil(lambda * opt.EltRatio);			// No. of elite Solutions per subpopulation
		int mu = floor(lambda * opt.MuLambdaRatio);			// No. of parents
		RowVectorXd weights = log(1 + mu) - RowVectorXd::LinSpaced(mu, 1, mu).array().log();
		weights = weights / weights.array().sum();
		int IterNo = 1;
		MatrixXd Xmean, Cini;
		VectorXd Smean, default_NormTabDis;
		initializeSubpop(archiveX, NormTabDis, NsubP, D_X, U_X, opt, rnd, Xmean, Smean, Cini, default_NormTabDis);
#ifdef OFEC_DATUM_RADIUS_REPULSION_H
		if (NormTabDis.size() > 0) {
			g_radius_repulsion.value = NormTabDis(0) * Smean(0);
			env->algorithm()->datumUpdated(env, g_radius_repulsion);
		}
#endif // OFEC_DATUM_RADIUS_REPULSION_H
		int noimpsize = 120 + floor(30.0 * D / lambda);		// No improvement(for termination criterion)
		int stallsize = 10 + floor(30.0 * D / lambda);		// Stall time(for termination criterion)
		std::vector<MatrixXd> eltX(NsubP), eltZ(NsubP), C(NsubP);
		std::vector<VectorXd> eltf(NsubP), eltS(NsubP);
		std::vector<VectorXd> bestf_NE(NsubP), medf_NE(NsubP);
		std::vector<VectorXi> superior_subpop(NsubP);
		MatrixXd bestXf(NsubP, D + 1);
		for (size_t k = 0; k < NsubP; ++k) {
			eltX[k] = Xmean.row(k).replicate(eltNo, 1);	// Repetition of the subpopulation mean for the initial elite solutions of the subpopulation  
			eltZ[k] = MatrixXd::Zero(eltNo, D);				// Initial value of the elite variation vector
			eltf[k] = VectorXd::Ones(eltNo) * 1e170;
			eltS[k] = VectorXd::Ones(eltNo) * Smean.mean();	// Initial value of the elite global step sizes
			C[k] = Cini;
			bestf_NE[k] = VectorXd::LinSpaced(noimpsize, noimpsize, 1) * 1e140;	// History of the best of non-elite solutions for each subpopulation
			medf_NE[k] = VectorXd::LinSpaced(noimpsize, noimpsize, 1) * 1e150;	// History of the median of non-elite solutions for each subpopulation
			superior_subpop[k] = VectorXi::LinSpaced(k, 0, k - 1);
			bestXf.row(k) << Xmean.row(k), 1e170;
		}
		size_t activeP_No = NsubP; // The number of subpopulations not terminated so far
		VectorXi converged_subpop; // Subpopulations that have converged
		VectorXi stalled_subpop;
		VectorXi active_subpop = VectorXi::LinSpaced(NsubP, 0, NsubP - 1);
		VectorXd used_iter = VectorXd::Zero(NsubP);
		// Main loop starts here
		while (IterNo <= opt.MaxIter && remained_eval >= lambda / 2 && activeP_No >= 1) {
			activeP_No = 0;
			for (int m : active_subpop) {
				SelfAdjointEigenSolver<MatrixXd> eigensolver(C[m]);
				MatrixXd Rot = eigensolver.eigenvectors();
				VectorXd Str2 = eigensolver.eigenvalues();
				Str2 = Str2.cwiseMax(0);
				VectorXd Str = Str2.array().sqrt();
				double largest_eig = Smean(m) * Str.maxCoeff();
				double condC = pow(Str.maxCoeff() / Str.minCoeff(), 2);
				double max_diff = peak2peak(bestf_NE[m].tail(stallsize));
				double min_imp_best = median(bestf_NE[m].segment(bestf_NE[m].size() - noimpsize + 1, 20)) - median(bestf_NE[m].tail(20));
				double min_imp_med = median(medf_NE[m](seq(indexing::last - noimpsize + 1, indexing::last - noimpsize + 20)))
					- median(medf_NE[m](seq(indexing::last - 19, indexing::last)));
				ArrayXb active_chk(6);	// Checking the termination criteria for the subpopulation
				active_chk << (largest_eig > opt.tolX), (condC < opt.MaxCondC), (remained_eval >= lambda),
					(max_diff > opt.TolHistFun), (min_imp_best > 0), (min_imp_med > 0);
				size_t used_evals;
				if (active_chk.minCoeff() == true) {	// Do not terminate this subpopulation
					used_iter(m) = IterNo;
					activeP_No++;
					evolveSubpop(Xmean(superior_subpop[m], indexing::all), default_NormTabDis(superior_subpop[m]),
						Xmean.row(m), Smean(m), lambda, mu, eltX[m], eltf[m], eltS[m], eltZ[m],
						D_X, U_X, Rot, Str, C[m], bestf_NE[m], medf_NE[m], env, rnd, 
						archiveX, archivef, NormTabDis, weights, opt, bestXf.row(m), used_evals, subpops[m]);
				}
				else {	// Terminate the subpopulation
					used_evals = 0;
					if ((active_chk(3) == true || active_chk(2) == true) && bestXf(m, indexing::last) < 1e100) {
						converged_subpop = join(converged_subpop, m);	// The subpopulation has converged or the evaluation budget is finished
						Xmean.row(m) = bestXf.row(m).head(D);
					}
					else {
						stalled_subpop = join(stalled_subpop, m);		// Diverged or stalled subpopulation
						Xmean.row(m) = bestXf.row(m).head(D);
					}
				}
				remained_eval -= used_evals;
			}	// One iteration of the restart was completed
			noimpsize = 120 + floor(30.0 * D / lambda + 0.2 * IterNo);	// Update this parameter, since it depends on the the iteration number
			active_subpop = setdiff(active_subpop, join(converged_subpop, stalled_subpop));	// Subpopulations that have not been terminated
			VectorXi ind0 = sort(bestXf(active_subpop, D));
			active_subpop = active_subpop(ind0);	// Reorder active subpopulations according to their fitness
			for (int m : active_subpop) {
				superior_subpop[m] = find(bestXf.col(D).array() < bestXf(m, D));	// Find the superior subpopulations wrt each subpopulation
			}
			IterNo++;
#ifdef OFEC_DATUM_MULTI_POP_H
			for (int i = 0; i < subpops.size(); ++i) subpops[i].setActive(false);
			for (int i : active_subpop) subpops[i].setActive(true);
			g_multi_pop.pops.clear();
			g_multi_pop.pops.resize(active_subpop.size());
			for (size_t k = 0; k < active_subpop.size(); ++k) {
				for (size_t i = 0; i < subpops[active_subpop[k]].size(); ++i) {
					if (CAST_CONOP(env->problem())->boundaryViolated(subpops[active_subpop[k]][i]))
						continue;
					g_multi_pop.pops[k].push_back(&subpops[active_subpop[k]][i]);
				}
			}
			env->algorithm()->datumUpdated(env, g_multi_pop);
#endif
		}
		all_converged = joinRow(all_converged, bestXf);
		VectorXi analyze_subpop;
		if (opt.WhichSubPopAnalyze == TypeWhichSubPop::kConverged)
			analyze_subpop = converged_subpop;
		else if (opt.WhichSubPopAnalyze == TypeWhichSubPop::kAll)
			analyze_subpop = VectorXi::LinSpaced(NsubP, 0, NsubP - 1);
		else
			throw Exception("The selected option for WhichSubPopAnalyze is not supported."
				"Valid options are kAll and kConverged");
		if (analyze_subpop.size() > 0) {
			size_t used_eval_isbasin;
			analyzeConvergedSol(bestXf(analyze_subpop, indexing::all), env, archiveXf,
				all_N_rep, NormTabDis, opt, remained_eval, used_eval_isbasin);
			remained_eval -= used_eval_isbasin;
		}
		average_used_iter = used_iter.mean();
	}
}