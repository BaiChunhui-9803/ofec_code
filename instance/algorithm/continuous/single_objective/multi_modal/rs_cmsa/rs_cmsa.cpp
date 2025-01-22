#include "rs_cmsa.h"
#include "../../../../../../core/problem/continuous/continuous.h"

#include "api/run_rs_cmsa.h"

namespace ofec {
	void RS_CMSA::initialize_(Environment *env) {
		Algorithm::initialize_(env);
		if (m_maximum_evaluations <= 0) {
			throw Exception("Maximum number of evaluations must be provided.");
		}
		m_remained_eval = m_maximum_evaluations;
		m_targetNsubP = 10;
		m_opt.initialPop = 3;
		m_opt.TargetNewNiche = 0.5;
		m_opt.DetectMultBudget = 10;
		m_opt.lambdaIncFactor = 1.01;
		m_opt.DisMetric = rs_cmsa::TypeDisMetric::kMahalanobis;
		m_opt.EltRatio = 0.15;
		m_opt.InitializeDisMult = 2;
		m_opt.MaxIniStepSize = 0.25;
		m_opt.TCovCoeff = 1;
		m_opt.iniC = rs_cmsa::TypeIniC::kDefault;
		m_opt.TSigmaCoeff = 1;
		m_opt.TolHistFun = 1e-5;
		m_opt.TargetTolFun = 1e-5;
		m_opt.InidhatPrctile = 25;
		m_opt.tolX = 0;
		m_opt.CriticThresh = .01;
		size_t D = CAST_CONOP(env->problem())->numberVariables();
		m_opt.RedRatio = pow(1.01, 1.0 / D);
		m_opt.d0hat = 1;
		m_opt.tau_d = pow(1.0 / D, 0.5);
		m_opt.MuLambdaRatio = 0.5;
		m_opt.MaxCondC = 1e14;
		m_opt.ArchiveConsiderTabooTol = 0;
		m_opt.WhichSubPopAnalyze = rs_cmsa::TypeWhichSubPop::kAll;
		m_opt.MaxIter = std::numeric_limits<size_t>::max();
		m_opt.SaveForResume = true;
	}
	
	void RS_CMSA::run_(Environment *env) {
		size_t D = CAST_CONOP(env->problem())->numberVariables();
		Eigen::RowVectorXd D_X(D), U_X(D);
		for (size_t j = 0; j < D; ++j) {
			D_X[j] = CAST_CONOP(env->problem())->range(j).first;
			U_X[j] = CAST_CONOP(env->problem())->range(j).second;
		}
		size_t pop0 = m_opt.initialPop * sqrt(D);
		Eigen::VectorXd NormTabDis;		// The normalized taboo distances of the archived solutions (vector)
		Eigen::VectorXd all_N_rep;		// The number of times a basin was identified by a subpopulation minus one (vector). Not used by the algorithm, only for monitoring the optimization process 
		Eigen::MatrixXd all_converged;	// The best solution of all subpopulations at the end of all restarts (matrix)
		Real lambda_exponent = 0;		// Used for increasing popultion size (scalar)
		size_t NsubP = m_targetNsubP;	// The number of subpopulations in the next restart(scalar)
		size_t restartNo = 0;			// The restart number(scalar)
		Eigen::MatrixXd archiveXf;		// archived solutions and their function values (matrix)	
		size_t lambda = pop0;			// Subpopulation size (Scalar)
		all_converged.resize(0, D + 1);
		archiveXf.resize(0, D + 1);
		
		while (!terminating()) {
			if (m_remained_eval < m_opt.DetectMultBudget * (archiveXf.rows() + 1) + 2 * lambda) {
				break;
			}
			lambda = floor(pop0 * pow(m_opt.lambdaIncFactor, lambda_exponent));
			restartNo++;
			Real average_used_iter;
			m_subpops.resize(NsubP, lambda, env, D);
			rs_cmsa::runRSCMSA(lambda, NsubP, D, D_X, U_X, m_remained_eval, all_converged, archiveXf, all_N_rep,
				NormTabDis, m_opt, env, m_random.get(), average_used_iter, m_subpops);
			auto max_nsub_p = floor(m_remained_eval / (m_opt.lambdaIncFactor * average_used_iter * lambda));
			if (max_nsub_p >= m_targetNsubP) {									// Sufficient budget
				NsubP = m_targetNsubP;
				lambda_exponent++;
			}
			else if (max_nsub_p >= (m_targetNsubP / m_opt.lambdaIncFactor)) {	// Increase lambda at a slower rate
				NsubP = m_targetNsubP;
				auto adjusted_lambda_mult = m_remained_eval / (NsubP * average_used_iter * lambda);
				lambda_exponent += log2(adjusted_lambda_mult);
			}
			else if (max_nsub_p < (m_targetNsubP / m_opt.lambdaIncFactor))		// Do not increase lambda, reduce NsubP
				NsubP = std::max<size_t>(1, floor(m_remained_eval / (average_used_iter * lambda)));
		}
	}
}