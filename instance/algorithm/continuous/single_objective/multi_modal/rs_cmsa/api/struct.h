#ifndef OFEC_RS_CMSA_ES_STRUCT_H
#define OFEC_RS_CMSA_ES_STRUCT_H

namespace ofec::rs_cmsa {
	enum class TypeDisMetric { kMahalanobis, kEuclidean };
	enum class TypeIniC { kEye, kDefault };
	enum class TypeWhichSubPop { kAll, kConverged };

	struct Opt {
		size_t initialPop;					// Coefficient for the initial subpopulation
		double TargetNewNiche;				// The desired fraction of identified global minima that should be new basins(The same as alpha_new in the  paper)
		size_t DetectMultBudget;			// Max evaluation budget for the hill - valley(DetectMultimodal) function
		double lambdaIncFactor;				// Population size increase factor
		TypeDisMetric DisMetric;			// The distance metric.Use 'Mahalanobis' or 'Euclidean'
		double EltRatio;					// The fraction of elite solutions
		double InitializeDisMult;			// Distance multiplier for the intialization
		double MaxIniStepSize;				// The start point for the initial global stepsize
		double TCovCoeff;					// Controls the learning rate for the covariance matrix
		TypeIniC iniC;						// 'eye' or 'default'
		double TSigmaCoeff;					// Controls the learning rate for the global Step size
		double TolHistFun;					// For termination and restarts. Use equal or smaller than the tightest desired function tolerance
		double TargetTolFun;				// Desired function tolerance
		double InidhatPrctile;				// Percentile for the normalized taboo distance of the current subpopulations
		double tolX;						// For termination criteria - not used by default
		double CriticThresh;				// The criticality threshold for taboo points
		double RedRatio;					// Ratio for temporary shrinkage of the taboo regions
		double d0hat;						// initla bvalue for the default normalized taboo distance
		double tau_d;						// Adaptation rate of the normalized taboo distances
		double MuLambdaRatio;				// Fraction of the parents
		double MaxCondC;					// Maximum condition number of the covariance matrix(A termination criterion)
		double ArchiveConsiderTabooTol;		// Tolerance for  being a better solution
		TypeWhichSubPop WhichSubPopAnalyze;	// Which subpopulations to be analyzed ? 'all' or 'Converged'
		size_t MaxIter;						// The maximum iterations per subpopulation in a restart
		bool SaveForResume;					// Allow for resuming the run if interrupted: Use 'yes' or 'no'
	};
}

#endif // !OFEC_RS_CMAS_ES_STRUCT_H
