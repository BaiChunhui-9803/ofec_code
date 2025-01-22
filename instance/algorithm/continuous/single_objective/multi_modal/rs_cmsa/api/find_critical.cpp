#include "find_critical.h"
#include "calc_dis.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	VectorXi findCritical(const MatrixXd &all_better, const VectorXd &all_better_NormTabDis, 
		const Ref<const RowVectorXd, 0, InnerStride<>> &Xmean,
		const VectorXd &Str, const MatrixXd &Cinv, double Smean, double criticality_treshhold, TypeDisMetric distance_metric)
	{
		VectorXi sorted_critical_ind;
		if (all_better_NormTabDis.size() > 0) {
			auto umean = (Str.array().pow(1.0 / Str.size())).prod();
			auto umax = Str.maxCoeff();
			auto L = calcDis(all_better, Xmean);
			VectorXd intU, intL;
			switch (distance_metric) {
			case TypeDisMetric::kMahalanobis:
				intU = (L + umean * Smean * all_better_NormTabDis) / (umax * Smean);
				intL = (L - umean * Smean * all_better_NormTabDis) / (umax * Smean);
				break;
			case TypeDisMetric::kEuclidean:
				intU = (L + umax * Smean * all_better_NormTabDis) / (umax * Smean);
				intL = (L - umax * Smean * all_better_NormTabDis) / (umax * Smean);
				break;
			default:
				throw Exception("Option for distance metric is not valid, aborting the run ...");
				break;
			}
			VectorXd criticality = normcdf(intU) - normcdf(intL);
			VectorXi ind = sort(criticality, false);
			int critical_number = (criticality.array() > criticality_treshhold).cast<int>().sum();
			if (critical_number > 0)
				sorted_critical_ind = ind.head(critical_number);
		}
		return sorted_critical_ind;
	}
}