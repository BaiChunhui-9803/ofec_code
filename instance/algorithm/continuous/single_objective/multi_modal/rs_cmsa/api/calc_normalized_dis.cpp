#include "calc_normalized_dis.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	VectorXd calcNormalizedDis(const Ref<const MatrixXd> &X, const Ref<const RowVectorXd, 0, InnerStride<>> &x,
		const RowVectorXd &str, const MatrixXd &Cinv, double S, const Opt &opt)
	{
		VectorXd dis;
		if (X.rows() == 0) {
			dis.resize(1);
			dis(0) = 1e100;
		}
		else {
			dis.resize(X.rows());
			if (opt.DisMetric == TypeDisMetric::kEuclidean) {
				auto mean_str = str.array().pow(1.0 / str.size()).prod();
				for (size_t k = 0; k < X.rows(); ++k) {
					auto u = X.row(k)- x;
					dis(k) = u.norm() / (S * mean_str);
				}
			}
			else if (opt.DisMetric == TypeDisMetric::kMahalanobis) {
				for (size_t k = 0; k < X.rows(); ++k) {
					auto u = X.row(k) - x;
					dis(k) = (u * Cinv * u.transpose()).cwiseSqrt()(0) / S;
				}
			}
			else
				throw Exception("Option for the distance metric is invalid, aborting ...");
		}
		return dis;
	}
}