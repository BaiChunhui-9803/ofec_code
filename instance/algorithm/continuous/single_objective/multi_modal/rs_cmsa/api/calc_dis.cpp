#include "calc_dis.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	VectorXd calcDis(const Ref<const MatrixXd> &X, const Ref<const RowVectorXd, 0, InnerStride<>> &x) {
		VectorXd dis;
		auto n = X.rows();
		if (n != 0) {
			dis.resize(n);
			for (size_t k = 0; k < n; ++k)
				dis[k] = (X.row(k) - x).norm();
		}
		return dis;
	}
}