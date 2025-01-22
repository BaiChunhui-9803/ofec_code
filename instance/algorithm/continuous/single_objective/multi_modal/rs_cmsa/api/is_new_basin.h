#ifndef OFEC_RS_CMA_ES_IS_NEW_BASIN_H
#define OFEC_RS_CMA_ES_IS_NEW_BASIN_H

#include "../../../../../../../utility/linear_algebra/eigen.h"
#include "../../../../../../../core/environment/environment.h"
#include "struct.h"

namespace ofec {
	class Problem;
	class Algorithm;

	namespace rs_cmsa {
		using namespace Eigen;

		void isNewBasin(RowVectorXd X1, RowVectorXd X4, double f1, double f4, Environment *env, const Opt &opt,
			bool &new_basin, int &iter);
	}
}

#endif // !OFEC_RS_CMA_ES_IS_NEW_BASIN_H