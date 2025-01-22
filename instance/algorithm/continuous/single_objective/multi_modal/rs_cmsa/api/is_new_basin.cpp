#include "is_new_basin.h"
#include "../../../../../../../core/problem/solution.h"
#include "../../../../../../../core/problem/problem.h"

namespace ofec::rs_cmsa {
	using namespace Eigen;

	void isNewBasin(RowVectorXd X1, RowVectorXd X4, double f1, double f4, Environment *env, const Opt &opt,
		bool &new_basin, int &iter)
	{
		new_basin = false; iter = 0;
		Solution<> tmp_sol(1, 0, X1.size());
		OptimizeMode opt_mode = env->problem()->optimizeMode(0);
		size_t D = X1.size();
		if ((X1 - X4).norm() > opt.tolX) {
			auto benchf = std::max(f1, f4);
			iter = 1;
			RowVectorXd dir = X4 - X1;	// search direction
			auto delta = dir.norm();
			dir = dir / delta / 2.618;
			RowVectorXd X2 = X1 + delta * dir;
			for (size_t j = 0; j < D; ++j)	tmp_sol.variable()[j] = X2(j);
			tmp_sol.evaluate(env);
			auto f2 = opt_mode == OptimizeMode::kMinimize ? tmp_sol.objective(0) : -tmp_sol.objective(0);
			if (f2 > benchf)
				new_basin = true;
			else if (iter < opt.DetectMultBudget) {
				RowVectorXd X3 = X1 + 1.618 * delta * dir;
				for (size_t j = 0; j < D; ++j)	tmp_sol.variable()[j] = X3(j);
				tmp_sol.evaluate(env);
				auto f3 = opt_mode == OptimizeMode::kMinimize ? tmp_sol.objective(0) : -tmp_sol.objective(0);
				iter++;	// revised 6 / 22 / 2016
				if (f3 > benchf)
					new_basin = true;
				else {
					while (iter < opt.DetectMultBudget && (X1 - X4).norm() > opt.tolX) {
						if (f2 < f3) {
							delta = (X4 - X2).norm() / 2.618;
							X1 = X2; f1 = f2;
							X2 = X3; f2 = f3;
							X3 = X1 + 1.618 * delta * dir;
							for (size_t j = 0; j < D; ++j)	tmp_sol.variable()[j] = X3(j);
							tmp_sol.evaluate(env);
							auto f3 = opt_mode == OptimizeMode::kMinimize ? tmp_sol.objective(0) : -tmp_sol.objective(0);
							iter++;
						}
						else if (f2 > f3) {
							delta = (X3 - X1).norm() / 2.618;
							X4 = X3; f4 = f3;
							X3 = X2; f3 = f2;
							X2 = X1 + delta * dir;
							for (size_t j = 0; j < D; ++j)	tmp_sol.variable()[j] = X2(j);
							tmp_sol.evaluate(env);
							auto f2 = opt_mode == OptimizeMode::kMinimize ? tmp_sol.objective(0) : -tmp_sol.objective(0);
							iter++;
						}
						else {
							delta = (X3 - X2).norm() / (2.618034);
							X1 = X2; f1 = f2;
							X4 = X3; f4 = f3;
							X2 = X1 + delta * dir;
							X3 = X4 - delta * dir;
							for (size_t j = 0; j < D; ++j)	tmp_sol.variable()[j] = X2(j);
							tmp_sol.evaluate(env);
							auto f2 = opt_mode == OptimizeMode::kMinimize ? tmp_sol.objective(0) : -tmp_sol.objective(0);
							for (size_t j = 0; j < D; ++j)	tmp_sol.variable()[j] = X3(j);
							tmp_sol.evaluate(env);
							auto f3 = opt_mode == OptimizeMode::kMinimize ? tmp_sol.objective(0) : -tmp_sol.objective(0);
							iter += 2;
						}
						if (std::max(f2, f3) > benchf) {	// Added on 6 / 22 / 2016
							new_basin = true;
							break;
						}
					}
				}	// f3 checked
			}	// f2 checked
		}	//  % tolX checked
	}
}