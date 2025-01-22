#include "evoler.h"
#include <numeric>
#include <Eigen/Dense>
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/environment/environment.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS
#include "../../../../../../datum/datum_inclusion.h"



namespace ofec {
	void EVOLER::addInputParameters() {
		m_input_parameters.add("discrete length", new RangedSizeT(m_discrete_length, 10, 5000, 20));
		m_input_parameters.add("sampling length", new RangedSizeT(m_sampling_length, 2, 200, 3));
	}

	void EVOLER::run_(Environment *env) {
		size_t M = m_discrete_length;
		size_t N = m_discrete_length;
		for (size_t s = 1; s <= m_sampling_length; ++s) {
			reconstructRepresentation(M, N, s, env);
		}
	}

	void EVOLER::reconstructRepresentation(size_t M, size_t N, size_t s, Environment *env) {
		using namespace Eigen;

		std::vector<size_t> index_c;
		std::vector<size_t> index_r;
		{
			std::vector<size_t> index(N);
			std::iota(index.begin(), index.end(), 0);
			m_random->uniform.shuffle(index.begin(), index.end(), s);
			index_c.assign(index.begin(), index.begin() + s);
		}
		{
			std::vector<size_t> index(M);
			std::iota(index.begin(), index.end(), 0);
			m_random->uniform.shuffle(index.begin(), index.end(), s);
			index_r.assign(index.begin(), index.begin() + s);
		}
		MatrixXd C(M, s);
		MatrixXd R(s, N);




		{
			std::unique_ptr<SolutionBase> tmp_sol(env->problem()->createSolution());
			auto &x = dynamic_cast<VariableVector<>&>(tmp_sol->variableBase());
			auto &range0 = CAST_CONOP(env->problem())->range(0);
			auto &range1 = CAST_CONOP(env->problem())->range(1);
			auto sep0 = (range0.second - range0.first) / (M - 1);
			auto sep1 = (range1.second - range1.first) / (N - 1);
			for (size_t i = 0; i < M; ++i) {
				x[0] = range0.first + i * sep0;
				for (size_t j = 0; j < s; ++j) {
					x[1] = range1.first + index_c[j] * sep1;
					tmp_sol->evaluate(env);
					C(i, j) = tmp_sol->objective(0);
				}
			}
			std::map<size_t, size_t> r_index_c;
			for (size_t j = 0; j < s; ++j) {
				r_index_c[index_c[j]] = j;
			}
			for (size_t i = 0; i < s; ++i) {
				x[0] = range0.first + index_r[i] * sep0;
				for (size_t j = 0; j < N; ++j) {
					if (r_index_c.count(j)) {
						R(i, j) = C(index_r[i], r_index_c[j]);
					}
					else {
						x[1] = range1.first + j * sep1;
						tmp_sol->evaluate(env);
						R(i, j) = tmp_sol->objective(0);
					}
				}
			}
		}
		MatrixXd U(s, s);
		for (size_t i = 0; i < s; ++i) {
			for (size_t j = 0; j < s; ++j) {
				U(i, j) = R(i, index_c[j]);
			}
		}
		BDCSVD<MatrixXd> svd(U, ComputeThinU | ComputeThinV);
		MatrixXd U_u = svd.matrixU();
		MatrixXd S_u = svd.singularValues().asDiagonal();
		MatrixXd V_u = svd.matrixV();
		MatrixXd Uu = U_u.transpose();
		size_t s_0 = s;
		CompleteOrthogonalDecomposition<MatrixXd> cod(S_u.block(0, 0, s_0, s_0));
		MatrixXd U1 = V_u.leftCols(s_0) * cod.pseudoInverse() * Uu.topRows(s_0);
		MatrixXd Z_est = C * U1 * R;
		m_grid_obj.resize(M);
		for (size_t i = 0; i < M; ++i) {
			m_grid_obj[i].resize(N);
			for (size_t j = 0; j < N; ++j) {
				m_grid_obj[i][j] = Z_est(i, j);
			}
		}
		//datumUpdated(env);
		//std::cout << "svd.singularValues():" << std::endl;
		//std::cout << svd.singularValues() << std::endl;
	}
}
