
#ifndef SELECTION_PROBLEM_MAP_H
#define SELECTION_PROBLEM_MAP_H

#include "../../../../core/global.h"
#include "../../../../core/problem/problem.h"
#include "../../../../core/problem/optima.h"
#include "../../../../core/problem/domain.h"
#include "../../../../core/problem/solution.h"
#include "../../../../utility/memory_record/memory_record.h"
#include "../../../../utility/linear_algebra/vector.h"
#include "../../../../utility/functional.h"
#include "../../../../utility/function/custom_function.h"
#include "../../../../core/problem/uncertainty/dynamic.h"
#include "../../../../core/problem/uncertainty/noisy.h"
#include "ttFun.h"
#include "ttFunPar.h"
#include <string>
#include <map>
#include <vector>
#include <array>
namespace ofec {
	namespace sp {
		enum class PosInfoIdx {
			WALL_IDX = 0,
			POS_FROM_IDX = 1,
			POS_TO_IDX = 2,
			COR_IDX = 3
		};
		struct Point {
			std::array<int, 4> m_pos_info =
			{ -1,-1,-1,-1 };
			DistributionFunTt m_angle_distr;
			DistributionFunTt m_radius_distr;
			DistributionFunTt m_price_distr;
			FeasibleDistributionFunTt m_feasible_distr;
		};
		// the feasible function can change the framework of the graph
		struct Edge {

			std::array<int, 3> m_pos_info =
			{ -1,-1,-1 };
			DistributionFunTt m_distance_distr;
			//	DistributionFunTt m_price_distr;
			FeasibleDistributionFunTt m_feasible_distr;
		};



		class Wall {
		public:
			FeasibleDistributionFunTt m_feasible_distr;
			Distribution3DFunTt m_mesh_price;
		protected:
			TypeWall m_type = TypeWall::SMALL_LOOP;
			std::array<int, 1> m_pos_info;
			int m_numPos = 0;
			double m_center;
			std::array<double, 3> m_radius_range;
			std::array<double, 3> m_angle_range;
			std::vector<int> m_center_ids;

			Vector getPosition(double radius, double angle) const {
				Vector pos(3, 0);
				pos[0] = m_center;
				pos[1] = radius * sin(angle);
				pos[2] = radius * cos(angle);
				return std::move(pos);
			}


			bool getCenterPos(double T, double t,
				const Point& p, double& radius, double& angle)const {
				if (m_type == TypeWall::NEEDLE &&
					p.m_pos_info[static_cast<int>(PosInfoIdx::POS_FROM_IDX)]
					== getCenterPosId(T, t)) {
					radius = m_radius_range[0];
					angle = (m_angle_range[0] + m_angle_range[1]) / 2.0 + OFEC_PI;
					return true;
				}
				else return false;
			}

			void setCurPosVal(double& curVal, double from, double to)const;
			inline void moidifyCurPos(double& radius, double& angle)const {
				setCurPosVal(radius, m_radius_range[0], m_radius_range[1]);
				if (m_type == TypeWall::NEEDLE) {
					setCurPosVal(angle, m_angle_range[0], m_angle_range[1]);
				}
			}
		public:

			std::array<int, 1>& getPosInfo() {
				return m_pos_info;
			}
			const std::array<double, 3>& radiusRange()const {
				return m_radius_range;
			}

			const std::array<double, 3>& angleRange()const {
				return m_angle_range;
			}

			inline double getCenter()const {
				return m_center;
			}
			Vector outerPosition(Real angle) const {
				return std::move(getPosition(m_radius_range[1], angle));
			}
			void PositionToXY(const Vector& p, double& x, double& y)const {
				x = p[1];
				y = p[2];
			}
			// limited: T \in  [0, 1] , t\in [0,1]
			int getCenterPosId(double T, double t) const {
				return m_center_ids[UTILITY::getCurSize(T, m_center_ids.size())];
			}
			inline int numPos()const {
				return m_numPos;
			}

			Vector getPosition(const Point& p, double T, double t, Random *rnd = -1) const {
				Real radius(0), angle(0);
				if (!getCenterPos(T, t, p, radius, angle)) {
					radius = p.m_radius_distr.getVal(T, t, rnd);
					angle = p.m_angle_distr.getVal(T, t, rnd);
					moidifyCurPos(radius, angle);
				}
				//radius = p.m_radius_distr.getVal(T, t, rnd);
				//angle = p.m_angle_distr.getVal(T, t, rnd);
				//moidifyCurPos(radius, angle);
				return std::move(getPosition(radius, angle));
			}
			//Vector getPosition(double T, double t, Random *rnd, const Point& p) const {				
			//	Real radius(0), angle(0);
			//	if (!getCenterPos(T, t, p, radius, angle)){ 
			//		radius=p.m_radius_distr.getVal(T, t,rand);
			//		angle= p.m_angle_distr.getVal(T, t, rnd);
			//		moidifyCurPos(radius, angle);
			//	}
			//	return std::move(getPosition(radius, angle));
			//}

			void setPositionStyle(Random *rnd,
				double center,
				TypeWall type,
				double max_radius = 2.0,
				double small_loop_radius = 0.1,
				double large_loop_inner = 0.1);
			int calNumPos(Random *rnd,
				TypeWall type,
				int num_min_points = 5,
				int num_small_points = 10,
				int num_large_points = 20,
				int num_circle_points = 50,
				int num_needle_points = 30);

			void setNumPos(Random *rnd,
				int num_min_points = 5,
				int num_small_points = 10,
				int num_large_points = 20,
				int num_circle_points = 50,
				int num_needle_points = 30);

			void setNumPos(const Wall& wall) {
				m_numPos = wall.m_numPos;
			}

			void initDynamicIdxs(Random *rnd, ParaInfo& para) {
				m_center_ids;
				para.m_dynamic_info.initDynamicIdxs(rnd, m_numPos, m_center_ids);
			}
		};

		struct TestInfo {

			std::vector<std::pair<std::string, double>> path_info_info;
			std::vector<std::pair<std::string, ofec::Vector>> path_info_pos;

			void clear() {
				path_info_info.clear();
				path_info_pos.clear();
			}
			void push_back(const std::string& name,
				double val
			) {
				path_info_info.emplace_back(name, val);
			}
			void push_back(const std::string& name,
				Vector val
			) {
				path_info_pos.emplace_back(name, val);
			}

		};
	//	extern TestInfo m_test_info;

		struct SolutionSPInfo {
			std::vector<Real> m_type_dis;
			std::vector<Real> m_type_val;
			Real m_sum_time;
			Vector m_cur_pos;
			std::pair<int, int> m_cur_dimIdx;
			Real m_objective;


			std::vector<TestInfo> m_test_infos;
		//	std::vector<std::pair<int, double>> m_pos_dis;
		//	std::vector<std::vector<std::pair<std::string, double>>> m_path_info;
		};


		class MapInfo {
		public:
			using solutionType = Solution<VariableVector<int>>;
			using variableType = VariableVector<int>;
		protected:
			std::vector<Wall> m_wall;
			std::vector<std::vector<Point>> m_stations;
			std::vector<std::vector<std::vector<Edge>>> m_edges;

			// the traffic jam scale = [0,1] 
			std::vector<std::vector<std::vector<double>>> m_edges_scale;
			std::vector<std::pair<double, double>> m_discount;
			//	Real m_feasible_pro = 0.9;
			double m_max_distance = 0.0;

			std::vector<double> m_dim_distance;
			//	Real m_radius_moving_speed = 3.0;
			//	Real m_angle_moving_speed = 2 * OFEC_PI;


			// for env noisy
			std::vector<solutionType> m_basic_sols;
			std::vector<std::vector<bool>> m_sols_factors_dims;
			std::vector<int> m_sols_factors_num_dim;

			double m_env_noisy_ratio = 0.0;
			double m_var_noisy_ratio = 0.0;

			// dijstra_info
//	std::vector<std::vector<SolutionSPInfo>> m_solution_spInfos;


			

			// 
			void calMaxDistance();
			double getPrice(int dim, int cur, double x, double y,
				double T, double t, Random *rnd = -1)const {
				return 	m_stations[dim][cur].m_price_distr.
					getVal(T, t, rnd)
					* m_wall[dim].m_mesh_price.getVal(x, y, T, t, rnd);
			}



			// add noisy 
			void transferXEnvNoisy(
				std::vector<int>& x,
				double T,
				Random *rnd = -1, bool flag_env_noisy = false)const;

			void transferXVarNoisy(
				std::vector<int>& x,
				double T,
				Random *rnd = -1, bool flag_var_noisy = false)const;


			void transferXWorkingSol(
				std::vector<int>& x,
				const std::vector<int>& working_x
			)const;

			void transferXNoisy(
				double noisy_ratio,
				bool flag_scene,
				int numDims,
				std::vector<int>& x,
				double T,
				Random *rnd = -1,
				bool flag_noisy = false)const;

		public:

			const solutionType& getBasicSols(double T)const
			{
				return m_basic_sols
					[UTILITY::getCurSize(T, m_basic_sols.size())];
			}
			const std::vector<bool>& getScenesDim(double T)const {
				return m_sols_factors_dims
					[UTILITY::getCurSize(T, m_sols_factors_dims.size())];
			}

			int getScenesNumDims(double T)const {
				return m_sols_factors_num_dim
					[UTILITY::getCurSize(T, m_sols_factors_num_dim.size())];
			}
			const std::vector<Wall>& getWall()const {
				return m_wall;
			}
			const std::vector<std::vector<Point>>& getStations() const {
				return m_stations;
			}
			const std::vector<std::vector<std::vector<Edge>>>&
				getEdges()const {
				return m_edges;
			}

			inline double getCenter(int dim)const {
				return m_wall[dim].getCenter();
			}
			int numWall()const {
				return m_wall.size();
			}
			int numCandidates(int dim)const {
				return m_wall[dim].numPos();
			}

			int getCorId(int dim, int from) const{
				return m_stations[dim][from].
					m_pos_info[static_cast<int>(PosInfoIdx::COR_IDX)];
			}
			inline Vector getPos(int dim, int cur, double T, double t, Random *rnd = -1) const {
				return m_wall[dim].getPosition(m_stations[dim][cur], T, t, rnd);
			}
			double getPrice(int dim, int cur,
				double T, double t, Random *rnd = -1)const {
				double x(0), y(0);
				m_wall[dim].PositionToXY(
					getPos(dim, cur, T, t, rnd), x, y);
				return getPrice(dim, cur, x, y,
					T, t, rnd);
			}
			double getPrice(int dim, int cur, const Vector& pos,
				double T, double t, Random *rnd = -1) const {
				double x(0), y(0);
				m_wall[dim].PositionToXY(pos, x, y);
				return getPrice(dim, cur, x, y,
					T, t, rnd);
			}

			double getDistance(int dim, int from, int to,
				const Vector& fromPos, const Vector& toPos,
				double T, double t, Random *rnd = -1) const{
				return fromPos.distance(toPos) *
					m_edges[dim][from][to].
					m_distance_distr.getVal(T, t, rnd);
			}
			double getDistance(int dim, int from, int to,
				double T, double t, Random *rnd = -1
			)const {
				Vector fromPos(getPos(dim, from, T, t, rnd));
				Vector toPos(getPos(dim + 1, to, T, t, rnd));
				return getDistance(dim, from, to, fromPos, toPos,
					T, t, rnd);
				///	return fromPos.distance(toPos) * m_edges[dim][from][to].getCurScale(T, t, rnd->uniform);
			}
			void updatePosInfo(Random *rnd);
			void initializeWallType(Random *rnd,
				const ParaInfo& para);
			void initSol(std::vector<int>& x, Random *rnd)const {
				x.resize(numWall());
				for (int i(0); i < x.size(); ++i) {
					x[i] = rnd->uniform.nextNonStd<int>(0, numCandidates(i));
				}
			}
			void initSol(std::vector<int>& x, Random *rnd)const {
				x.resize(numWall());
				for (int i(0); i < x.size(); ++i) {
					x[i] = rnd->uniform.nextNonStd<int>(0, numCandidates(i));
				}
			}

			
			void initEnvDynamics(Random *rnd, ParaInfo& para) {
				m_basic_sols.resize(std::max<size_t>(1,para.m_dynamic_info.m_num_basic_sols));
				initSol(m_basic_sols.front().variable().vect(), rnd);
				{

					auto pro(para.m_dynamic_info.m_change_env_factors_raito);
					int valnumWall(numWall());
					for (size_t idT(1); idT < m_basic_sols.size(); ++idT) {
						auto& x = m_basic_sols[idT].variable().vect();
						x = m_basic_sols[idT - 1].variable().vect();
						while (--valnumWall && rnd->uniform.next() < pro) {
							int idD = rnd->uniform.nextNonStd<int>(0, valnumWall);
							x[idD] =
								rnd->uniform.nextNonStd<int>(0, numCandidates(idD));
						}
					}
				}

				auto& factors_dims(m_sols_factors_dims);
				factors_dims.resize(std::max<size_t>(1, para.m_dynamic_info.m_num_feasible_basic_sols_states));
				for (auto& it : factors_dims) {
					it.resize(numWall());

					for (int idx(0); idx < it.size(); ++idx) {
						it[idx] = false;
					}
					//for (auto& it2 : it) it2 = false;
				}
				auto& factors_numdims(m_sols_factors_num_dim);
				factors_numdims.resize(factors_dims.size());
				{
					{
						auto& x = factors_dims.front();
						auto& pro = para.m_dynamic_info.m_init_env_factors_dim_raito;
						for (int idx(0); idx < numWall(); ++idx) {
							if (rnd->uniform.next() < pro) {
								x[idx] = true;
							}
						}
					}
	
					auto pro(para.m_dynamic_info.m_change_env_factors_dim_raito);
					int valNumWall(numWall());
					for (size_t idT(1); idT < factors_dims.size(); ++idT) {
						auto& x = factors_dims[idT];
						x = factors_dims[idT - 1];
						while (--valNumWall && rnd->uniform.next() < pro) {
							int idD = rnd->uniform.nextNonStd<int>(0, valNumWall);
							x[idD] = !x[idD];
						}
					}

					for (size_t idT(0); idT < factors_dims.size(); ++idT) {
						int numDim(0);
						for (int idx(0); idx < factors_dims[idT].size(); ++idx) {
							if (factors_dims[idT][idx]) ++numDim;
						}
						//for (auto& it : factors_dims[idT])if (it) ++numDim;
						factors_numdims[idT] = numDim;
					}
				}

				m_env_noisy_ratio = para.m_dynamic_info.m_env_noisy_ratio;
			}
			void initializeDynamicInfo(Random *rnd, ParaInfo& para);
			void initialize(Random *rnd, ParaInfo& para);
			double get_heuristic_info(int dim, int from, int to,
				double T, Random *rnd = -1) {
				if (dim == -1) {
					double cur_time = 0;
					return getDistance(0, 0, 0, T, cur_time, rnd) *
						//getCurPrice(Random *rnd,double T , double tint dim,int cur) 
						getPrice(dim + 1, to, T, cur_time, rnd);
				}
				else {
					double cur_time = 1.0 * (dim + 1) / static_cast<double>(m_wall.size());
					return getDistance(dim, from, to, T, cur_time, rnd) *
						getPrice(dim, from, T, cur_time, rnd);
				}
			}


			void initSpInfo(SolutionSPInfo& info)const {
				info.m_type_dis.resize(m_discount.size());
				info.m_type_val.resize(m_discount.size());
				for (auto& it : info.m_type_dis) it = 0;
				for (auto& it : info.m_type_val) it = 0;
				info.m_sum_time = 0;
				info.m_cur_dimIdx.first = -1;
				info.m_objective = 0;
			//	m_test_info.clear();
				info.m_test_infos.clear();
				
			}
			void transferXfactors(std::vector<int>& x, double T)const;
			void addSpInfoPos(SolutionSPInfo& info,
				const std::vector<int>& cx,
				int to_dim,
				double T,
				Random *rnd ,
				bool flag_correlation,
				bool flag_obj_noisy,
				bool flag_var_noisy,
				bool flag_env_noisy,
				bool flag_change_online,
				const std::vector<int>* const working_sol=nullptr
			) const {
				if (to_dim <= 0)return;
				if (!flag_obj_noisy) {
					rnd = -1;
				}


				if (info.m_cur_dimIdx.first == -1) {
					info.m_cur_dimIdx.first = 0;
					info.m_cur_dimIdx.second = cx[0];
					info.m_cur_pos = getPos(
						info.m_cur_dimIdx.first,
						info.m_cur_dimIdx.second,
						T, info.m_sum_time, rnd
					);
				//	m_test_info.clear();
				//	m_test_info.path_info_pos.emplace_back("1", info.m_cur_pos);
				//	info.m_test_infos.push_back(m_test_info);
					
				}

				std::vector<int> x = cx;
				//transferXfactors(x, T);
				transferXVarNoisy(x, T, rnd, flag_var_noisy);
				transferXEnvNoisy(x, T, rnd, flag_env_noisy);
				if (working_sol != nullptr) {
					transferXWorkingSol(x, *working_sol);
				}


				Vector cur_pos;
				Real sum_dis(0), next_dis(0), next_val(0);
				for (int idx(info.m_cur_dimIdx.first + 1); idx < to_dim; ++idx) {

				//	m_test_info.clear();
				//	m_test_info.push_back("idx\t", info.m_cur_dimIdx.first);
				//	m_test_info.push_back("pos\t", info.m_cur_dimIdx.second);
					cur_pos = getPos(idx, x[idx],
						T, info.m_sum_time, rnd);
					next_dis = getDistance(
						info.m_cur_dimIdx.first,
						info.m_cur_dimIdx.second,
						x[idx],
						info.m_cur_pos, cur_pos,
						T, info.m_sum_time, rnd);
					//m_test_info.push_back("getDistance:first", info.m_cur_dimIdx.first);
					//m_test_info.push_back("getDistance:second", info.m_cur_dimIdx.second);
					//m_test_info.push_back("getDistance:x[idx]", x[idx]);
					//m_test_info.push_back("getDistance:T", T);
					//m_test_info.push_back("getDistance:m_sum_time", info.m_sum_time);
					//m_test_info.push_back("getDistance:rnd", rnd);


					//m_test_info.push_back("from_pos", info.m_cur_pos);
					//m_test_info.push_back("cur_pos", cur_pos);
					//m_test_info.push_back("next_dis\t", next_dis);
					next_val = getPrice(
						idx, x[idx], cur_pos,
						T, info.m_sum_time, rnd);
					//	m_test_info.push_back("next_val\t", next_val);
					info.m_type_dis
						[getCorId(info.m_cur_dimIdx.first, info.m_cur_dimIdx.second)] += next_dis;
					info.m_type_val
						[getCorId(info.m_cur_dimIdx.first, info.m_cur_dimIdx.second)] += next_dis * next_val;
					info.m_cur_dimIdx.first = idx;
					info.m_cur_dimIdx.second = x[idx];
							//info.m_pos_dis.push_back(std::make_pair(x[idx], info.m_objective));
					if (flag_change_online)
						info.m_sum_time += next_dis / m_max_distance;;


				//	m_test_info.push_back("obj offset\t", next_dis * next_val);
					double origin_obj(next_dis * next_val + info.m_objective);
					
					calSPInfoObj(info, flag_correlation);
				//	m_test_info.push_back("obj\t", info.m_objective);
					std::swap(info.m_cur_pos, cur_pos);
			//		info.m_test_infos.emplace_back(m_test_info);
					if (abs(origin_obj - info.m_objective)>1e-6) {
						int stop = -1;
					}
				}
				calSPInfoObj(info, flag_correlation);
			}

			

			void addSpInfoPos(SolutionSPInfo& info,
				const std::vector<int>& x, int to_dim,
				double T,
				Random *rnd = -1,
				bool online_flag = true,
				bool corelation_flag = true
			) const{
				if (to_dim <= 0)return;
				if (info.m_cur_dimIdx.first == -1) {
					info.m_cur_dimIdx.first = 0;
					info.m_cur_dimIdx.second = x[0];
					info.m_cur_pos = getPos(
						info.m_cur_dimIdx.first,
						info.m_cur_dimIdx.second,
						T, info.m_sum_time, rnd
					);
					///info.m_pos_dis.push_back(std::make_pair(x[0], 0));
				}
				Vector cur_pos;
				//std::pair<std::string, double> info2;
				
				Real sum_dis(0), next_dis(0), next_val(0);
				for (int idx(info.m_cur_dimIdx.first + 1); idx < to_dim; ++idx) {

					//m_test_info.clear();
					//m_test_info.push_back("idx\t", info.m_cur_dimIdx.first);
					//m_test_info.push_back("pos\t", info.m_cur_dimIdx.second);
					cur_pos = getPos(idx, x[idx],
						T, info.m_sum_time, rnd);
					next_dis = getDistance(
						info.m_cur_dimIdx.first,
						info.m_cur_dimIdx.second,
						x[idx],
						info.m_cur_pos, cur_pos,
						T, info.m_sum_time, rnd);
					//m_test_info.push_back("next_dis\t", next_dis);
					next_val = getPrice(
						idx, x[idx], cur_pos,
						T, info.m_sum_time, rnd);
				//	m_test_info.push_back("next_val\t", next_val);
					info.m_type_dis
						[getCorId(info.m_cur_dimIdx.first, info.m_cur_dimIdx.second)] += next_dis;
					info.m_type_val
						[getCorId(info.m_cur_dimIdx.first, info.m_cur_dimIdx.second)] += next_dis * next_val;
					info.m_cur_dimIdx.first = idx;
					info.m_cur_dimIdx.second = x[idx];
				//	calSPInfoObj(info, corelation_flag);
				//	m_test_info.push_back("obj\t", info.m_objective);
				//	info.m_pos_dis.push_back(std::make_pair(x[idx], info.m_objective));
					if (online_flag)
						info.m_sum_time += next_dis/ m_max_distance;;
					
					std::swap(info.m_cur_pos, cur_pos);
					//info.m_test_infos.emplace_back(m_test_info);
				}


				calSPInfoObj(info, corelation_flag);
			}

			void calSPInfoObj(SolutionSPInfo& info,
				bool flag_correlation = true) const{
				info.m_objective = 0;
				if (flag_correlation) {
					for (int cor_id(0); cor_id < info.m_type_dis.size(); ++cor_id) {
						if (info.m_type_dis[cor_id] >= m_discount[cor_id].first) {
							info.m_objective += info.m_type_val[cor_id] * m_discount[cor_id].second;
						}
						else {
							info.m_objective += info.m_type_val[cor_id];
						}
					}
				}
				else {
					for (int cor_id(0); cor_id < info.m_type_dis.size(); ++cor_id) {
						info.m_objective += info.m_type_val[cor_id];
					}
				}

			}

			//	Real get_effective_mean_value(double T ,const std::vector<int>& x);
			void get_optimal(int time_omega, std::vector<std::pair<int, std::vector<int>>>& sols, Real min_dis);
		};
	}
}



#endif