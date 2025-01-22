#include "sp_map.h"

//ofec::sp::TestInfo ofec::sp::m_test_info;

void ofec::sp::Wall::setCurPosVal(double& curVal, double from, double to) const {
	double gap(to - from);
	int times(curVal / (gap * 2.0));
	curVal -= 2.0 * times * gap;
	if (curVal < 0) {
		curVal += gap * 2;
	}
	if (curVal < gap) {
		curVal += from;
	}
	else {
		curVal = to - (curVal - gap);
	}
}
void ofec::sp::Wall::setPositionStyle(Random *rnd,
	double center,
	TypeWall type,
	double max_radius,
	double small_loop_radius,
	double large_loop_inner)
{
	m_center = center;
	m_type = type;
	double& from_theta(m_angle_range[0]);
	double& to_theta(m_angle_range[1]);
	from_theta = 0;
	to_theta = 2 * OFEC_PI;
	double& inner_radius(m_radius_range[0]);
	double& outer_radius(m_radius_range[1]);
	switch (m_type)
	{
	case TypeWall::SMALL_LOOP: {
		inner_radius = 0;
		outer_radius = rnd->uniform.nextNonStd<double>(0, small_loop_radius);
		break;
	}
	case TypeWall::LARGE_LOOP: {
		outer_radius = rnd->uniform.nextNonStd<double>(large_loop_inner, max_radius);
		inner_radius = outer_radius - large_loop_inner;
		break;
	}
	case TypeWall::NEEDLE: {
		outer_radius = rnd->uniform.nextNonStd<double>(large_loop_inner + 0.1, max_radius);
		Real loop_r = rnd->uniform.nextNonStd<double>(0, outer_radius);
		inner_radius = outer_radius - loop_r;

		from_theta = rnd->uniform.nextNonStd<double>(0, 2 * OFEC_PI);
		to_theta = rnd->uniform.nextNonStd<double>(0, 2 * OFEC_PI);

		if (from_theta > to_theta) {
			std::swap(from_theta, to_theta);
		}
		//	generate_wall_points_needle(m_walls[dim_idx].m_center, m_walls[dim_idx].m_inner_radius, m_walls[dim_idx].m_outer_radius, m_walls[dim_idx].m_from_theta, m_walls[dim_idx].m_to_theta, m_walls[dim_idx].m_numPos, m_points[dim_idx]);
		break;
	}
	default: {
		inner_radius = 0;
		outer_radius = rnd->uniform.nextNonStd<double>(small_loop_radius, max_radius);
		break;
	}
	}
	m_angle_range[2] = m_angle_range[1] - m_angle_range[0];
	m_radius_range[2] = m_radius_range[1] - m_radius_range[0];
}



int ofec::sp::Wall::calNumPos(Random *rnd,
	TypeWall type,
	int num_min_points,
	int num_small_points,
	int num_large_points,
	int num_circle_points,
	int num_needle_points) {

	int numPos(0);
	switch (type)
	{
	case TypeWall::SMALL_LOOP: {
		numPos = rnd->uniform.nextNonStd<int>(num_min_points, num_small_points);
		//	generate_wall_points(m_walls[dim_idx].m_center, m_walls[dim_idx].m_inner_radius, m_walls[dim_idx].m_outer_radius, m_walls[dim_idx].m_numPos, m_points[dim_idx]);
		break;
	}
	case TypeWall::LARGE_LOOP: {
		numPos = rnd->uniform.nextNonStd<int>(num_min_points, num_large_points);
		break;
	}
	case TypeWall::NEEDLE: {
		numPos = rnd->uniform.nextNonStd<int>(num_min_points, num_needle_points);
		break;
	}
	default: {
		numPos = rnd->uniform.nextNonStd<int>(num_min_points, num_circle_points);
		break;
	}
	}
	return numPos;
}

void ofec::sp::Wall::setNumPos(Random *rnd,
	int num_min_points,
	int num_small_points,
	int num_large_points,
	int num_circle_points,
	int num_needle_points)
{
	m_numPos = calNumPos(rnd, m_type, num_min_points,
		num_small_points, num_large_points, num_circle_points, num_needle_points);
}



void ofec::sp::MapInfo::calMaxDistance()
{
	m_max_distance = 0;
	double maxScale(0);
	m_dim_distance.resize(m_wall.size());
	for (int dim(0); dim + 1 < m_wall.size(); ++dim) {
		maxScale = 0;
		for (auto& dimIter : m_edges[dim]) {
			for (auto& edgeIter : dimIter) {
				maxScale = std::max(maxScale, edgeIter.m_distance_distr.getMaxVal());
			}
		}
	//	m_dim_distance[]
		// for test
	//	auto pos1(m_wall[dim].outerPosition((dim % 2) * OFEC_PI));
	//	auto pos2(m_wall[dim + 1].outerPosition(((dim + 1) % 2) * OFEC_PI));
		m_max_distance += maxScale *
			m_wall[dim].outerPosition((dim % 2) * OFEC_PI)
			.distance(m_wall[dim + 1].outerPosition(((dim + 1) % 2) * OFEC_PI));
	}
}

void  ofec::sp::MapInfo::initializeWallType(Random *rnd, const ParaInfo& para) {

	std::vector<TypeWall> random_needle(
		para.m_static_info.m_wall_num[static_cast<int>(TypeWall::CIRCLE)]
		+ para.m_static_info.m_wall_num[static_cast<int>(TypeWall::NEEDLE)], TypeWall::CIRCLE);
	for (int i(0); i < para.m_static_info.m_wall_num[static_cast<int>(TypeWall::NEEDLE)]; ++i) {
		random_needle[i] = TypeWall::NEEDLE;
	}

	rnd->uniform.shuffle(random_needle.begin(), random_needle.end());
	std::vector<TypeWall> modility(
		para.m_static_info.m_wall_num[static_cast<int>(TypeWall::SMALL_LOOP)]
		+ para.m_static_info.m_wall_num[static_cast<int>(TypeWall::LARGE_LOOP)], TypeWall::LARGE_LOOP);
	for (int i(0); i < para.m_static_info.m_wall_num[static_cast<int>(TypeWall::LARGE_LOOP)]; ++i) {
		modility[i] = TypeWall::LARGE_LOOP;
	}
	rnd->uniform.shuffle(modility.begin(), modility.end());
	std::vector<TypeWall> type;

	int type_idx(rnd->uniform.nextNonStd<int>(0, random_needle.size()));
	type.insert(type.begin(), random_needle.begin(), random_needle.begin() + type_idx);
	type.insert(type.begin() + type.size(), modility.begin(), modility.begin() + modility.size());
	//type.emplace_back(modility);
	type.insert(type.begin() + type.size(), random_needle.begin() + type_idx, random_needle.begin() + random_needle.size());

	int variable_size(type.size());
	m_wall.resize(variable_size);
	for (int dim_idx(0); dim_idx < variable_size; ++dim_idx) {
		Real wall_center(0);
		if (dim_idx) {
			wall_center = m_wall[dim_idx - 1].getCenter() + rnd->uniform.nextNonStd<Real>(para.m_static_info.m_width_range.first, para.m_static_info.m_width_range.second);
		}
		m_wall[dim_idx].setPositionStyle(rnd, wall_center, type[dim_idx], para.m_static_info.m_max_radius, para.m_static_info.m_small_loop_radius, para.m_static_info.m_large_loop_inner);
		m_wall[dim_idx].setNumPos(rnd, para.m_static_info.m_num_min_points, para.m_static_info.m_num_small_points, para.m_static_info.m_num_large_points, para.m_static_info.m_num_circle_points, para.m_static_info.m_num_needle_points);
		//	m_wall[dim_idx].m_center_posId = rnd->uniform.nextNonStd<int>(0, m_wall_types[dim_idx].m_numPos);
		//	m_wall_types[dim_idx].m_noisy_ratio = para.m_static_info.m_wallinfo_noisy_ratio;

	}
}


void  ofec::sp::MapInfo::initialize(Random *rnd, ParaInfo& para) {

	//para.m_fun_para.initialize(para.m_static_info.m_T);
	para.initialize();
	initializeWallType(rnd, para);
	m_discount.resize(para.m_static_info.m_station_type_num);
	for (auto& it : m_discount) {
		it.first = rnd->uniform.nextNonStd<Real>(para.m_static_info.m_discount_level.first, para.m_static_info.m_discount_level.second);
		it.second = rnd->uniform.nextNonStd<Real>(para.m_static_info.m_discount_ratio.first, para.m_static_info.m_discount_ratio.second);
	}

	m_stations.resize(m_wall.size());
	for (int wall_id(0); wall_id < m_wall.size(); ++wall_id) {
		m_stations[wall_id].resize(m_wall[wall_id].numPos());

	}
	m_edges.resize(m_wall.size());
	for (size_t wall_id(0); wall_id + 1 < m_wall.size(); ++wall_id) {
		m_edges[wall_id].resize(m_wall[wall_id].numPos());
		for (int from_id(0); from_id < m_wall[wall_id].numPos(); ++from_id) {
			m_edges[wall_id][from_id].resize(m_wall[wall_id + 1].numPos());
		}
	}

	updatePosInfo(rnd);
	initializeDynamicInfo(rnd, para);
	calMaxDistance();
}

void ofec::sp::MapInfo::transferXfactors(std::vector<int>& x, double T)const {
	for (int idx(0); idx < numWall(); ++idx) {
		if (getScenesDim(T)[idx]) {
			x[idx] = getBasicSols(T).variable()[idx];
		}
	}
}




void ofec::sp::MapInfo::transferXEnvNoisy(
	std::vector<int>& x,
	double T,
	Random *rnd , bool flag_env_noisy ) const {

	transferXNoisy(m_env_noisy_ratio, true, getScenesNumDims(T), x, T, rnd, flag_env_noisy);
	/*if (flag_env_noisy&&rnd!=-1&&
		rnd->uniform.next() < m_env_noisy_ratio) {
		int change_dim(0);
		change_dim = rnd->uniform.nextNonStd<int>(0,);
		for (int idx(0); idx < numWall(); ++idx) {
			if (getScenesDim(T)[idx]) {
				if (--change_dim == 0) {
					x[idx] = rnd->uniform.nextNonStd<int>(0, numCandidates(idx));
					break;
				}
			}
		}
	}*/
}


void ofec::sp::MapInfo::transferXNoisy(
	double noisy_ratio,
	bool flag_scene,
	int numDims,
	std::vector<int>& x,
	double T,
	Random *rnd ,
	bool flag_noisy )const {

	if (flag_noisy && rnd != -1 &&
		rnd->uniform.next() < noisy_ratio) {
		int change_dim(0);
		change_dim = rnd->uniform.nextNonStd<int>(0, numDims);
		for (int idx(0); idx < numWall(); ++idx) {
			if (getScenesDim(T)[idx]==flag_scene) {
				if (change_dim-- == 0) {
					x[idx] = rnd->uniform.nextNonStd<int>(0, numCandidates(idx));
					break;
				}
			}
		}
	}
}
void ofec::sp::MapInfo::transferXVarNoisy(
	std::vector<int>& x,
	double T,
	Random *rnd , bool flag_var_noisy )const {

//	transferXNoisy(m_var_noisy_ratio,)

	transferXNoisy(m_var_noisy_ratio, false, 
		numWall() - getScenesNumDims(T), x, T, rnd, flag_var_noisy);

	//if (flag_var_noisy && rnd != -1 &&
	//	rnd->uniform.next() < m_var_noisy_ratio) {
	//	int change_dim(0);
	//	change_dim = rnd->uniform.nextNonStd<int>(0, numWall() - getScenesNumDims(T));
	//	for (int idx(0); idx < numWall(); ++idx) {
	//		if (!getScenesDim(T)[idx]) {
	//			if (--change_dim == 0) {
	//				x[idx] = rnd->uniform.nextNonStd<int>(0, numCandidates(idx));
	//				break;
	//			}
	//		}
	//	}
	//}
}


void ofec::sp::MapInfo::transferXWorkingSol(
	std::vector<int>& x,
	const std::vector<int>& working_x)const {
	for (int idx(0); idx < working_x.size(); ++idx) {
		x[idx] = working_x[idx];
	}
}


void ofec::sp::MapInfo::updatePosInfo(Random *rnd) {
	for (int id_wall(0); id_wall < m_wall.size(); ++id_wall) {
		m_wall[id_wall].getPosInfo()[static_cast<int>(PosInfoIdx::WALL_IDX)] = id_wall;
		for (int id_pos(0); id_pos < m_stations[id_wall].size(); ++id_pos) {
			m_stations[id_wall][id_pos].m_pos_info[static_cast<int>(PosInfoIdx::WALL_IDX)] = id_wall;
			m_stations[id_wall][id_pos].m_pos_info[static_cast<int>(PosInfoIdx::POS_FROM_IDX)] = id_pos;
			m_stations[id_wall][id_pos].m_pos_info[static_cast<int>(PosInfoIdx::POS_TO_IDX)] = -1;
			m_stations[id_wall][id_pos].m_pos_info[static_cast<int>(PosInfoIdx::COR_IDX)] =
				rnd->uniform.nextNonStd<int>(0, m_discount.size());
			if (id_wall + 1 < m_wall.size()) {
				for (int id_pos_to(0); id_pos_to < m_stations[id_wall + 1].size(); ++id_pos_to) {
					m_edges[id_wall][id_pos][id_pos_to].
						m_pos_info[static_cast<int>(PosInfoIdx::WALL_IDX)]
						= id_wall;
					m_edges[id_wall][id_pos][id_pos_to].
						m_pos_info[static_cast<int>(PosInfoIdx::POS_FROM_IDX)]
						= id_pos;
					m_edges[id_wall][id_pos][id_pos_to].
						m_pos_info[static_cast<int>(PosInfoIdx::POS_TO_IDX)]
						= id_pos_to;
				}
			}

		}
	}

}
void ofec::sp::MapInfo::initializeDynamicInfo(Random *rnd, ParaInfo& para) {
	for (auto& it : m_wall) {
		it.initDynamicIdxs(rnd, para);
	}
	para.m_dynamic_info.initFeasiblePar(rnd);
	for (auto& it : m_stations) {
		for (auto& it2 : it) {
			para.m_dynamic_info.initPar(rnd);
			it2.m_feasible_distr.initialize(rnd, para.m_dynamic_info.m_par);
		}
	}

	para.m_dynamic_info.initPriceVarPar(rnd);
	for (auto& it : m_stations) {
		for (auto& it2 : it) {
			para.m_dynamic_info.initPar(rnd);
			it2.m_price_distr.initialize(rnd, para.m_dynamic_info.m_par);
		}
	}
	//std::array<double, 3> radius = {0.3,2.0,1.7};
	//std::array<double, 3> angle = { 0,2 * OFEC_PI, 2.0*OFEC_PI  };
	for (int id_wall(0); id_wall < m_wall.size(); ++id_wall) {
		//para.m_dynamic_info.initDynamicVarStation(rnd, radius, true);
		//para.m_dynamic_info.initFunTestRadius(rnd);
		para.m_dynamic_info.initDynamicVarStation(rnd, m_wall[id_wall].radiusRange(), true);
		for (auto& it : m_stations[id_wall]) {
			para.m_dynamic_info.initPar(rnd);
			it.m_radius_distr.initialize(rnd, para.m_dynamic_info.m_par);
		}
		//para.m_dynamic_info.initDynamicVarStation(rnd, angle, false);
		//para.m_dynamic_info.initFunTestAngle(rnd);
		para.m_dynamic_info.initDynamicVarStation(rnd, m_wall[id_wall].angleRange(), false);
		for (auto& it : m_stations[id_wall]) {
			para.m_dynamic_info.initPar(rnd);
			it.m_angle_distr.initialize(rnd, para.m_dynamic_info.m_par);
		}
	}
	para.m_dynamic_info.initFeasiblePar(rnd);
	for (auto& it : m_edges) {
		for (auto& it2 : it) {
			for (auto& it3 : it2) {
				para.m_dynamic_info.initPar(rnd);
				it3.m_feasible_distr.initialize(rnd, para.m_dynamic_info.m_par);
			}
		}
	}

	para.m_dynamic_info.initRunningTimeVarPar(rnd);
	for (auto& it : m_edges) {
		for (auto& it2 : it) {
			for (auto& it3 : it2) {
				para.m_dynamic_info.initPar(rnd);
				it3.m_distance_distr.initialize(rnd, para.m_dynamic_info.m_par);
			}
		}
	}




	para.m_dynamic_info.initMeshSeedFun(rnd);
	for (auto& it : m_wall) {
		para.m_dynamic_info.initPar(rnd);
		it.m_mesh_price.initialize(rnd, para.m_dynamic_info.m_par);
	}

	initEnvDynamics(rnd, para); 
}

void  ofec::sp::MapInfo::get_optimal(int time_omega, std::vector<std::pair<int, std::vector<int>>>& sols, Real min_dis) {

	// as time is infinit, we can not find the optimal solutions by dp
//
//	// use dp to find the optimal solutions;
//	std::vector<solution_type> sols;
//	std::vector<std::vector<solution_type>> dp_time_pos_sols;
//	std::vector<std::vector<solution_type>> dp_time_pos_after_sols;
//	solution_type cur_sol;
//	cur_sol.resize_objective(2);
//	cur_sol.variable().vect().clear();
////	dp_sols.push_back(std::vector<solution_type>());
//	std::vector<solution_type> v_cur_sol;
	//
	//{
	//	int variable_size(m_wall_types.size());
	//	
	//	int dim(0);
	////	dp_sols.resize(m_wall_types[dim].get_numPos());
	//	for (int idx(0); idx < m_wall_types[dim].get_numPos(); ++idx) {
	//		dp_sols[idx].resize(1);
	//		dp_sols[idx].front().variable().vect().push_back(dim);
	//		dp_sols[idx].front().resize_objective(2);
	//	}

	//	for (++dim; dim < variable_size; ++dim) {
	//		Real cur_val(0);
	//	//	dp_after_sols.resize(m_wall_types[dim].get_numPos());
	//		for (int idx(0); idx < m_wall_types[dim].get_numPos(); ++idx) {
	//			
	//		}
	//	}
	//}
	//

}