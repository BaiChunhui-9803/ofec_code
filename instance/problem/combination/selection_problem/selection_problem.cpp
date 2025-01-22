#include"selection_problem.h"
using namespace ofec;
const int OFEC::SP::wall_type_num = 4;
const Real OFEC::SP::SelectionProblem::m_eps = 1e-6;

void OFEC::SP::wall_info::set_position_style(Random *rnd,Real* center, wall_type* type, Real max_radius, Real small_loop_radius, Real large_loop_inner)
{
	if (center != nullptr)
		m_center = *center;
	if (type != nullptr)
		m_type = *type;

	switch (m_type)
	{
	case wall_type::small_loop: {
		m_inner_radius = 0;
		m_outer_radius = rnd->uniform.nextNonStd<Real>(0, small_loop_radius);
		break;
	}
	case wall_type::large_loop: {
		
		m_outer_radius = rnd->uniform.nextNonStd<Real>(large_loop_inner, max_radius);
		m_inner_radius = m_outer_radius - large_loop_inner;
		break;
	}

	case wall_type::needle: {
		m_outer_radius = rnd->uniform.nextNonStd<Real>(large_loop_inner + 0.1, max_radius);
		Real loop_r = rnd->uniform.nextNonStd<Real>(0, m_outer_radius);
		m_inner_radius = m_outer_radius - loop_r;

		m_from_theta = rnd->uniform.nextNonStd<Real>(0, 2 * OFEC_PI);
		m_to_theta = rnd->uniform.nextNonStd<Real>(0, 2 * OFEC_PI);

		if (m_from_theta > m_to_theta) {
			std::swap(m_from_theta, m_to_theta);
		}
		//	generate_wall_points_needle(m_walls[dim_idx].m_center, m_walls[dim_idx].m_inner_radius, m_walls[dim_idx].m_outer_radius, m_walls[dim_idx].m_from_theta, m_walls[dim_idx].m_to_theta, m_walls[dim_idx].m_numPos, m_points[dim_idx]);
		break;
	}
	default: {
		m_inner_radius = 0;
		m_outer_radius = rnd->uniform.nextNonStd<Real>(small_loop_radius, max_radius);
		break;
	}

	}
}



int OFEC::SP::wall_info::cal_numPos(Random *rnd,wall_type type, int num_small_points,
	int num_large_points,
	int num_circle_points,
	int num_needle_points) {

	int numPos(0);
	switch (type)
	{
	case wall_type::small_loop: {
		numPos = rnd->uniform.nextNonStd<int>(5, num_small_points);
		//	generate_wall_points(m_walls[dim_idx].m_center, m_walls[dim_idx].m_inner_radius, m_walls[dim_idx].m_outer_radius, m_walls[dim_idx].m_numPos, m_points[dim_idx]);
		break;
	}
	case wall_type::large_loop: {
		numPos = rnd->uniform.nextNonStd<int>(5, num_large_points);
		break;
	}

	case wall_type::needle: {
		numPos = rnd->uniform.nextNonStd<int>(5, num_needle_points);
		break;
	}
	default: {
		numPos = rnd->uniform.nextNonStd<int>(5, num_circle_points);
		break;
	}
	}
	return numPos;
}

void OFEC::SP::wall_info::set_numPos(Random *rnd,int num_small_points, int num_large_points, int num_circle_points, int num_needle_points)
{

	m_numPos = cal_numPos(rnd,m_type, num_small_points, num_large_points, num_circle_points, num_needle_points);
}



void  OFEC::SP::map_info::initialize_wall_type(Random *rnd) {
	//m_points.clear();

	//problem::resize(variable_size);
	std::vector<wall_type> random_needle(m_parameter.m_wall_num[static_cast<int>(wall_type::circle)] + m_parameter.m_wall_num[static_cast<int>(wall_type::needle)], wall_type::circle);
	for (int i(0); i < m_parameter.m_wall_num[static_cast<int>(wall_type::needle)]; ++i) {
		random_needle[i] = wall_type::needle;
	}
	
	rnd->uniform.shuffle(random_needle.begin(), random_needle.end());
	std::vector<wall_type> modility(m_parameter.m_wall_num[static_cast<int>(wall_type::small_loop)] + m_parameter.m_wall_num[static_cast<int>(wall_type::large_loop)], wall_type::small_loop);
	for (int i(0); i < m_parameter.m_wall_num[static_cast<int>(wall_type::large_loop)]; ++i) {
		modility[i] = wall_type::large_loop;
	}
	rnd->uniform.shuffle(modility.begin(), modility.end());
	std::vector<wall_type> type;

	int type_idx(rnd->uniform.nextNonStd<int>(0, random_needle.size()));
	type.insert(type.begin(), random_needle.begin(), random_needle.begin() + type_idx);
	type.insert(type.begin() + type.size(), modility.begin(), modility.begin() + modility.size());
	//type.emplace_back(modility);
	type.insert(type.begin() + type.size(), random_needle.begin() + type_idx, random_needle.begin() + random_needle.size());

	int variable_size(type.size());
	m_wall_types.resize(variable_size);
	for (int dim_idx(0); dim_idx < variable_size; ++dim_idx) {
		Real wall_center(0);
		if (dim_idx) {
			wall_center = m_wall_types[dim_idx - 1].m_center + rnd->uniform.nextNonStd<Real>(m_parameter.m_width_range.first, m_parameter.m_width_range.second);
		}
		m_wall_types[dim_idx].set_position_style(rnd,&wall_center, &type[dim_idx], m_parameter.m_max_radius, m_parameter.m_small_loop_radius, m_parameter.m_large_loop_inner);
		m_wall_types[dim_idx].set_numPos(rnd,m_parameter.m_num_small_points, m_parameter.m_num_large_points, m_parameter.m_num_circle_points, m_parameter.m_num_needle_points);
		m_wall_types[dim_idx].m_bag_value_seed = rnd->uniform.next();
		m_wall_types[dim_idx].m_center_id = rnd->uniform.nextNonStd<int>(0, m_wall_types[dim_idx].m_numPos);
		m_parameter.update_price_basic_par();
		m_wall_types[dim_idx].m_bag_basic_value.initialize(m_parameter.m_price_basic_par);
	}
}


void  OFEC::SP::map_info::initialize(Random *rnd) {

	m_parameter.initialize_edge_par();
	m_parameter.initialize_feasible_par();
	m_parameter.initialize_price_basic_par();
	m_parameter.initialize_station_basic_par();
	

	initialize_wall_type(rnd);
	m_discount.resize(m_parameter.m_station_type_num);
	for (auto& it: m_discount) {
		it.first = rnd->uniform.nextNonStd<Real>(m_parameter.m_ratio_level.first, m_parameter.m_ratio_level.second);
		it.second= rnd->uniform.nextNonStd<Real>(m_parameter.m_discount.first, m_parameter.m_discount.second);
	}

	m_stations.resize(m_wall_types.size());
	for (int wall_id(0); wall_id < m_wall_types.size(); ++wall_id) {
		m_stations[wall_id].resize(m_wall_types[wall_id].m_numPos);
		for (int pos_id(0); pos_id < m_wall_types[wall_id].m_numPos; ++pos_id) {
			auto& cur_station(m_stations[wall_id][pos_id]);

			cur_station.m_wall = &m_wall_types[wall_id];
			cur_station.m_station_wall_id = pos_id;
			cur_station.m_cor_id = rnd->uniform.nextNonStd<int>(0, m_parameter.m_station_type_num);
			
			cur_station.m_static_flag =  rnd->uniform.next() > m_parameter.m_moving_station_ratio;

			m_parameter.update_station_basic_par(cur_station.m_static_flag);
			cur_station.m_radius.initialize(m_parameter.m_station_basic_par);
			m_parameter.update_station_basic_par(cur_station.m_static_flag);
			cur_station.m_angle.initialize(m_parameter.m_station_basic_par);
			m_parameter.update_feasible_par();
			cur_station.m_feasible.initialize(m_parameter.m_feasible_par);
			
		}
	}


	m_edges.resize(m_wall_types.size());
	for (size_t wall_id(0); wall_id+1 < m_wall_types.size(); ++wall_id) {
		m_edges[wall_id].resize(m_wall_types[wall_id].m_numPos);
		for (int from_id(0); from_id < m_wall_types[wall_id].m_numPos; ++from_id) {
			m_edges[wall_id][from_id].resize(m_wall_types[wall_id+1].m_numPos);
			for (int to_id(0); to_id < m_wall_types[wall_id+1].m_numPos; ++to_id) {
				m_edges[wall_id][from_id][to_id].m_station_from = from_id;
				m_edges[wall_id][from_id][to_id].m_station_to = to_id;
				m_parameter.update_edge_par();
				m_edges[wall_id][from_id][to_id].m_distance_t_scale.initialize(m_parameter.m_edge_par);

				m_parameter.update_feasible_par();
				m_edges[wall_id][from_id][to_id].m_feasible.initialize(m_parameter.m_feasible_par);
			}
		}
	}
	m_max_distance = 0;
	for (int dim(0); dim + 1 < m_wall_types.size(); ++dim) {
		m_max_distance += m_wall_types[dim].get_outer_position((dim % 2) * OFEC_PI).distance(m_wall_types[dim + 1].get_outer_position(((dim + 1) % 2) * OFEC_PI));
	}
	m_max_distance *= (1.0 + m_parameter.m_edge_noisy_time.second);
}


Real OFEC::SP::map_info::get_effective_mean_value(int T,const std::vector<int>& x) {
	if (x.size() != m_wall_types.size()) {
		return 0.0;
	}
	else {
		std::vector<int> nei(x);
		bool feasible_tag = true;
		Real x_obj(get_mean_value(T,x,feasible_tag));
		Real y_obj(0);
		Real eff_val(0);
		int variable_size(m_wall_types.size());
		for (int dim(0); dim < variable_size; ++dim) {
			Real cur_val(0);
			for (int idx(0); idx < m_wall_types[dim].get_numPos(); ++idx) {
				if (idx == x[dim]) {
					cur_val += x_obj;
				}
				else {
					nei[dim] = idx;
					cur_val += get_mean_value(T, nei, feasible_tag);

				}
			}
			nei[dim] = x[dim];
			eff_val += cur_val / static_cast<double>(m_wall_types[dim].get_numPos());
		}
		eff_val /= static_cast<double>(variable_size);
		return eff_val;
	}
}


void  OFEC::SP::map_info::get_optimal(int time_omega, std::vector<std::pair<int, std::vector<int>>>& sols,Real min_dis) {

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

bool OFEC::SP::SelectionProblem::same(const SolutionBase& s1, const SolutionBase& s2) const
{
	auto& x1 = dynamic_cast<const Solution<VariableVector<int>>&>(s1).variable();
	auto& x2 = dynamic_cast<const Solution<VariableVector<int>>&>(s2).variable();
	return x1.vect() == x2.vect();
}

Real OFEC::SP::SelectionProblem::variableDistance(const SolutionBase& s1, const SolutionBase& s2) const
{
	auto& x1 = dynamic_cast<const Solution<VariableVector<int>>&>(s1).variable();
	auto& x2 = dynamic_cast<const Solution<VariableVector<int>>&>(s2).variable();

	if (x1.size() != x2.size())return x1.size();

	Real dis(0);
	Real maxDis(m_number_variables * m_number_variables);

	for (int i(0); i < m_number_variables; ++i) {
		if (x1[i] == x2[i]) {
			dis += 1;
			if (i + 1 < m_number_variables && x1[i + 1] == x2[i + 1]) {
				dis += m_number_variables;
			}
		}
	}
	return (maxDis - dis) / maxDis;
}

Real OFEC::SP::SelectionProblem::variableDistance(const VariableBase& s1, const VariableBase& s2) const
{
	auto& x1 = dynamic_cast<const variable_type&>(s1).vect();
	auto& x2 = dynamic_cast<const variable_type&>(s2).vect();



	if (x1.size() != x2.size())return x1.size();

	Real dis(0);
	Real maxDis(m_number_variables * m_number_variables);
	for (int i(0); i < m_number_variables; ++i) {
		if (x1[i] == x2[i]) {
			dis += 1;
			if (i + 1 < m_number_variables && x1[i + 1] == x2[i + 1]) {
				dis += m_number_variables;
			}
		}
	}
	return (maxDis - dis) / maxDis;
}



void OFEC::SP::SelectionProblem::evaluate_(SolutionBase& s, bool effective) {
	auto& x = dynamic_cast<solution_type&>(s).variable();
	auto& obj = dynamic_cast<solution_type&>(s).objective();
	bool feasible_flag(false);
	//obj.resize(1);
	obj[0] = m_map.get_value(m_T, x.vect(), feasible_flag);
}
//
//evaluation_tag OFEC::SP::SelectionProblem::evaluate_(solution_base& s, caller call, bool effective, bool initialized)
//{
//
//	auto& x = dynamic_cast<solution_type&>(s).variable();
//	auto& obj = dynamic_cast<solution_type&>(s).objective();
//	bool feasible_flag(false);
//	obj.resize(1);
//	obj[0] = m_map.get_value(m_T, x.vect(), feasible_flag);
//
//	if (initialized) {
//		if (effective)		m_evaluations++;
//		if (call == caller::Algorithm && global::ms_global->m_algorithm && global::ms_global->m_algorithm->terminating())
//			return evaluation_tag::Terminate;
//		if (is_change())change();
//	}
//	return evaluation_tag::Normal;
//}

//void OFEC::SP::SelectionProblem::initialize_solution(solution_base& s) const
//{
//	auto& x = dynamic_cast<solution_type&>(s).variable();
//	x.resize(m_number_variables);
//	for (int i(0); i < x.size(); ++i) {
//		x[i] = rnd->uniform.nextNonStd<int>(0, m_map.m_wall_types[i].m_numPos);
//	}
//	
//}


void OFEC::SP::SelectionProblem::initializeSolution(SolutionBase& s, Random *rnd) const  {

	auto& x = dynamic_cast<solution_type&>(s).variable();
	x.resize(m_number_variables);
	for (int i(0); i < x.size(); ++i) {
		x[i] =rnd->uniform.nextNonStd<int>(0, m_map.m_wall_types[i].m_numPos);
		//rnd.
	//	x[i] = rnd->uniform.nextNonStandard<int>(0, m_map.m_wall_types[i].m_numPos);
	}

}

void OFEC::SP::SelectionProblem::initialize_()
{
	Problem::initialize_();
	Noisy::initialize_();
	Dynamic::initialize_();
	/*
	* 		m_total_evaluations = 0; 
		m_number_objectives = 0;
		m_number_constraints = 0;
		m_objective_accuracy = -1;
		m_variable_accuracy = -1;
		m_variable_niche_radius = -1;
		m_optimize_mode.clear();
		m_constraint.clear();
		m_params.clear();
		m_initialized = false;
	*/
	addTag(ProblemTag::DOP);
	addTag(ProblemTag::NoisyOP);


	m_T = 0;
	m_map.initialize(m_random.get());
	m_number_variables = m_map.m_wall_types.size();

	m_number_objectives = 1;
	m_number_constraints = 0;
	//m_frequency = 1000;
	m_optimize_mode.resize(m_number_objectives);
	m_optimize_mode.front()=OptimizeMode::kMinimize;

}

void OFEC::SP::SelectionProblem::copy(const Problem& rP)
{
}


//void OFEC::SP::SelectionProblem::generate_HLS_samples(const SolutionBase& s, std::vector<std::unique_ptr<SolutionBase>>& samples, Random *rnd) {
//	auto& x = dynamic_cast<const solution_type&>(s).variable();
//	auto& obj = dynamic_cast<const solution_type&>(s).objective();
//	int sample_id(0);
//	std::vector<int> dim_sample(m_number_variables,0);
//	generate_samples_numbers(dim_sample, samples.size(), rnd);
//	std::vector<int> pos_sample;
//	for (int dim(0); dim < m_number_variables; ++dim) {
//		pos_sample.resize(m_map.m_wall_types[dim].get_numPos());
//		generate_samples_numbers(pos_sample, dim_sample[dim], rnd);
//		for (int pos_idx(0); pos_idx < pos_sample.size(); ++pos_idx) {
//			while (pos_sample[pos_idx]--) {
//				samples[sample_id].reset(new solution_type());
//				solution_type& cur(dynamic_cast<solution_type&>(*samples[sample_id]));
//				cur.variable() = x;
//				cur.objective().resize(1);
//				cur.variable()[dim] = pos_idx;
//				++sample_id;
//			}
//		}
//	}
//
//	//std::vector<int> dim_idxs(m_number_variables);
//	//for (int idx(0); idx < dim_idxs.size(); ++idx) {
//	//	dim_idxs[idx] = idx;
//	//}
//	//rnd->uniform.shuffle(dim_idxs.begin(), dim_idxs.end());
//
//	//int sample_num(static_cast<double>(samples.size()) / m_number_variables);
//	//int sample_left(samples.size() - sample_num * m_number_variables);
//	//int cur_sample(0), cur_dim(0);
//	//for (int dim_idx(0); dim_idx < m_number_variables; ++dim_idx) {
//	//	cur_sample = sample_num;
//	//	if (dim_idx < sample_left) ++cur_sample;
//	//	if (cur_sample) {
//	//		std::vector<int> pos_idxs(m_map.m_wall_types[dim_idx].get_numPos());
//	//		for (int idx(0); idx < pos_idxs.size(); ++idx) {
//	//			pos_idxs[idx] = idx;
//	//		}
//	//		rnd->uniform.shuffle(pos_idxs.begin(), pos_idxs.end());
//
//	//		int inner_sample_num(cur_sample/ m_map.m_wall_types[dim_idx].get_numPos());
//	//		int inner_sample_left(samples.size() - inner_sample_num * m_map.m_wall_types[dim_idx].get_numPos());
//
//	//		for (int pos_idx(0); pos_idx < m_number_variables; ++pos_idx) {
//	//			cur_sample = sample_num;
//	//			if (dim_idx < sample_left) ++cur_sample;
//	//			if (cur_sample)
//	//		}
//	//	}
//	//	else break;
//	//}
//}




Real OFEC::SP::SelectionProblem::evaluate_mean_value(SolutionBase& s)
{
	auto& x = dynamic_cast<solution_type&>(s).variable();
	bool feasible_flag(false);
	return m_map.get_mean_value(m_T, x.vect(), feasible_flag);
}
void OFEC::SP::SelectionProblem::generate_samples_numbers(std::vector<int>& idx_samples, int sample_num, Random *rnd) {
	int size(idx_samples.size());
	std::vector<int> idxs(size);
	for (int idx(0); idx < idxs.size(); ++idx) {
		idxs[idx] = idx;
	}
	rnd->uniform.shuffle(idxs.begin(), idxs.end());

	int s_size(static_cast<double>(sample_num) / size);
	int s_left(sample_num - s_size * size);
	int cur_sample(0), cur_dim(0);
	for (int idx(0); idx < size; ++idx) {
		cur_sample = s_size;
		if (idx < s_left) ++cur_sample;
		idx_samples[idxs[idx]] = cur_sample;
	}
}
void OFEC::SP::SelectionProblem::showInfomations(Algorithm *alg)
{
	if (!alg->candidates().empty()) {
		printfSolution(alg, *alg->candidates().front());
	}
}
void OFEC::SP::SelectionProblem::printfSolution(Algorithm *alg,const SolutionBase& sol)
{
	auto& x = dynamic_cast<const solution_type&>(sol).variable();
	std::cout << "effective_eval \t" << alg->evaluations()
		<< "\t objective value\t" << sol.objective()[0]
		<< "\t true objective value\t" << std::endl;
	std::cout << "solution encoding\t";
	for (auto& it : x.vect()) std::cout << it << "\t";
	std::cout << std::endl;
}
std::vector<int> OFEC::SP::SelectionProblem::generate_samples_numbers(int size, int sample_num, Random *rnd)
{
	std::vector<int> idxs(size);
	for (int idx(0); idx < idxs.size(); ++idx) {
		idxs[idx] = idx;
	}
	rnd->uniform.shuffle(idxs.begin(), idxs.end());

	int s_size(static_cast<double>(sample_num) / size);
	int s_left(sample_num - s_size * size);
	int cur_sample(0), cur_dim(0);
	std::vector<int> each_num(size, 0);
	for (int idx(0); idx < size; ++idx) {
		cur_sample = s_size;
		if (idx < s_left) ++cur_sample;
		each_num[idxs[idx]] = cur_sample ;
	}
	return std::move(each_num);
}

void OFEC::SP::SelectionProblem::generate_level_samples(std::vector<solution_type>& samples, Random *rnd)
{
	samples.clear();
	int max_nodes(0);
	std::vector<std::vector<int>> solus(m_map.m_wall_types.size());

	int variable_size(m_map.m_wall_types.size());
	for (int dim(0); dim < m_map.m_wall_types.size(); ++dim) {
		solus[dim].resize(m_map.m_wall_types[dim].m_numPos);
		max_nodes = std::max(max_nodes, m_map.m_wall_types[dim].m_numPos);
		for (int idx(0); idx < solus[dim].size(); ++idx) {
			solus[dim][idx] = idx;
		}
		rnd->uniform.shuffle(solus[dim].begin(), solus[dim].end());
		//global::ms_global->m_uniform[call]->shuffle(solus[dim].begin(), solus[dim].end());	
	}

	samples.resize(max_nodes);
	for (auto& it : samples) {
		it.variable().resize(variable_size);
		it.objective().resize(1);
	}

	for (int idx(0); idx < max_nodes; ++idx) {
		for (int dim(0); dim < variable_size; ++dim) {
			samples[idx].variable()[dim] = solus[dim][idx % m_map.m_wall_types[dim].m_numPos];
		}
	}
}


Real OFEC::SP::SelectionProblem::evaluate_effective_value(SolutionBase& s) {
	solution_type neighbors(dynamic_cast<const solution_type&>(s));
	solution_type cur_sol(dynamic_cast<const solution_type&>(s));
	Real eff_val(0);
	int variable_size(m_map.m_wall_types.size());
	Real cur_sol_val(evaluate_mean_value(cur_sol));
	for (int dim(0); dim < variable_size; ++dim) {
		Real cur_val(0);
		for (int idx(0); idx < m_map.m_wall_types[dim].get_numPos(); ++idx) {
			if (idx == cur_sol.variable()[dim]) {
				cur_val += cur_sol_val;
			}
			else {
				neighbors.variable()[dim] = idx;
				cur_val += evaluate_mean_value(neighbors);
			}
		}
		neighbors.variable()[dim] = cur_sol.variable()[dim];
		eff_val+=cur_val / static_cast<double>(m_map.m_wall_types[dim].get_numPos());
	}
	eff_val /= static_cast<double>(variable_size);
	return eff_val;
}

//evaluation_tag OFEC::SP::SelectionProblem::evaluate_effective_value(solution_base& s, caller call, bool effective, bool initialized)
//{
//	evaluation_tag rf_tag(evaluate_(s, call, effective, initialized));
//	if (rf_tag == evaluation_tag::Terminate) return rf_tag;
//	solution_type neighbors(dynamic_cast<const solution_type&>(s));
//	solution_type cur_sol(dynamic_cast<const solution_type&>(s));
//	Real eff_val(0);
//	int variable_size(m_map.m_wall_types.size());
//	for (int dim(0); dim < variable_size; ++dim) {
//		Real cur_val(0);
//		for (int idx(0); idx < m_map.m_wall_types[dim].get_numPos(); ++idx) {
//			if (idx == cur_sol.variable()[dim]) {
//				cur_val += cur_sol.objective()[0];
//			}
//			else {
//				neighbors.variable()[dim] = idx;
//				rf_tag=evaluate_(neighbors, call, effective, initialized);
//				if (rf_tag == evaluation_tag::Terminate) return rf_tag;
//				cur_val += neighbors.objective()[0];
//
//			}
//		}
//		neighbors.variable()[dim] = cur_sol.variable()[dim];
//		eff_val+=cur_val / static_cast<double>(m_map.m_wall_types[dim].get_numPos());
//	}
//	eff_val /= static_cast<double>(variable_size);
//	cur_sol.constraint_value().front() = eff_val;
//	return rf_tag;
//}
