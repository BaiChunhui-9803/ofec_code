#ifndef SELECTIOIN_PROBLEM_TTFUNPAR_H
#define SELECTIOIN_PROBLEM_TTFUNPAR_H

#include"ttFun.h"


namespace ofec {
	namespace sp {
		enum class TypeWall {
			SMALL_LOOP = 0,
			LARGE_LOOP = 1,
			CIRCLE = 2,
			NEEDLE = 3
		};
		struct ParamMapInfo {
			int m_T = 1000;

			// 
			std::vector<int> m_wall_num = { 0,0,20,0 };
			//! parameters
			std::pair<double, double> m_width_range = { 1.0,2.0 };
			// wall character
			double m_max_radius = 2;
			double m_small_loop_radius = 0.1;
			double m_large_loop_inner = 0.1;
			int m_num_min_points = 5;
			int m_num_small_points = 10;
			int m_num_large_points = 20;
			int m_num_circle_points = 30;
			int m_num_needle_points = 20;


			int m_station_type_num = 10;
			std::pair<double, double> m_discount_level = { 0.6,0.9 };
			std::pair<double, double> m_discount_ratio = { 0.8,1.0 };
			Real m_moving_pos_ratio = 1.0;
			Real m_wallinfo_noisy_ratio = 0.1;


			const double m_movingt_speed_ratio = 0.3;
			const std::pair<double, double> m_angle_speed = { OFEC_PI * m_movingt_speed_ratio,OFEC_PI * 2.0 };
			const std::pair<double, double> m_radius_speed = { m_movingt_speed_ratio,2.0 };
			                                                                                                                                                                                                                    

			const double m_angle_noisy_offset = OFEC_PI * 0.01;
			const double m_radius_noisy_offset = 2.0 * 0.01;



		

		};


		struct FunsTtPar {
			std::vector<fun_gen::elementaryFunctionParameters> m_pars;
			FunTt::FunType m_fun_type = FunTt::ADD;
			std::array<bool, 2> m_fun_exist = { true,false };
			double m_obj_noisy_ratio = 0.1;
			Random::random_type m_random_type = Random::NORMAL_TYPE;
			FunsTt::DtriRangeType  m_range_type = FunsTt::ORIGIN;
			double m_feasible_threshold = 0.1;
			std::array<double, 2> m_scale_xy 
				= { 0.1,0.1 };
			std::array<double, 2> m_offset_xy 
				= { 1.0,0 };

			static double mc_obj_noisy_ratio;
			

			void setDefault() {
				m_fun_type = FunTt::ADD;
				m_fun_exist = { true,true };
				m_obj_noisy_ratio = mc_obj_noisy_ratio;
				m_random_type = Random::NORMAL_TYPE;
				m_range_type = FunsTt::ORIGIN;
				m_feasible_threshold = 0.1;
				m_scale_xy = { 0.1,0.1 };
			}
		};

		//struct TestPointPos {
		//	std::pair<double, double> real_radius_range = { 0.0,2.0 };
		//	std::pair<double, double> radius_range = { 0.8,2.0 };
		//	std::pair<double, double> real_angle_range = { 0,2.0 * OFEC_PI };
		//	std::pair<double, double> angle_range = { 0 * OFEC_PI ,2.0 * OFEC_PI };
		//	int sampleT = 5;
		//	int samplet = 5;

		//	double t_ratio = 0.1;
		//};


		//extern TestPointPos g_testWall;
		struct RandonFunParGen {

			//bool m_flag_change_execution = false;



			int m_edge_max_cos_fun = 4;
			const int mc_T = 1000;
			int m_curT = 1000;
			double sampletRatio = 0.2;
			double sampleEdgeScaleRatio = 0.05;
			double sampleTRatio = 0.3;
			double samplePRatio = 0.1;

			int m_samplePTimes = 8;
			const std::vector<std::pair<double, double>> m_x_from = 
			{ {0,0.5},{0,1},{0,2} };
			std::vector<std::pair<double, double>> m_fun_ranges=
			{ {0,0.5},{0,1},{0,1} };
			std::vector<int> m_fun_num;
			std::vector<int> m_sample_num;
			std::vector<double>  m_numPI;
			std::vector<fun_gen::elementary_fun_type> m_fun_flag;
			//double m_noisy_ratio = 0.01;
			Random::random_type m_random_type = Random::NORMAL_TYPE;
		//	FunTt::FunType m_fun_type;

			std::pair<double, double> m_radius_range = { 0,2.0 };
			std::pair<double, double> m_angle_range = { 0 * OFEC_PI ,2.0 * OFEC_PI };

			


			int m_sampleT = 10;
			int m_samplet = 10;
			double m_range_scale = 4.0;
			double m_point_t_moving_range = 0.5;

			double m_env_noisy_ratio = 0.1;
		

			std::pair<double, double> m_num_dynamic_centers_ratio
			= {0.01,0.03};
			
			FunsTtPar m_par;		

			size_t m_num_basic_sols = 1;
			size_t m_num_feasible_basic_sols_states = 2;
			double m_change_env_factors_raito = 0.0;
			double m_init_env_factors_dim_raito = 0.0;
			double m_change_env_factors_dim_raito = 0.0;

			void initialize(int cur_T) {
				m_curT = cur_T;
				m_fun_ranges.resize(3);
	/*			m_radius_noisy_range = m_radius_range;
				m_radius_noisy_range.second *= (1.0 + m_noisy_ratio* FunsTt::Three_Sigma_Limit);
				m_angle_noisy_range = m_angle_range;
				m_angle_noisy_range.second *= (1.0 + m_noisy_ratio * FunsTt::Three_Sigma_Limit);
	*/		}

			void initDynamicIdxs(Random *rnd,
				int num_centers,std::vector<int>& center_idxs) {
				int center_size = std::max<double>(
					1,
					std::round(m_curT * rnd->uniform.nextNonStd<double>
					(m_num_dynamic_centers_ratio.first, m_num_dynamic_centers_ratio.second)));
				center_idxs.resize(center_size);
				for (auto& it : center_idxs) {
					it = rnd->uniform.nextNonStd<int>(0, num_centers);
				}
			}



			void initDynamicVarStationTest(
				const std::pair<double, double>& range1,
				const std::pair<double, double>& range2,
				int samplet, int sampleT
			) {
				//m_fun_type = FunTt::FunType::Add;
				//m_x_from = { {0,1},{0,1},{0,2} };
				m_fun_ranges = { {0,0.5} ,range1,range2 };
				m_fun_num = { 0,0,0 };
				m_sample_num =
				{ sampleT ,sampleT,samplet };
				m_numPI = { 0,0,3 };
				m_fun_flag = {
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::random_fun };
			}


			//void initFunTestAngle(Random *rnd) {
			//	std::pair<double, double> range1, range2;
			//	range1.first = 0;
			//	range1.second = m_test_para.angle_range.second - m_test_para.angle_range.first;
			//	range1.second *= 2;
			//	range2 = range1;
			//	range2.second *= m_test_para.t_ratio;
			//	initDynamicVarStationTest(range1, range2,
			//		m_test_para.samplet, m_test_para.sampleT);
			////	range1.first = range1.second = 1.0;
			////	range2.first = range2.second = 1.0;
			//	m_par.setDefault();
			//	m_par.m_fun_type = FunTt::FunType::ADD;
			//	initPar(rnd);
			//}


			//void initFunTestRadius(Random *rnd) {
			//	std::pair<double, double> range1, range2;
			//	range1.first = 0;
			//	range1.second = m_test_para.radius_range.second - m_test_para.radius_range.first;
			//	range1.second *= 2;
			//	range2 = range1;
			//	range2.second *= m_test_para.t_ratio;
			//	//range1.first = range1.second = 1.0;
			//	//range2.first = range2.second = 1.0;
			//	initDynamicVarStationTest(range1, range2,
			//		m_test_para.samplet, m_test_para.sampleT);
			//	m_par.setDefault();
			//	m_par.m_fun_type = FunTt::FunType::ADD;
			//	initPar(rnd);
			//}


			void initDynamicVarStation(
				Random *rnd,
				const std::array<double, 3>& cur_range,
				bool flag_radius = true) {
				std::pair<double, double> range;
				if (flag_radius) {
					range = m_radius_range;
				}
				else {
					range = m_angle_range;
				}
				double cur_gap = cur_range[1] - cur_range[0];
				double gap = range.second - range.first;
				double scale = cur_gap / gap;

				int sample_T = std::max<int>(2,std::round(m_sampleT *
					double(m_curT) / double(mc_T) * scale));
				int sample_t = std::max<int>(2,std::round(m_samplet * scale));
				m_fun_ranges[0] = { 0,0.5 };
				m_fun_ranges[1].first = 0;
				m_fun_ranges[1].second = cur_gap;
				m_fun_ranges[1].second *= m_range_scale;
				m_fun_ranges[2] = m_fun_ranges[1];
				m_fun_ranges[2].second *= m_point_t_moving_range;

				m_fun_num = { 0,0,0 };
				m_sample_num =
				{ sample_T ,sample_T,sample_t };
				m_numPI = { 0,0,3 };
				m_fun_flag = {
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::random_fun };

				m_par.setDefault();
				m_par.m_fun_type = FunTt::FunType::ADD;
				initPar(rnd);
			}
			void initFeasiblePar(
				Random *rnd
			) {
				int sampelt = m_x_from.back().second / sampletRatio;
				int sampelT = sampleTRatio * m_curT;
				m_fun_ranges = { {0,0.5} ,{0.8,1.0},
					{0.8,1.0} };
				m_fun_num = { 0,0,m_edge_max_cos_fun };
				m_sample_num =
				{ sampelT ,sampelT,sampelt };
				m_numPI = { 0,0,3 };
				m_fun_flag = {
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::cos_fun };


				m_par.setDefault();
				m_par.m_fun_type = FunTt::FunType::TIMES;
				initPar(rnd);
			}


			void initPriceVarPar() {
				int sampelt = m_x_from.back().second / sampletRatio;
				int sampelT = sampleTRatio * m_curT;
				m_fun_ranges = { {0,0.5} ,{1.0,1.5},{1.0,1.5} };
				//	m_fun_ranges = { {0,0} ,{1.0,2.0}, {1.0,1.0} };
				m_fun_num = { 0,0,0 };
				m_sample_num =
				{ sampelT ,sampelT,sampelt };
				m_numPI = { 0,0,0 };
				m_fun_flag = {
			fun_gen::elementary_fun_type::constant_fun,
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::constant_fun };


				m_par.setDefault();
				m_par.m_fun_type = FunTt::FunType::TIMES;
			}

			void initPriceVarPar(Random *rnd) {
				initPriceVarPar();
				initPar(rnd);
			}

			
			std::pair<double, double> getFunTtRange(bool noisy_offset=false) {
				auto range = FunTt::getRange(m_fun_ranges, m_par.m_fun_type);
				if (noisy_offset) {
					range.second *= (1.0 + 3.0 * FunsTt::Three_Sigma_Limit);
				}
				return FunTt::getRange(m_fun_ranges, m_par.m_fun_type);
			}

			void initRunningTimeVarPar() {
				int sampelt = m_x_from.back().second / sampleEdgeScaleRatio;
				int sampelT = sampleTRatio * m_curT;
				m_fun_ranges = { {0,0.5} ,{1.0,1.5},{1.0,1.5} };
				m_fun_num = { 0,0,m_edge_max_cos_fun };
				m_sample_num =
				{ sampelT ,sampelT,sampelt };
				m_numPI = { 0,0,3 };
				m_fun_flag = {
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::cos_fun };

				m_par.setDefault();
				m_par.m_fun_type = FunTt::FunType::TIMES;
			}


			void initRunningTimeVarPar(Random *rnd) {
				initRunningTimeVarPar();
				initPar(rnd);
			}


			void initMeshSeedFun(Random *rnd) {
			//	m_fun_type = FunTt::FunType::Add;
			//	m_x_from = { {0,mc_T},{0,mc_T},{0,2} };
				int sampelt = m_samplePTimes;
				int sampelT = sampleTRatio * mc_T;

				//	double m_noisy_offset = 0.1;
				m_fun_ranges = { {0,0.5} ,{0,100},
					{0, 10} };
				m_fun_num = { 0,0,m_edge_max_cos_fun };
				m_sample_num =
				{ sampelT ,sampelT,sampelt };
				m_numPI = { 0,0,3 };
				m_fun_flag = {
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::random_fun,
			fun_gen::elementary_fun_type::cos_fun };

				m_par.setDefault();
				m_par.m_fun_type = FunTt::FunType::ADD;
				initPar(rnd);
			}



			void initPar(Random *rnd = -1) {
				m_par.m_pars.resize(FunTt::m_numFun);
				auto fun_ranges = m_fun_ranges;
				auto fun_num = m_fun_num;
				if (rnd != -1) {
					for (int idx(0); idx < FunTt::m_numFun; ++idx) {
						if (m_fun_ranges[idx].first != m_fun_ranges[idx].second)
							fun_ranges[idx].second = rnd->uniform.nextNonStd<double>(m_fun_ranges[idx].first, m_fun_ranges[idx].second);
					}

					for (int idx(0); idx < FunTt::m_numFun; ++idx) {
						if (fun_num[idx] != 0) {
							fun_num[idx] = rnd->uniform.nextNonStd<int>(1, m_fun_num[idx]);
						}
					}
				}

				for (int idx(0); idx < FunTt::m_numFun; ++idx) {
					m_par.m_pars[idx].set(m_fun_flag[idx], m_x_from[idx],
						fun_ranges[idx], m_sample_num[idx],
						fun_num[idx], m_numPI[idx]);
				}
			}



		};

		struct ParaInfo {
			int mc_T = 100;
			ParamMapInfo m_static_info;
			RandonFunParGen m_dynamic_info;
			//ParFunTt m_fun_para;
			//std::vector<fun_gen::elementaryFunctionParameters> m_pars;
			//TestPointPos m_test_para;
			void initialize() {
				m_static_info.m_T = mc_T;
				m_dynamic_info.initialize(mc_T);
			}
		};


		//struct ParFunTt {
		//	int m_edge_max_cos_fun = 4;
		//	int mc_T = 1000;
		//	std::vector<std::pair<double, double>> m_x_from = { {0,mc_T},{0,mc_T},{0,1} };
		//	std::vector<std::pair<double, double>> m_fun_ranges;
		//	std::vector<int> m_fun_num;
		//	std::vector<int> m_sample_num;
		//	std::vector<double>  m_numPI;
		//	std::vector<fun_gen::elementary_fun_type> m_fun_flag;
		//	double sampletRatio = 0.2;
		//	double sampleTRatio = 0.3;
		//	double samplePRatio = 0.5;
		//	//double m_noisy_offset = 0.3;
		//	FunTt::FunType m_fun_type = FunTt::FunType::Add;

		//	void initialize(int T) {
		//		mc_T = T;
		//	}




		//	void initEdgeNoisePar() {
		//		m_fun_type = FunTt::FunType::Times;
		//		m_x_from = { {0,mc_T},{0,mc_T},{0,2} };
		//		//	double m_noisy_offset = 0.01;
		//		int sampelt = m_x_from.back().second / sampletRatio;
		//		int sampelT = sampleTRatio * mc_T;
		//		m_fun_ranges = { {0,0.5} ,{0,1.0},
		//			{0,m_noisy_offset} };
		//		m_fun_num = { 0,0,0 };
		//		m_sample_num =
		//		{ sampelT ,sampelT,sampelt };
		//		m_numPI = { 0,0,0 };
		//		m_fun_flag = {
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun };
		//	}

		//	void initEdgeVarPar() {
		//		m_x_from = { {0,mc_T},{0,mc_T},{0,2} };
		//		m_fun_type = FunTt::FunType::Times;
		//		int sampelt = m_x_from.back().second / sampletRatio;
		//		int sampelT = sampleTRatio * mc_T;
		//		m_fun_ranges = { {0,0.5} ,{},
		//			{0,1} };
		//		m_fun_num = { 0,0,m_edge_max_cos_fun };
		//		m_sample_num =
		//		{ sampelT ,sampelT,sampelt };
		//		m_numPI = { 0,0,3 };
		//		m_fun_flag = {
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::cos_fun };
		//	}

		//	std::array<double, 2> getRange(double noisy_offset = 0.0) {
		//		std::array<double, 2> range;
		//		if (m_fun_type == FunTt::FunType::Times) {
		//			range[0] = m_fun_ranges[1].first * m_fun_ranges[2].first;
		//			range[1] = m_fun_ranges[1].second * m_fun_ranges[2].second;
		//		}
		//		else if (m_fun_type == FunTt::FunType::Add) {
		//			range[0] = m_fun_ranges[1].first + m_fun_ranges[2].first;
		//			range[1] = m_fun_ranges[1].second + m_fun_ranges[2].second;
		//		}
		//		range[1] *= (1.0 + 3.0 * m_noisy_offset);

		//		return std::move(range);
		//	}
		//	std::array<double, 2> getEdgeVarRange(double noisy_offset = 0.0) {
		//		initEdgeVarPar();
		//		return std::move(getRange(noisy_offset));
		//	}

		//	std::array<double, 2> getPriceVarRange(double noisy_offset = 0.0) {
		//		initPriceVarPar();
		//		return std::move(getRange(noisy_offset));
		//	}






		//	void initDynamicVarStationTest(
		//		const std::pair<double, double>& range1,
		//		const std::pair<double, double>& range2,
		//		int samplet, int sampleT
		//	) {
		//		m_fun_type = FunTt::FunType::Add;
		//		m_x_from = { {0,1},{0,1},{0,2} };
		//		m_fun_ranges = { {0,0.5} ,range1,range2 };
		//		m_fun_num = { 0,0,0 };
		//		m_sample_num =
		//		{ sampleT ,sampleT,samplet };
		//		m_numPI = { 0,0,3 };
		//		m_fun_flag = {
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun };
		//	}



		//	void initDynamicVarStationTest(
		//		const std::pair<double, double>& range) {
		//		m_fun_type = FunTt::FunType::Add;
		//		m_x_from = { {0,mc_T},{0,mc_T},{0,2} };
		//		int sampelt = 5;
		//		int sampelT = 300;
		//		m_fun_ranges = { {0,0.5} ,range,range };
		//		m_fun_num = { 0,0,0 };
		//		m_sample_num =
		//		{ sampelT ,sampelT,sampelt };
		//		m_numPI = { 0,0,3 };
		//		m_fun_flag = {
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun };
		//	}





		//	void initDynamicVarStation(
		//		const std::pair<double, double>& funt,
		//		const std::pair<double, double>& funT,
		//		const std::pair<double, double>& speed,
		//		FunTt::FunType fun_type = FunTt::FunType::Add) {
		//		m_fun_type = fun_type;
		//		m_x_from = { {0,mc_T},{0,mc_T},{0,2} };
		//		int sampelt = std::round(m_x_from.back().second * speed.first);
		//		int sampelT = std::round(m_x_from.front().second * speed.second);
		//		m_fun_ranges = { {0,0.5} ,funt,funT };
		//		m_fun_num = { 0,0,0 };
		//		m_sample_num =
		//		{ sampelT ,sampelT,sampelt };
		//		m_numPI = { 0,0,3 };
		//		m_fun_flag = {
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun };
		//	}

		//	void initStaticVarStation(
		//		const std::pair<double, double>& funt,
		//		const std::pair<double, double>& funT) {
		//		m_fun_type = FunTt::FunType::Add;
		//		m_x_from = { {0,mc_T},{0,mc_T},{0,2} };
		//		int sampelt = 1;
		//		int sampelT = 1;
		//		m_fun_ranges = { {0,0.5} ,{funt.first,funt.second},
		//			{funT.first,funT.second} };
		//		m_fun_num = { 0,0,0 };
		//		m_sample_num =
		//		{ sampelT ,sampelT,sampelt };
		//		m_numPI = { 0,0,3 };
		//		m_fun_flag = {
		//	fun_gen::elementary_fun_type::constant_fun,
		//	fun_gen::elementary_fun_type::constant_fun,
		//	fun_gen::elementary_fun_type::constant_fun };
		//	}


		//	//void initDynamicNoisyStation(double noisy_offset) {
		//	//	m_fun_type = FunTt::FunType::Add;
		//	//	m_x_from = { {0,mc_T},{0,mc_T},{0,2} };
		//	//	int sampelt = m_x_from.back().second / sampletRatio;
		//	//	int sampelT = sampleTRatio * mc_T;

		//	////	double m_noisy_offset = 0.1;
		//	//	m_fun_ranges = { {0,0.5} ,{0,noisy_offset},
		//	//		{0,noisy_offset} };
		//	//	m_fun_num = { 0,0,m_edge_max_cos_fun };
		//	//	m_sample_num =
		//	//	{ sampelT ,sampelT,sampelt };
		//	//	m_numPI = { 0,0,3 };
		//	//	m_fun_flag = {
		//	//fun_gen::elementary_fun_type::random_fun,
		//	//fun_gen::elementary_fun_type::random_fun,
		//	//fun_gen::elementary_fun_type::random_fun };
		//	//}




		//	void initMeshSeedFun() {
		//		m_fun_type = FunTt::FunType::Add;
		//		m_x_from = { {0,mc_T},{0,mc_T},{0,2} };
		//		int sampelt = m_x_from.back().second / samplePRatio;
		//		int sampelT = sampleTRatio * mc_T;

		//		//	double m_noisy_offset = 0.1;
		//		m_fun_ranges = { {0,0.5} ,{0,20},
		//			{0, 5} };
		//		m_fun_num = { 0,0,m_edge_max_cos_fun };
		//		m_sample_num =
		//		{ sampelT ,sampelT,sampelt };
		//		m_numPI = { 0,0,3 };
		//		m_fun_flag = {
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::random_fun,
		//	fun_gen::elementary_fun_type::cos_fun };
		//	}




		//	void initPar(
		//		std::vector<fun_gen::elementaryFunctionParameters>& pars,
		//		Random *rnd = -1) {
		//		pars.resize(FunTt::m_numFun);
		//		auto fun_ranges = m_fun_ranges;
		//		auto fun_num = m_fun_num;
		//		if (rnd != -1) {
		//			for (int idx(0); idx < FunTt::m_numFun; ++idx) {
		//				if (m_fun_ranges[idx].first != m_fun_ranges[idx].second)
		//					fun_ranges[idx].second = rnd->uniform.nextNonStd<double>(m_fun_ranges[idx].first, m_fun_ranges[idx].second);
		//			}

		//			for (int idx(0); idx < FunTt::m_numFun; ++idx) {
		//				if (fun_num[idx] != 0) {
		//					fun_num[idx] = rnd->uniform.nextNonStd<int>(1, m_fun_num[idx]);
		//				}
		//			}
		//		}

		//		for (int idx(0); idx < FunTt::m_numFun; ++idx) {
		//			pars[idx].set(m_fun_flag[idx], m_x_from[idx],
		//				fun_ranges[idx], m_sample_num[idx],
		//				fun_num[idx], m_numPI[idx]);
		//		}
		//	}



		//	void initFun(
		//		const std::vector<fun_gen::elementaryFunctionParameters>& pars,
		//		FunTt& funs, 
		//		Random *rnd) {
		//		//	funs.m_funs.resize(FunTt::m_numFun);
		//		for (int idx(0); idx < FunTt::m_numFun; ++idx) {
		//			pars[idx].resetFun(funs.m_funs[idx]);
		//		}
		//		for (auto& it : funs.m_funs) {
		//			it->initialize(rnd->uniform);
		//		}
		//		funs.m_type = m_fun_type;
		//	}
		//};
	}
}

#endif 

