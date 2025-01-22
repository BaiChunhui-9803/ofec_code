#ifndef SELECTIOIN_PROBLEM_TTFUN_H
#define SELECTIOIN_PROBLEM_TTFUN_H

#include "../../../../utility/function_generator/elementary_function.h"
#include "../../../../utility/function_generator/random_smooth_function/noisy_fun.h"
#include "../../../../utility/random/newran.h"
#include<array>
#include<memory>
#include<vector>

namespace ofec {
	namespace sp {
		extern struct FunsTtPar;
		// limited: T \in  [0, 1] , t\in [0,1]
		// times: val =f2(t+f0(T))*f1(T), add:  val =f2(t+f0(T))+f1(T) 
		struct FunTt {
			static const int m_numFun;
			enum FunType { ADD, TIMES };
			FunType m_type = ADD;
			std::array<std::unique_ptr<fun_gen::elementary_function>,3> 
				m_funs = {nullptr,nullptr,nullptr};
		
			double getVal(int idx, double T, double t)const {
				if (idx < 2) {
					return m_funs[idx]->get_value(T);
				}
				else return m_funs[idx]->get_value(t);
			}
			double getVal(double T, double t) const{
				switch (m_type)
				{
				case ADD: {
					return m_funs[2]->get_value(t +
						m_funs[0]->get_value(T)) + m_funs[1]->get_value(T);
					break;
					}
				case TIMES: {
					return m_funs[2]->get_value(t +
						m_funs[0]->get_value(T)) * m_funs[1]->get_value(T);
					break;
				}
				default:
					return 0;
					break;
				}
			}
			//
			static std::pair<double, double> getRange(const std::vector<std::pair<double, double>>& origin_range, FunType type) {
				std::pair<double, double> range;
				switch (type)
				{
				case ADD: {
					range.first = origin_range[1].first + origin_range[2].first;
					range.second = origin_range[1].second+ origin_range[2].second;
					break;
				}
				case TIMES: {
					range.first = origin_range[1].first * origin_range[2].first;
					range.second = origin_range[1].second * origin_range[2].second;

					break;
				}
				default:
					range.first = 0;
					range.second = 0;
					break;
				}
				return std::move(range);
			}
			
			void getRange(std::pair<double, double>& range)const {
				switch (m_type)
				{
				case ADD: {
					range.first = m_funs[1]->get_to_y_range().first + m_funs[2]->get_to_y_range().first;
					range.second = m_funs[1]->get_to_y_range().second + m_funs[2]->get_to_y_range().second;
					break;
				}
				case TIMES: {
					range.first = m_funs[1]->get_to_y_range().first * m_funs[2]->get_to_y_range().first;
					range.second = m_funs[1]->get_to_y_range().second * m_funs[2]->get_to_y_range().second;
					break;
				}
				default:
					range.first = 0;
					range.second = 0;
					break;
				}
			}
			double getRange()const {
				std::pair<double, double> range;
				getRange(range);
				return range.second;
				//switch (m_type)
				//{
				//case ADD: {
				//	return m_funs[1]->get_to_y_range().second + m_funs[2]->get_to_y_range().second;
				//	break;
				//}
				//case TIMES: {
				//	return m_funs[1]->get_to_y_range().second * m_funs[2]->get_to_y_range().second;
				//	break;
				//}
				//default:
				//	return 0;
				//	break;
				//}
			}
		};

		//struct FunTt

		class FunsTt {
		public:
			enum DtriRangeType {
				ORIGIN = 0,
				WITHIN3SIGMA = 1,
				ZERO = 2
			};
			// 3σ准则又称为拉依达准则
			static const double Three_Sigma_Limit;
			static double getRandomWithin3Sigma(double rand);
		protected:


		protected:
			std::array<std::unique_ptr<sp::FunTt>, 2> m_funs;
			double m_noisy_ratio = 0.01;
			Random::random_type m_random_type = Random::NORMAL_TYPE;
			DtriRangeType  m_range_type = DtriRangeType::ORIGIN;
			double getRandVal(Random *rnd)const {
				switch (m_range_type)
				{
				case ORIGIN: {
					return rnd.next(m_random_type);
					break;
				}
				case WITHIN3SIGMA: {
					return getRandomWithin3Sigma(rnd.next(m_random_type));
					break;
				}
				default:
					return 0;
					break;
				}
			}
		public:
			void initialize(Random *rnd, const FunsTtPar& par);
		};
		class DistributionFunTt: public FunsTt {
		protected:
			double getVal(double T, double t, double rand_num) const {
				double mu(0), sigma(0);
				getDistribution(T, t, mu, sigma);
				return mu + sigma * rand_num * m_noisy_ratio;
			}
		public:
			void getDistribution(double T, double t, double& mu, double& sigma)const {
				mu = sigma = 0;
				if (m_funs[0] != nullptr) {
					mu = m_funs[0]->getVal(T, t);
				}
				if (m_funs[1] != nullptr) {
					sigma = m_funs[1]->getVal(T, t);
				}
			}
			double getVal(double T, double t, Random *rnd)const {
				return getVal(T, t, getRandVal(rand));
			}
			double getVal(double T, double t, Random *rnd = -1)const {
				if (rnd == -1) {
					double val(0), sigma(0);
					getDistribution(T, t, val, sigma);
					return val;
				}
				else return getVal(T, t, rnd);
			}

			double getMaxVal()const{
				double mu(0), sigma(0);
				if (m_funs[0] != nullptr) {
					mu = m_funs[0]->getRange();
				}
				if (m_funs[1] != nullptr) {
					sigma = m_funs[1]->getRange();
				}
				return mu + sigma * Three_Sigma_Limit * m_noisy_ratio;
			}
		};

		// f_val< m_feasible_threshold , the edge is unfeasible
		class FeasibleDistributionFunTt : public DistributionFunTt {
		protected:
			double m_feasible_threshold = 0.1;
			bool isFeasible(double val)const {
				return val >= m_feasible_threshold;
			}
		public:
	
			bool isFeasible(double T, double t, Random *rnd) const {
				return isFeasible(getVal(T, t, rnd));
			}
			bool isFeasible(double T, double t, Random *rnd=-1) const {
				return isFeasible(getVal(T, t, rnd));
			}
			void initialize(Random *rnd, const FunsTtPar& par);
		};
		class Distribution3DFunTt : public FunsTt {
		protected:
			std::array<double, 2> m_scale_xy = { 0.1,0.1 };
			std::array<double, 2> m_offset_xy = { 1.0,0.0};

			double getVal(double x, double y, double T, double t, double rand_num) const {
				double mu(0), sigma(0);
				getDistribution(x, y, T, t, mu, sigma);
				return mu + sigma * rand_num * m_noisy_ratio;
			}
			void scaleXY(double& x, double& y) const {
				x *= m_scale_xy[0];
				y *= m_scale_xy[1];
			}
		public:
			void getDistribution( double x,double y, double T, double t, double& mu, double& sigma)const {
				scaleXY(x, y);
				mu = m_offset_xy[0];
				sigma = m_offset_xy[1];
				if (m_funs[0] != nullptr) {
					mu+= noisy_fun::m_ofNoise(x, y, m_funs[0]->getVal(T, t));
				}
				if (m_funs[1] != nullptr) {
					sigma += noisy_fun::m_ofNoise(x, y, m_funs[1]->getVal(T, t));
				}
			}
			double getVal(double x, double y, double T, double t, Random *rnd )const {
				return getVal(x, y, T, t, getRandVal(rand));
			}
			double getVal(double x, double y, double T, double t, Random *rnd=-1)const {
				if (rnd == -1) {
					double val(0), sigma(0);
					getDistribution(x, y, T, t,  val, sigma);
					return val;
				}
				else {
					return getVal(x, y, T, t, getRandVal(rnd));
				}
			}
			void initialize(Random *rnd, const FunsTtPar& par);
		};
	}
}
#endif