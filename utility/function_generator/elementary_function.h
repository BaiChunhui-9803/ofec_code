#ifndef OFEC_ELEMENT_FUNCTION_H
#define OFEC_ELEMENT_FUNCTION_H

#include <cmath>
#include <vector>
#include <string>
#include <memory>
#include"../cubic_spline/cubic_spline.h"
#include "../random/newran.h"
#include "../../core/definition.h"

namespace ofec {

	namespace fun_gen {

		enum class elementary_fun_type { constant_fun, random_fun, cos_fun };

		// map to [0,1]
		class  elementary_function {
		protected:
			// map the x to [0,1] 
			std::pair<double, double> m_from_x_range = { 0,1 };
			double m_fromXLen = 1.0;
			std::pair<double, double> m_to_y_range = { 0,1 };
			double m_toYLen = 1.0;
		protected:
			inline double toX(double x) {
				return (x - m_from_x_range.first) / m_fromXLen;
			}
		public:
			elementary_function(
				const std::pair<double, double>& from_x_range,
				const std::pair<double, double>& to_y_range)
				:m_from_x_range(from_x_range), m_to_y_range(to_y_range) {
				m_fromXLen = m_from_x_range.second - m_from_x_range.first;
				m_toYLen = m_to_y_range.second - m_to_y_range.first;

			};
			virtual double get_value(double x) = 0;
			virtual void set_from_x_range(double x_from, double x_to) {
				m_from_x_range.first = x_from;
				m_from_x_range.second = x_to;
				m_fromXLen = m_from_x_range.second - m_from_x_range.first;
			}
			//virtual void set_to_x_range(double x_from, double x_to) {
			//	m_to_x_range.first = x_from;
			//	m_to_x_range.second = x_to;
			//	m_toXLen = m_to_x_range.second - m_to_x_range.first;
			//}

			virtual void set_y_range(double y_from, double y_to) {
				m_to_y_range.first = y_from;
				m_to_y_range.second = y_to;
				m_toYLen = m_to_y_range.second - m_to_y_range.first;
			}
			const std::pair<double, double>& get_to_y_range() {
				return m_to_y_range;
			}
			virtual void initialize(RandBase & randomBase) {};
		};

		class constant_fun :public elementary_function {
		protected:
			double m_mean_value;
		public:
	//		constant_fun(const ParamMap& v) :elementary_function(v) {}
			template<typename ... Args>
			constant_fun(Args&& ... args) : elementary_function(std::forward<Args>(args)...) {
				m_mean_value = (m_to_y_range.first + m_to_y_range.second) / 2.0;

			};
			virtual void set_y_range(double y_from, double y_to) {
				elementary_function::set_y_range(y_from, y_to);
				m_mean_value = (m_to_y_range.first + m_to_y_range.second) / 2.0;

		  	}

			virtual double get_value(double x)override {
				return m_mean_value;
			}
			virtual void initialize(RandBase& randomBase)override {
			};
		};


		class random_fun :public elementary_function {
			CubicSpline<double> m_fun;
			int m_sample_num;
		public:
			template<typename ... Args>
			random_fun(int sample_num, Args&& ... args) : elementary_function(std::forward<Args>(args)...), m_sample_num(sample_num) {};
			
			virtual double get_value(double x)override {
				double val(m_fun.Interpolate_static_inte(x));
				if (val < m_to_y_range.first) val = m_to_y_range.first;
				else if (val > m_to_y_range.second)val = m_to_y_range.second;
				return val;
			}
			virtual void initialize(RandBase& randomBase)override {
				std::vector<double> x(m_sample_num);
				std::vector<double> y(m_sample_num);
				double x_offset(static_cast<double>(m_from_x_range.second - m_from_x_range.first) / static_cast<double>(m_sample_num));
				double from_x(m_from_x_range.first);
				double len(m_to_y_range.second - m_to_y_range.first);
				for (int sample_id(0); sample_id < m_sample_num; ++sample_id) {
					x[sample_id] = from_x + x_offset * sample_id;
					y[sample_id] = randomBase.next() * len + m_to_y_range.first;
				}
				m_fun.Initialize(x, y);
				m_fun.set_interval(x_offset, from_x);
			};
		};


		class cos_fun :public elementary_function {

			// fourier series  f(x)=\sum Ansin(n*omega*x+ phi )
		protected:
			std::vector<double> m_an;
			std::vector<double> m_pi;
			std::pair<double, double> m_real_y;
			double m_realYLen = 0;
	
			int m_fun_num = 1;
			double m_omega = 2 * OFEC_PI;
			int m_sample_times = 100;

		protected:
			//\sum An*sin(n*omega*x+ phi )
			inline double ox2y(double x) {
				double y(0);
				for (int idx(0); idx < m_fun_num; ++idx) {
					y += m_an[idx] * sin(m_omega * (idx + 1) * x + m_pi[idx]);
				}
				return y;
			}
			inline double oy2y(double y) {
				return (y - m_real_y.first) / m_realYLen * m_toYLen + m_to_y_range.first;
			}

		public:
			template<typename ... Args>
			cos_fun(int fun_num,double numPI,Args&& ... args) : 
				elementary_function(std::forward<Args>(args)...),
				m_fun_num(fun_num), m_omega(OFEC_PI*numPI), m_an(fun_num,0),m_pi(fun_num,0){};

			virtual double get_value(double x) override {
				double y(ox2y(toX(x)));
				if (y > m_real_y.second) y = m_real_y.second;
				else if (y < m_real_y.first) y = m_real_y.first;
				return oy2y(y);
			}

			virtual void initialize(RandBase& randomBase)override {
				for (auto& it : m_pi) {
					it = randomBase.next()*2 * OFEC_PI;
				}
				m_an.front() = 1.0;
				for (size_t idx(1); idx < m_an.size(); ++idx) {
					m_an[idx] = m_an[idx - 1] * randomBase.next() * 0.2 + 0.4;
				}
				{
					double maX = 1.0 * 2* OFEC_PI /m_omega;
					double x_offset(maX / double(m_sample_times));
					m_real_y = std::make_pair(1e9, -1e9);
					double cur_x(0), cur_y(0);
					while (cur_x < maX) {
						cur_y = ox2y(cur_x);
						m_real_y.first = std::min(m_real_y.first, cur_y);
						m_real_y.second = std::max(m_real_y.second, cur_y);
						cur_x += x_offset;
					}
					m_realYLen = m_real_y.second - m_real_y.first;
				}
			}
		};
		struct elementaryFunctionParameters {
			elementary_fun_type m_type;
			std::pair<double, double> m_from_x_range = { 0,1 };
			std::pair<double, double> m_to_y_range = { 0,1 };
			int m_sample_num;
			int m_fun_num;
			double m_numPI;
			void set(const elementary_fun_type& type,
				const std::pair<double, double>& rangex,
				const std::pair<double, double>& rangey,
				int sample_num,
				int fun_num,
				double numPI) {
				m_type = type;
				m_from_x_range = rangex;
				m_to_y_range = rangey;
				m_sample_num = sample_num;
				m_fun_num = fun_num;
				m_numPI = numPI;
			}
			void resetFun(std::unique_ptr<elementary_function>& fun)const ;
		};

		extern void elementary_function_initialize
		(std::unique_ptr<elementary_function>& fun, elementaryFunctionParameters par, RandBase& randomBase);



		//extern void elementary_function_initialize(std::unique_ptr<elementary_function>& fun,
		//	const elementaryFunctionParameters& par
		//);

		//extern void elementary_function_initialize(std::unique_ptr<elementary_function>& fun,
		//	elementary_function::elementary_fun_type type,
		//	const std::pair<double, double>& from_x, const std::pair<double, double>& to_x,
		//	const std::pair<double, double>& to_y,
		//	Random *rnd, int fun_num, int sample_num);


		//extern bool elementary_function_generator(std::unique_ptr<elementary_function>& fun, elementary_function::elementary_fun_type type,
		//	const std::pair<double, double>& from_x, const std::pair<double, double>& to_x,
		//	const std::pair<double, double>& to_y,
		//	int* rnd = nullptr, int* fun_num = nullptr, int* sample_num = nullptr
		//);

	}


	


}

#endif