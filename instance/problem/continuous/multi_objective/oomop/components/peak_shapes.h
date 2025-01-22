
/*
* these classes construct peaks with different shapes 
*/

/*
* created by tanqingshan on June 19th, 2022.
*/

#ifndef OFEC_PEAK_SHAPES_H
#define OFEC_PEAK_SHAPES_H

#include<vector>
#include"../../../../../../utility/functional.h"
#include"../../../../../../utility/linear_algebra/matrix.h"
#include"../../../../../../utility/linear_algebra/vector.h"

namespace ofec {
	class Peak {
	protected:
		std::vector<Real> m_variables;
		Real m_H;
		std::vector<Real> m_dim_norm;
		Real m_exp = 2;
		Real m_rotation = 0;
	public:
		Peak(const std::vector<Real>& x, Real h, const std::vector<Real>& dim_norm) :m_variables(x), m_H(h), m_dim_norm(dim_norm) {}
		Peak(const std::vector<Real>& x, Real h, const std::vector<Real>& dim_norm, Real exp) :\
			m_variables(x), m_H(h), m_dim_norm(dim_norm), m_exp(exp) {}
		Peak(const std::vector<Real>& x, Real h, const std::vector<Real>& dim_norm, Real exp, Real rotation) :\
			m_variables(x), m_H(h), m_dim_norm(dim_norm), m_exp(exp), m_rotation(rotation) {}
		std::vector<Real> getPeakPos() { return m_variables; }
		Real getHeight() { return m_H; }
		std::vector<Real> getDimNorm() { return m_dim_norm; }
		Real getExp() { return m_exp; }
		Real getRotation() { return m_rotation; }

		void setPeakPos(const std::vector<Real>& x) { m_variables = x; }
		void setHeight(Real h) { m_H = h; }
		void setDimNorm(const std::vector<Real>& norm) { m_dim_norm = norm; }
		void setExp(Real exp) { m_exp = exp; }
		void setRotation(Real r) { m_rotation = r; }

		virtual Real calPeakValue(const std::vector<Real>& sol) { return 0; }
	};

	class Linear_peak :public Peak {
	private:
		Real m_slope;
	public:
		Linear_peak(const std::vector<Real>& x, Real h, Real s, const std::vector<Real>& dim_norm, Real exp) :\
			Peak(x, h, dim_norm, exp), m_slope(s) {}
		Linear_peak(const std::vector<Real>& x, Real h, Real s, const std::vector<Real>& dim_norm, Real exp, Real rotation) :\
			Peak(x, h, dim_norm, exp, rotation), m_slope(s) {}

		Real getSlope() { return m_slope; }
		void setSlope(Real s) { m_slope = s; }

		Real calPeakValue(const std::vector<Real>& sol) {
			if (m_rotation == 0) {
				Real dist = 0.;
				for (size_t i = 0; i < sol.size(); ++i) {
					dist = dist + std::pow(std::fabs(sol[i] - m_variables[i]), m_dim_norm[i]);
				}
				return m_H - m_slope * std::pow(dist, 1 / m_exp);
			}
			else {
				Vector temp_sol(sol);
				Vector temp_peak(m_variables);
				Matrix rotate_matrix(2, 2);
				rotate_matrix[0][0] = std::cos(-1 * m_rotation);
				rotate_matrix[0][1] = -1 * std::sin(-1 * m_rotation);
				rotate_matrix[1][0] = std::sin(-1 * m_rotation);
				rotate_matrix[1][1] = std::cos(-1 * m_rotation);
				temp_sol = temp_sol - temp_peak;
				Matrix rotate_point(temp_sol.size(), 1);
				for (size_t i = 0; i < temp_sol.size(); ++i) {
					rotate_point[i][0] = temp_sol[i];
				}
				if (sol.size() == 2) {
					temp_sol = (rotate_matrix * rotate_point).transpose().data()[0];
				}
				temp_sol = temp_sol + m_variables;
				Real dist = 0.;
				for (size_t i = 0; i < sol.size(); ++i) {
					dist = dist + std::pow(std::fabs(temp_sol[i] - m_variables[i]), m_dim_norm[i]);
				}
				return m_H - m_slope * std::pow(dist, 1 / m_exp);
			}
		}
	};

	class Nonliear_peak :public Peak {
	private:
		Real m_width;
	public:
		Nonliear_peak(const std::vector<Real>& x, Real h, Real w, const std::vector<Real>& dim_norm, Real exp) :\
			Peak(x, h, dim_norm, exp), m_width(w) {}
		Nonliear_peak(const std::vector<Real>& x, Real h, Real w, const std::vector<Real>& dim_norm, Real exp, Real rotation) :\
			Peak(x, h, dim_norm, exp, rotation), m_width(w) {}

		Real getWidth() { return m_width; }
		void setWidth(Real w) { m_width = w; }

		Real calPeakValue(const std::vector<Real>& sol) {
			if (m_rotation == 0) {
				Real dist = 0.;
				for (size_t i = 0; i < sol.size(); ++i) {
					dist = dist + std::pow(std::fabs(sol[i] - m_variables[i]), m_dim_norm[i]);
				}
				return m_H / (1 + m_width * std::pow(dist, 1 / m_exp));
			}
			else {
				Vector temp_sol(sol);
				Vector temp_peak(m_variables);
				Matrix rotate_matrix(2, 2);
				rotate_matrix[0][0] = std::cos(-1 * m_rotation);
				rotate_matrix[0][1] = -1 * std::sin(-1 * m_rotation);
				rotate_matrix[1][0] = std::sin(-1 * m_rotation);
				rotate_matrix[1][1] = std::cos(-1 * m_rotation);
				temp_sol = temp_sol - temp_peak;
				Matrix rotate_point(temp_sol.size(), 1);
				for (size_t i = 0; i < temp_sol.size(); ++i) {
					rotate_point[i][0] = temp_sol[i];
				}
				if (sol.size() == 2) {
					temp_sol = (rotate_matrix * rotate_point).transpose().data()[0];
				}
				temp_sol = temp_sol + m_variables;
				Real dist = 0.;
				for (size_t i = 0; i < sol.size(); ++i) {
					dist = dist + std::pow(std::fabs(temp_sol[i] - m_variables[i]), m_dim_norm[i]);
				}
				return m_H / (1 + m_width * std::pow(dist, 1 / m_exp));
			}
		}
	};
}

#endif // !OFEC_PEAK_SHAPES_H


