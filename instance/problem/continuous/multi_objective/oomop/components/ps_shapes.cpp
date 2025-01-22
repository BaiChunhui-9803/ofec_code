#include"ps_shapes.h"

namespace ofec {
	Real PS_shape::calOneDimValue(const Real& offset, size_t idx) {
		auto value = calValue(offset);
		return value[idx];
	}

	bool PS_shape::if_in_range(std::vector<Real>& sol) {
		std::vector<Real> temp_sol;
		for (size_t i = 0; i < sol.size(); ++i) {
			if (i == 0) {
				if ((sol[i] - (std::min(m_points[0][0], m_points[1][0]) - 10e-6)) * (sol[i] - (std::max(m_points[1][0], m_points[0][0]) + 10e-6)) > 0) {
					return false;
				}
				else {
					temp_sol = calValue(sol[0] - m_points[0][0]);
				}
			}
			else {
				if (temp_sol[i] < temp_sol[i] + m_offset[i]) {
					if ((sol[i] - (temp_sol[i] - 10e-6)) * (sol[i] - temp_sol[i] - m_offset[i] - 10e-6) > 0) {
						return false;
					}
				}
				else {
					if ((sol[i] - (temp_sol[i] + 10e-6)) * (sol[i] - (temp_sol[i] + m_offset[i] - 10e-6)) > 0) {
						return false;
					}
				}
			}
		}
		return true;
	}

	//attached points and distance
	std::pair<std::vector<Real>, Real> PS_shape::dist_from_PS(const std::vector<Real>& sol, const std::vector<Real>& p_norms) {
		std::vector<Real> dim_dist;
		std::vector<Real> attach_point;
		bool in_range = false;//第一维是否在范围内
		if ((sol[0] - m_points[0][0]) * (sol[0] - m_points[1][0]) <= 0) {
			in_range = true;
		}
		for (size_t i = 0; i < sol.size(); ++i) {
			if (in_range) {
				if (i == 0) {
					attach_point.push_back(sol[i]);
					dim_dist.push_back(0.);
				}
				else {
					Real ref_point = calOneDimValue(sol[0] - m_points[0][0], i);
					if ((sol[i] < ref_point) && (sol[i] < ref_point + m_offset[i])) {
						attach_point.push_back(std::min(ref_point, ref_point + m_offset[i]));
						dim_dist.push_back(sol[i] - attach_point.back());
					}
					else if ((sol[i] > ref_point) && (sol[i] > ref_point + m_offset[i])) {
						attach_point.push_back(std::max(ref_point, ref_point + m_offset[i]));
						dim_dist.push_back(sol[i] - attach_point.back());
					}
					else {
						attach_point.push_back(sol[i]);
						dim_dist.push_back(0.);
					}
				}
			}
			else {
				if (m_points[0][0] < m_points[1][0]) {
					if (i == 0) {
						if (sol[i] < m_points[0][0]) {
							attach_point.push_back(m_points[0][0]);
						}
						else {
							attach_point.push_back(m_points[1][0]);
						}
						dim_dist.push_back(sol[i] - attach_point.back());
					}
					else {
						if (sol[i] < m_points[0][0]) {
							if ((sol[i] < m_points[0][i]) && (sol[i] < m_points[0][i] + m_offset[i])) {
								attach_point.push_back(std::min(m_points[0][i], m_points[0][i] + m_offset[i]));
								dim_dist.push_back(sol[i] - attach_point.back());
							}
							else if ((sol[i] > m_points[0][i]) && (sol[i] > m_points[0][i] + m_offset[i])) {
								attach_point.push_back(std::max(m_points[0][i], m_points[0][i] + m_offset[i]));
								dim_dist.push_back(sol[i] - attach_point.back());
							}
							else {
								attach_point.push_back(m_points[0][i]);
								dim_dist.push_back(0.);
							}
						}
						else {
							if ((sol[i] < m_points[1][i]) && (sol[i] < m_points[1][i] + m_offset[i])) {
								attach_point.push_back(std::min(m_points[1][i], m_points[1][i] + m_offset[i]));
								dim_dist.push_back(sol[i] - attach_point.back());
							}
							else if ((sol[i] > m_points[1][i]) && (sol[i] > m_points[1][i] + m_offset[i])) {
								attach_point.push_back(std::max(m_points[1][i], m_points[1][i] + m_offset[i]));
								dim_dist.push_back(sol[i] - attach_point.back());
							}
							else {
								attach_point.push_back(m_points[1][i]);
								dim_dist.push_back(0.);
							}
						}
					}
				}
				else {
					if (i == 0) {
						if (sol[i] < m_points[1][0]) {
							attach_point.push_back(m_points[1][0]);
						}
						else {
							attach_point.push_back(m_points[0][0]);
						}
						dim_dist.push_back(sol[i] - attach_point.back());
					}
					else {
						if (sol[0] < m_points[1][0]) {
							if ((sol[i] < m_points[1][i]) && (sol[i] < m_points[1][i] + m_offset[i])) {
								attach_point.push_back(std::min(m_points[1][i], m_points[1][i] + m_offset[i]));
								dim_dist.push_back(sol[i] - attach_point.back());
							}
							else if ((sol[i] > m_points[1][i]) && (sol[i] > m_points[1][i] + m_offset[i])) {
								attach_point.push_back(std::max(m_points[1][i], m_points[1][i] + m_offset[i]));
								dim_dist.push_back(sol[i] - attach_point.back());
							}
							else {
								attach_point.push_back(m_points[1][i]);
								dim_dist.push_back(0.);
							}
						}
						else {
							if ((sol[i] < m_points[0][i]) && (sol[i] < m_points[0][i] + m_offset[i])) {
								attach_point.push_back(std::min(m_points[0][i], m_points[0][i] + m_offset[i]));
								dim_dist.push_back(sol[i] - attach_point.back());
							}
							else if ((sol[i] > m_points[0][i]) && (sol[i] > m_points[0][i] + m_offset[i])) {
								attach_point.push_back(std::max(m_points[0][i], m_points[0][i] + m_offset[i]));
								dim_dist.push_back(sol[i] - attach_point.back());
							}
							else {
								attach_point.push_back(m_points[0][i]);
								dim_dist.push_back(0.);
							}
						}
					}
				}

			}
		}
		Real sum = 0;
		for (size_t i = 0; i < dim_dist.size(); i++) {
			sum = sum + std::pow(std::fabs(dim_dist[i]), p_norms[i]);
		}
		sum = std::sqrt(sum);
		return std::make_pair<>(attach_point, sum);
	}


	std::vector<Real> Line_shape::calValue(const Real& offset) {
		std::vector<Real> value;
		value.push_back(m_points[0][0]+offset);
		for (size_t i = 1; i < m_points[0].size(); ++i) {
			value.push_back(m_points[0][i] + (m_points[1][i]-m_points[0][i]) * (offset/m_offset[0]));
		}
		return value;
	}

	std::vector<Real> Trigo_shape::calValue(const Real& offset) {
		std::vector<Real> value;
		value.push_back(m_points[0][0] + offset);
		for (size_t i = 1; i < m_points[0].size(); ++i) {
			value.push_back(m_points[0][i] + (m_points[1][i] - m_points[0][i]) * (offset / m_offset[0]));
		}
		Real dist = sign(offset)*euclideanDistance(value.begin(),value.end(),m_points[0].begin());
		value.back() = value.back() + m_amp * std::sin(m_fre*dist);
		return value;
	}

	std::vector<Real> Power_shape::calValue(const Real& offset) {
		std::vector<Real> value;
		value.push_back(m_points[0][0] + offset);
		for (size_t i = 1; i < m_points[0].size(); ++i) {
			value.push_back(m_points[0][i] + (m_points[1][i] - m_points[0][i]) * (offset / m_offset[0]));
		}
		Real dist = euclideanDistance(value.begin(), value.end()-1, m_points[0].begin());
		value.back() = m_points[0].back() + m_scale * std::pow(dist,m_power);
		return value;
	}
}