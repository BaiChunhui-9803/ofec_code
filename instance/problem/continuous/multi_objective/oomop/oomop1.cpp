#include"oomop1.h"
#include <fstream>
#include<cmath>
#include"../../../../../utility/functional.h"

namespace ofec {
	void OOMOP1::initialize_() {
		OOMOP::initialize_();
		// assign variables
		auto& v = *m_param;
#ifdef OFEC_DEMO
		auto pri_var = v.get<std::string>("vector of private variables");
		m_num_pri_var = str2vector(pri_var, 1);
		auto pri_peak_num = v.get<std::string>("vector of private peaks");
		m_pri_peak_num = str2vector(pri_peak_num, 1);
		auto pri_peak_height_rate = v.get<std::string>("private peak height ratio");
		m_pri_peak_height_rate = str2vector(pri_peak_height_rate, 1.);
		auto pri_peak_slope_range = v.get<std::string>("private peak slope range");
		m_pri_peak_slope_range = str2vector(pri_peak_slope_range, 1.);
		auto pri_norm = v.get<std::string>("private dim norm");
		m_pri_norm = str2vector(pri_norm, 1.);
		auto pub_norm = v.get<std::string>("public dim norm");
		m_pub_norm = str2vector(pub_norm, 1.);
		auto global_ps_span = v.get<std::string>("vector of global PS span");
		m_global_ps_span = str2vector(global_ps_span, 1.);
		auto ps_decay_rate = v.get<std::string>("PS decay rate");
		m_ps_decay_rate = str2vector(ps_decay_rate, 1.);
		auto ps_radiate_slope = v.get<std::string>("PS radiate slope");
		m_ps_radiate_slope = str2vector(ps_radiate_slope, 1.);
		auto ps_type = v.get<std::string>("global PS shape type");
		m_ps_type = str2vector(ps_type, 1);
		auto pf_type = v.get<std::string>("global PF shape type");
		m_pf_type = str2vector(pf_type, 1);
#endif // ofec
		
#ifndef OFEC_DEMO
		m_num_pri_var = std::get<std::vector<int>>(v.at("vector of private variables"));
		m_pri_peak_num= std::get<std::vector<int>>(v.at("vector of private peaks"));
		m_pri_peak_height_rate= std::get<std::vector<Real>>(v.at("private peak height ratio"));
		m_pri_peak_slope_range = std::get<std::vector<Real>>(v.at("private peak slope range"));
		m_pri_norm= std::get<std::vector<Real>>(v.at("private dim norm"));
		m_pub_norm= std::get<std::vector<Real>>(v.at("public dim norm"));
		m_global_ps_span= std::get<std::vector<Real>>(v.at("vector of global PS span"));
		m_ps_decay_rate= std::get<std::vector<Real>>(v.at("PS decay rate"));
		m_ps_radiate_slope = std::get<std::vector<Real>>(v.at("PS radiate slope"));
		m_ps_type= std::get<std::vector<int>>(v.at("global PS shape type"));
		m_pf_type = std::get<std::vector<int>>(v.at("global PF shape type"));
#endif // ofec_demo

		m_num_ps = v.get<int>("number of PS shapes");
		m_num_global_ps = v.get<int>("number of global PS shapes");
		m_num_pub_var = v.get<int>("number of public variables");

		size_t total_var=0;
		for (size_t i = 0; i < m_num_pri_var.size(); ++i) {
			total_var += m_num_pri_var[i];
		}
		resizeVariable(total_var+m_num_pub_var);//recomend n=10
		// set domains
		setDomain(-1., 1.);

		/***********************************
		      construct private space
		************************************/ 
		size_t h = 20;
		//construct mpb for each obj
		size_t count = 0;
		for (size_t i = 0; i < numberObjectives(); ++i) {
			Matrix peak_pos(m_pri_peak_num[i], m_num_pri_var[i]);//pos of peaks
			for (size_t j = 0; j < m_num_pri_var[i]; ++j) {
				for (size_t k = 0; k < m_pri_peak_num[i]; ++k) {
					// random position
					peak_pos[k][j] = m_domain[count].limit.first + m_random->uniform.next() * (m_domain[count].limit.second - m_domain[count].limit.first);
					// designated position
				}
				count++;
			}
			std::vector<std::shared_ptr<Peak>> pri_peaks;
			for (size_t j = 0; j < m_pri_peak_num[i]; ++j) {
				Real pri_peak_h;
				if (m_num_pri_var[i] != 0) {
					if (j == 0) {
						pri_peak_h=h;
					}
					else {
						pri_peak_h=m_pri_peak_height_rate[i] * h + (1 - m_pri_peak_height_rate[i]) * h * m_random->uniform.next();
					}
				}
				else {
					pri_peak_h=0.;
				}
				Real slope = m_pri_peak_slope_range[0]+(m_pri_peak_slope_range[1]- m_pri_peak_slope_range[0])* m_random->uniform.next();
				Vector dim_norm(m_num_pri_var[i], 1);
				dim_norm = m_pri_norm[i] * dim_norm;
				Real exp_pri = 2;
				Real rotate_angle = 0;
				pri_peaks.emplace_back(std::make_shared<Linear_peak>(peak_pos[j].vect(), pri_peak_h, slope, dim_norm.vect(), exp_pri, rotate_angle));
			}
			m_mpb_pri.emplace_back(std::make_shared<Mpb_class>(pri_peaks));
		}

		/***********************************
		       construct public space
		************************************/
		std::vector<std::pair<Real, size_t>> ps_record;//record decay rate and class of ps
		for (size_t i = 0; i < m_num_ps; ++i) {
			std::pair<Real, size_t> ps_recordi;
			if (i < m_num_global_ps) {
				ps_recordi.first = 1.;
				ps_recordi.second = i;
			}
			else {
				size_t depend_idx = std::floor(m_num_global_ps * m_random->uniform.next());
				ps_recordi.first = m_ps_decay_rate[depend_idx]+(1- m_ps_decay_rate[depend_idx])* m_random->uniform.next();
				ps_recordi.second = depend_idx;
			}
			ps_record.emplace_back(ps_recordi);
		}
		m_ps_record = ps_record;
		std::vector<std::pair<Real, Real>> pub_dim_bounds;
		for (size_t j = 0; j < m_num_pub_var; ++j) {
			std::pair<Real, Real> dim_bound;
			dim_bound.first = m_domain[count + j].limit.first;
			dim_bound.second = m_domain[count + j].limit.second;
			pub_dim_bounds.emplace_back(dim_bound);
		}
		for (size_t i = 0; i < m_num_ps; ++i) {
			//offset of all public dims
			std::vector<Real> offset;
			for (size_t i = 0; i < m_num_pub_var; ++i) {
				if (i == 0) {
					offset.push_back(m_global_ps_span[ps_record[i].second]);
				}
				else {
					offset.push_back(0.);
				}
			}
			//boundary points
			std::vector<std::vector<Real>> boundary_points;
			Vector first_point;
			Vector second_point;
			//set ps decay ratio
			Real ps_ratio = ps_record[i].first;
			// pre generate boundary points
			for (size_t j = 0; j < m_num_pub_var; ++j) {
				first_point.pushBack(m_domain[count + j].limit.first + m_random->uniform.next() * (m_domain[count + j].limit.second - m_domain[count + j].limit.first));
			}
			boundary_points.emplace_back(first_point.vect());
			// determine slope and the second point by the first point,no sequential between start and end points;
			second_point = first_point;
			if (m_num_pub_var > 1) {
				// check boundary
				do {
					for (size_t j = 0; j < m_num_pub_var; ++j) {
						second_point[j] = m_domain[count + j].limit.first + m_random->uniform.next() * (m_domain[count + j].limit.second - m_domain[count + j].limit.first);
					}
					//scale length into set span and update the first offset
					auto diff_vector = second_point - first_point;
					second_point = first_point + m_global_ps_span[ps_record[i].second] / diff_vector.length() * diff_vector;
				} while (!checkBoundary(pub_dim_bounds, second_point.vect()));
				offset[0] = second_point[0] - first_point[0];
			}
			else {
				if (first_point[0] + offset[0] > m_domain[count].limit.second) {
					second_point[0] = first_point[0] - offset[0];
				}
				else {
					second_point[0] = first_point[0] + offset[0];
				}
			}
			boundary_points.emplace_back(second_point.vect());
			Real amplitude,temp_max,temp_min,period,temp_dist,power,scale;
			switch (m_ps_type[ps_record[i].second]) {
			case 1://line
				//set offset, check boundary
				if (m_num_pub_var > 2 || m_number_objectives>2) {
					std::vector<size_t> selected_dim;
					selected_dim.push_back(0);
					size_t off_num = std::min(m_number_objectives-1,m_num_pub_var);
					for (size_t j = 1; j < off_num; ++j) {
						//determin offset dim
						size_t idx = 0;
						do {
							idx = std::floor(1 + (m_num_pub_var - 1) * m_random->uniform.next());
						} while (std::find(selected_dim.begin(), selected_dim.end(), idx) != selected_dim.end());
						selected_dim.push_back(idx);
						offset[idx] = offset[0];//offset same dist
						//check boundary
						if (boundary_points[0][idx] + offset[idx] > pub_dim_bounds[idx].second || boundary_points[1][idx] + offset[idx] > pub_dim_bounds[idx].second) {
							offset[idx] = -1*offset[0];
						}
					}
				}
				m_ps_shapes.emplace_back(std::make_shared<Line_shape>(boundary_points, offset, ps_ratio));
				break;
			case 2://sin
				// revise end points according to amplitude
				amplitude = 0.1;
				temp_max = std::max(boundary_points[0].back(), boundary_points[1].back());
				temp_min = std::min(boundary_points[0].back(), boundary_points[1].back());
				if (temp_max + amplitude > pub_dim_bounds.back().second) {
					boundary_points[0].back() = boundary_points[0].back() - amplitude;
					boundary_points[1].back() = boundary_points[1].back() - amplitude;
				}
				if (temp_max - amplitude < pub_dim_bounds.back().first) {
					boundary_points[0].back() = boundary_points[0].back() + amplitude;
					boundary_points[1].back() = boundary_points[1].back() + amplitude;
				}
				// period
				period = 2. * OFEC_PI /euclideanDistance(boundary_points[0].begin(), boundary_points[0].end(),boundary_points[1].begin());
				//set offset,
				if (m_num_pub_var > 2 || m_number_objectives > 2) {
					std::vector<size_t> selected_dim;
					selected_dim.push_back(0);
					size_t off_num = std::min(m_number_objectives - 1, m_num_pub_var);
					for (size_t j = 1; j < off_num; ++j) {
						//determin offset dim
						size_t idx = 0;
						do {
							idx = std::floor(1 + (m_num_pub_var - 1) * m_random->uniform.next());
						} while (std::find(selected_dim.begin(), selected_dim.end(), idx) != selected_dim.end());
						selected_dim.push_back(idx);
						offset[idx] = offset[0];//offset same dist
						//check boundary
						if (idx == m_num_pub_var - 1) {
							Real temp_max = std::max(boundary_points[0][idx], boundary_points[1][idx]);
							if (temp_max + offset[idx] + amplitude > pub_dim_bounds[idx].second) {
								offset[idx] = -1 * offset[0];
							}
						}
						else {
							if (boundary_points[0][idx] + offset[idx] > pub_dim_bounds[idx].second || boundary_points[1][idx] + offset[idx] > pub_dim_bounds[idx].second) {
								offset[idx] = -1 * offset[0];
							}
						}
						
					}
				}
				m_ps_shapes.emplace_back(std::make_shared<Trigo_shape>(boundary_points,amplitude,period,offset, ps_ratio));
				break;
			case 3://power
				// calculate scale according to end points
				temp_dist = euclideanDistance(boundary_points[0].begin(), boundary_points[0].end()-1, boundary_points[1].begin());
				power = 4.;
				scale = (boundary_points[1].back() - boundary_points[0].back()) / std::pow(temp_dist,power);
				//set offset,
				if (m_num_pub_var > 2 || m_number_objectives > 2) {
					std::vector<size_t> selected_dim;
					selected_dim.push_back(0);
					size_t off_num = std::min(m_number_objectives - 1, m_num_pub_var);
					for (size_t j = 1; j < off_num; ++j) {
						//determin offset dim
						size_t idx = 0;
						do {
							idx = std::floor(1 + (m_num_pub_var - 1) * m_random->uniform.next());
						} while (std::find(selected_dim.begin(), selected_dim.end(), idx) != selected_dim.end());
						selected_dim.push_back(idx);
						offset[idx] = offset[0];//offset same dist
						//check boundary
						if (idx == m_num_pub_var - 1) {
							Real temp_max = std::max(boundary_points[0][idx], boundary_points[1][idx]);
							if (temp_max + offset[idx] > pub_dim_bounds[idx].second) {
								offset[idx] = -1 * offset[0];
							}
						}
						else {
							if (boundary_points[0][idx] + offset[idx] > pub_dim_bounds[idx].second || boundary_points[1][idx] + offset[idx] > pub_dim_bounds[idx].second) {
								offset[idx] = -1 * offset[0];
							}
						}
					}
				}
				m_ps_shapes.emplace_back(std::make_shared<Power_shape>(boundary_points, power, scale, offset, ps_ratio));
				break;
			default:
				break;
			}
		}

		/******************************************/
		/*   construct PS peaks in PS regions     */
		/******************************************/
		Real max_h = 20;
		Real min_h = 2;
		Real h_step = (max_h - min_h) / m_num_global_ps;
		for (size_t i = 0; i < m_num_ps; ++i) {
			//set anchor points
			std::vector<std::vector<Real>> anchor_points;
			for (size_t j = 0; j < m_ps_shapes[i]->getOffset().size();++j) {
				if (j == 0) {
					anchor_points.emplace_back(m_ps_shapes[i]->getPoints()[0]);
					anchor_points.emplace_back(m_ps_shapes[i]->getPoints()[1]);
				}
				else {
					if (m_ps_shapes[i]->getOffset()[j] != 0) {
						size_t temp_size = anchor_points.size();
						for (size_t k = 0; k < temp_size; ++k) {
							auto temp = anchor_points[k];
							temp[j] = temp[j] + m_ps_shapes[i]->getOffset()[j];
							anchor_points.push_back(temp);
						}
					}
				}
			}
			if (m_ps_shapes[i]->getRatio()==1.) {
				m_global_anchor_points.emplace_back(anchor_points);
			}
			//the max distance for each anchor point
			std::vector<Real> max_span;
			for (size_t j = 0; j < m_number_objectives-1; ++j) {
				std::vector<Real> dist;
				for (size_t k = 0; k < anchor_points.size(); ++k) {
					dist.push_back(euclideanDistance(anchor_points[j].begin(),anchor_points[j].end(),anchor_points[k].begin()));
				}
				max_span.push_back(maxElement(dist));
			}
			std::vector<std::shared_ptr<Peak>> ps_peaks;
			for (size_t j = 0; j < m_number_objectives - 1; ++j) {
				//peak location
				auto peak_pos = anchor_points[j];
				//peak height and slope
				Real ps_peak_h;
				Real ps_peak_slope;
				if ( ps_record[i].first== 1.) {
					ps_peak_h = max_h - h_step * i;
					ps_peak_slope = h_step / max_span[j];
				}
				else {
					ps_peak_h = m_ps_peaks[ps_record[i].second][j]->getHeight()* ps_record[i].first;
					ps_peak_slope = std::dynamic_pointer_cast<Linear_peak>(m_ps_peaks[ps_record[i].second][j])->getSlope();
				}
				std::vector<Real> ps_peak_norm=m_pub_norm;
				Real ps_exp = 2;
				Real rotate_angle = 0.;
				ps_peaks.emplace_back(std::make_shared<Linear_peak>(peak_pos,ps_peak_h,ps_peak_slope,ps_peak_norm,ps_exp,rotate_angle));
			}
			m_ps_peaks.emplace_back(ps_peaks);
		}

		/*************************************
		         generate PF samples
		*************************************/
		generatePF();
	}

	void OOMOP1::generatePF() {
		// get optima in private space
		std::vector<std::vector<Real>> optima_in_pri;
		for (size_t i = 0; i < m_number_objectives; ++i) {
			if (m_num_pri_var[i] > 0) {
				auto idx = m_mpb_pri[i]->getMaxPeakIndex();
				for (size_t j = 0; j < idx.size(); ++j) {
					auto x_pos = m_mpb_pri[i]->getPeak(idx[j])->getPeakPos();
					optima_in_pri.emplace_back(x_pos);
				}
			}
		}
		// get optima in public space
		size_t dim_num = 150;
		std::vector<std::vector<Real>> optima_in_pub=samplePubOptima(dim_num,m_num_global_ps);
		Real max_value = evaluate_opt_pub_sol(optima_in_pub,m_num_global_ps);
		setMaxValueLastObj(max_value);
		
		// composite sols and evaluate
		std::vector<Real> pri_optima;
		for (size_t j = 0; j < optima_in_pri.size(); ++j) {
			auto temp = optima_in_pri[j];
			for (size_t k = 0; k < temp.size(); ++k) {
				pri_optima.push_back(temp[k]);
			}
		}
		for (size_t i = 0; i < optima_in_pub.size(); ++i) {
			auto optima_sol = pri_optima;
			for (const auto& j : optima_in_pub[i]) {
				optima_sol.push_back(j);
			}
			std::vector<Real> temp_obj(m_number_objectives,0);
			evaluateObjective(optima_sol.data(), temp_obj);
			m_optima->appendObj(temp_obj);
		}
		m_optima->setObjectiveGiven(true);
		m_optima->setVariableGiven(true);

		for (size_t i = 0; i < m_number_objectives; ++i) {
			Real min_v = 1. * 10e14;
			Real max_v = -1. * 10e14;
			std::pair<Real, Real> obj_range;
			for (size_t j = 0; j < m_optima->numberObjectives(); ++j) {
				if (min_v > m_optima->objective(j)[i]) {
					min_v = m_optima->objective(j)[i];
				}
				if (max_v < m_optima->objective(j)[i]) {
					max_v = m_optima->objective(j)[i];
				}
			}
			obj_range.first = min_v;
			obj_range.second = max_v;
			m_range_optima.emplace_back(obj_range);
		}
	}

	void OOMOP1::evaluateObjective(Real* x, std::vector<Real>& obj) {
		/********************************************
		                 group variables
		*********************************************/
		std::vector<std::vector<Real>> seperated_variable;
		size_t count = 0;
		for (size_t i = 0; i < m_num_pri_var.size(); ++i) {
			std::vector<Real> pri_vars_for_obj;
			if (m_num_pri_var[i] > 0) {
				for (size_t j = 0; j < m_num_pri_var[i]; ++j) {
					pri_vars_for_obj.push_back(x[count]);
					count++;
				}
			}
			else {
				pri_vars_for_obj.push_back(0);
			}
			seperated_variable.push_back(pri_vars_for_obj);
		}
		std::vector<Real> pub_vars_for_obj;
		for (size_t i = 0; i < m_num_pub_var; ++i) {
			pub_vars_for_obj.push_back(x[count]);
			count++;
		}
		seperated_variable.push_back(pub_vars_for_obj);
		
		/***********************************************
		       calculate values in private space
		************************************************/
		std::vector<Real> pri_objs;
		std::vector<Real> max_peak_in_pri;
		for (size_t i = 0; i < m_number_objectives; ++i) {
			if (m_num_pri_var[i] > 0) {
				max_peak_in_pri.push_back(m_mpb_pri[i]->getPeak(m_mpb_pri[i]->getMaxPeakIndex()[0])->getHeight());
				pri_objs.push_back(m_mpb_pri[i]->calMpbValue(seperated_variable[i]));
			}
			else {
				max_peak_in_pri.push_back(0);
				pri_objs.push_back(0);
			}
		}
		
		/*****************************************
		     calculate values in public space 
		*****************************************/
		std::vector<Real> pub_objs(m_number_objectives,0.);
		std::vector<Real> max_peak_in_pub;
		for (size_t i = 0; i < m_number_objectives-1; ++i) {
			max_peak_in_pub.push_back(m_ps_peaks[0][i]->getHeight());
		}
		calPubValue(seperated_variable.back(),pub_objs);

		/*******************************************
		               composite objs
		*******************************************/
		for (size_t i = 0; i < m_number_objectives; ++i) {
			if (m_num_pri_var[i] > 0) {
				if (i != m_number_objectives - 1) {
					obj[i] = 2 - pri_objs[i] / max_peak_in_pri[i] - pub_objs[i] / max_peak_in_pub[i];
				}
				else {
					obj[i] = 2 - pri_objs[i] / max_peak_in_pri[i] - pub_objs[i] / getMaxValueLastObj();
				}
			}
			else{
				if (i != m_number_objectives - 1) {
					obj[i] = 1 - pub_objs[i] / max_peak_in_pub[i];
				}
				else {
					obj[i] = 1 - pub_objs[i] / getMaxValueLastObj();
				}
			}
		}
	}

	Real OOMOP1::evaluate_opt_pub_sol(std::vector<std::vector<Real>>& sols,size_t num_global_ps) {
		std::vector<Real> last_obj;
		std::vector<Real> max_peak_in_pub;
		for (size_t i = 0; i < m_number_objectives - 1; ++i) {
			max_peak_in_pub.push_back(m_ps_peaks[0][i]->getHeight());
		}
		for (size_t i = 0; i < sols.size(); ++i) {
			std::vector<Real> pub_objs;
			int gps_index = -1;
			for (size_t j = 0; j < num_global_ps; ++j) {
				if (m_ps_shapes[j]->if_in_range(sols[i])) {
					gps_index = j;
					break;
				}
			}
			for (size_t j = 0; j < m_number_objectives; ++j) {
				if (j != m_number_objectives - 1) {
					if (gps_index == -1) {
						MyExcept("solution is not in the designed PS");
					}
					pub_objs.push_back(m_ps_peaks[gps_index][j]->calPeakValue(sols[i]));
				}
				else {
					//calculate composite variable
					Real composite_x = 0;
					for (size_t k = 0; k < m_number_objectives - 1; ++k) {
						composite_x = composite_x + std::pow(pub_objs[k] / max_peak_in_pub[k], 2);
					}
					composite_x = std::sqrt(composite_x);
					Real alpha=0;
					Real beta=0;
					Real temp_k=0;
					switch (m_pf_type[gps_index])
					{
					case 1://linear:
						pub_objs.push_back(20 * (std::sqrt(m_number_objectives-1) - composite_x));
						break;
					case 2://convex:
						pub_objs.push_back(40 / (1 + std::exp(-5 * (std::sqrt(m_number_objectives - 1) - composite_x))) - 20);
						break;
					case 3://extreConvex:
						pub_objs.push_back(20 - 10 / (20 * (std::sqrt(m_number_objectives - 1) - composite_x) + 0.5));
						break;
					case 4://concave:
						pub_objs.push_back(40 - 38 / (1 + std::exp(-5 * composite_x)));
						break;
					case 5://extreConcave:
						pub_objs.push_back(10 / (20 * composite_x + 0.5));
						break;
					case 6://mix:
						pub_objs.push_back(20 / (1 + std::exp(-14 * (std::sqrt(m_number_objectives - 1)/2 - composite_x))));
						break;
					case 7://multiple_mix:
						alpha = 5;
						beta = 0.5;
						temp_k = 2 * (composite_x / beta - std::ceil(composite_x / beta)) + 1;
						pub_objs.push_back(20 * (1 - beta / 2 * sign(temp_k) * std::pow(std::fabs(temp_k), alpha) - beta * (std::ceil(composite_x / beta) - 1. / 2)));
						break;
					case 8://discontinuity:
						pub_objs.push_back(15 * (1.2 + 0.15 * std::cos(15 * composite_x) - composite_x/ std::sqrt(m_number_objectives - 1)));
						break;
					default:
						break;
					}
				}
			}
			last_obj.push_back(pub_objs.back());
		}
		return maxElement(last_obj);
	}

	void OOMOP1::calPubValue(std::vector<Real>& x, std::vector<Real>& pub_objs) {
		std::vector<Real> max_peak_in_pub;
		for (size_t i = 0; i < m_number_objectives - 1; ++i) {
			max_peak_in_pub.push_back(m_ps_peaks[0][i]->getHeight());
		}
		//determine the attached points and dist for each PS region
		int in_ps_index = -1;
		std::vector<std::vector<Real>> compare_obj;
		for (size_t i = 0; i < m_ps_shapes.size(); ++i) {
			auto temp_slope = m_ps_radiate_slope;
			std::vector<Real> p_norms=m_pub_norm;//norm=2 in the public space
			std::pair<std::vector<Real>, Real> attached_point =m_ps_shapes[i]->dist_from_PS(x, p_norms);
			//calculate objective values of attached points
			std::vector<Real> ref_objs;
			for (size_t j = 0; j < m_number_objectives; ++j) {
				if (j != m_number_objectives - 1) {
					ref_objs.push_back(m_ps_peaks[i][j]->calPeakValue(attached_point.first));
				}
				else {
					//calculate composite variable
					Real composite_x = 0;
					for (size_t k = 0; k < m_number_objectives - 1; ++k) {
						composite_x = composite_x + std::pow(ref_objs[k] / max_peak_in_pub[k], 2);
					}
					composite_x = std::sqrt(composite_x);
					if (m_ps_shapes[i]->getRatio() != 1.) {
						composite_x /= m_ps_shapes[i]->getRatio();
					}
					Real alpha = 5;
					Real beta = 0.5;
					Real temp_k = 0;
					switch (m_pf_type[m_ps_record[i].second])
					{
					case 1://linear:
						ref_objs.push_back(20 * (std::sqrt(m_number_objectives-1) - composite_x));
						break;
					case 2://convex:
						ref_objs.push_back(40 / (1 + std::exp(-5 * (std::sqrt(m_number_objectives - 1) - composite_x))) - 20);
						break;
					case 3://extreConvex:
						ref_objs.push_back(20 - 10 / (20 * (std::sqrt(m_number_objectives - 1) - composite_x) + 0.5));
						break;
					case 4://concave:
						ref_objs.push_back(40 - 38 / (1 + std::exp(-5 * composite_x)));
						break;
					case 5://extreConcave:
						ref_objs.push_back(10 / (20 * composite_x + 0.5));
						break;
					case 6://mix:
						ref_objs.push_back(20 / (1 + std::exp(-14 * (std::sqrt(m_number_objectives - 1)/2 - composite_x))));
						break;
					case 7://multiple_mix:
						alpha = 5;
						beta = 0.5;
						temp_k = 2 * (composite_x / beta - std::ceil(composite_x / beta)) + 1;
						ref_objs.push_back(20 * (1 - beta / 2 * sign(temp_k) * std::pow(std::fabs(temp_k), alpha) - beta * (std::ceil(composite_x / beta) - 1. / 2)));
						break;
					case 8://discontinuity:
						ref_objs.push_back(15 * (1.2 + 0.15 * std::cos(15 * composite_x) - composite_x/ std::sqrt(m_number_objectives - 1)));
						break;
					default:
						break;
					}
					if (m_ps_shapes[i]->getRatio() != 1.) {
						ref_objs[j] = ref_objs[j] * m_ps_shapes[i]->getRatio();
					}
				}
			}
			std::vector<Real> ps_obj(m_number_objectives);
			for (size_t j = 0; j < m_number_objectives; ++j) {
				ps_obj[j] = ref_objs[j] - temp_slope[i] * attached_point.second;
			}
			compare_obj.emplace_back(ps_obj);
			if (m_ps_shapes[i]->if_in_range(x)) {
				in_ps_index = i;
				break;
			}
		}
		if (in_ps_index != -1) {
			pub_objs = compare_obj[in_ps_index];
		}
		else if (m_ps_shapes.size() > 1) {
			size_t select_idx = 0;
			std::vector<Real> max_value_in_obj;//select the max value of objs as the baseline
			for (size_t i = 0; i < compare_obj.size(); ++i) {
				max_value_in_obj.push_back(maxElement(compare_obj[i]));
			}
			Real max_value = maxElement(max_value_in_obj);
			for (size_t j = 0; j < max_value_in_obj.size(); ++j) {
				if (max_value == max_value_in_obj[j]) {
					select_idx = j;
				}
			}
			pub_objs = compare_obj[select_idx];
		}
		else {
			pub_objs = compare_obj.back();
		}
	}

	std::vector<std::vector<Real>> OOMOP1::samplePubOptima(size_t dim_num, size_t num_g_ps) {
		std::vector<std::vector<Real>> optima_in_pub;
		for (size_t i = 0; i < num_g_ps; ++i) {
			std::vector<std::vector<Real>> optima_in_ps;
			auto points = m_ps_shapes[i]->getPoints();
			for (size_t j = 0; j < dim_num; ++j) {
				std::vector<Real> temp_point = points[0];
				for (size_t k = 0; k < temp_point.size(); ++k) {
					temp_point[k] = temp_point[k] + j * (points[1][k] - points[0][k]) / (dim_num - 1);
				}
				optima_in_ps.push_back(temp_point);
			}
			//
			switch (m_ps_shapes[i]->getShapeFlag())
			{
			case 1:

				break;
			case 2:
				for (size_t j = 1; j < optima_in_ps.size(); ++j) {
					Real temp_dist = sign(optima_in_ps[j][0]- optima_in_ps[0][0]) * euclideanDistance(optima_in_ps[j].begin(),optima_in_ps[j].end(),optima_in_ps[0].begin());
					Real temp_a = std::dynamic_pointer_cast<Trigo_shape>(m_ps_shapes[i])->getAmp();
					Real temp_f = std::dynamic_pointer_cast<Trigo_shape>(m_ps_shapes[i])->getFre();
					optima_in_ps[j].back() = optima_in_ps[j].back() + temp_a * std::sin(temp_f * temp_dist);
				}
				break;
			case 3:
				for (size_t j = 1; j < optima_in_ps.size(); ++j) {
					Real temp_dist = euclideanDistance(optima_in_ps[j].begin(), optima_in_ps[j].end()-1, optima_in_ps[0].begin());
					Real temp_scale = std::dynamic_pointer_cast<Power_shape>(m_ps_shapes[i])->getScale();
					Real temp_power = std::dynamic_pointer_cast<Power_shape>(m_ps_shapes[i])->getPower();
					optima_in_ps[j].back() = optima_in_ps[0].back() + temp_scale*std::pow(temp_dist,temp_power);
				}
				break;
			default:
				break;
			}
			Real step = euclideanDistance(points[0].begin(), points[0].end(), points[1].begin()) / (dim_num - 1);
			// scale in other dimension
			auto offset = m_ps_shapes[i]->getOffset();
			for (size_t j = 1; j < offset.size(); ++j) {
				if (std::fabs(offset[j]) > 0.0001) {
					size_t num = std::ceil(std::fabs(offset[j]) / step);
					sample2HighDim(optima_in_ps, offset[j], num, j);
				}
			}
			for (size_t j = 0; j < optima_in_ps.size(); ++j) {
				optima_in_pub.emplace_back(optima_in_ps[j]);
			}
		}
		return optima_in_pub;
	}
}