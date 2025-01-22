#include "uncertianty_continuous.h"
#include "../../../../../core/global.h"
#include "../../../../../utility/linear_algebra/vector.h"
#include "../../../../../core/problem/solution.h"
#include <algorithm>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif // OFEC_DEMO

namespace ofec {
	const std::vector<std::string> UncertaintyContinuous::ms_type = { "SmallStep", "LargeStep", "Random", "Recurrent", "Chaotic", "RecurrentNoisy" };

	void UncertaintyContinuous::resizeVariable(size_t num_vars)	{
		Continuous::resizeVariable(num_vars);
		for (int i(0); i < m_peak.size();++i) {
			m_peak[i].resize(m_number_variables);
			m_pre_peak[i].resize(m_number_variables);
			m_ini_peak[i].resize(m_number_variables);
		}
	}

	void UncertaintyContinuous::setNumChange(Real rRatio) {
		if (rRatio < 0 || rRatio>1) {
			throw MyExcept("the ratio of changing peaks is invalid@UncertaintyContinuous::set_num_change");
		}

		m_ratio_changing_peak = rRatio;
		m_num_changing_peak = (int)(m_peak.size()*m_ratio_changing_peak) > 1 ? (int)(m_peak.size()*m_ratio_changing_peak) : 1;
		updateNumChange();

		m_params["changePeakRatio"] = m_ratio_changing_peak;
	}

	void UncertaintyContinuous::updateNumChange() {
		if (m_num_changing_peak == m_peak.size()) {
			for (int i = 0; i < m_peak.size(); i++) m_flag_change[i] = true;
			return;
		}
		std::vector<int> a(m_peak.size());
		for (int i = 0; i < m_peak.size(); ++i) a[i] = i;
		m_random->uniform.shuffle(a.begin(), a.end());
//		global::ms_global->m_uniform[caller::Problem]->shuffle(a.begin(), a.end());
		// make sure the global optimum changes always
		int gopt = 0;
		for (int i = 0; i < m_peak.size(); i++) {
			if (m_flag_global_optima[i]) {
				gopt = i;
				break;
			}
		}
		int gidx;
		for (int i = 0; i < m_peak.size(); i++) {
			if (a[i] == gopt) {
				gidx = i;
				break;
			}
		}
		int t = a[0];
		a[0] = a[gidx];
		a[gidx] = t;

		for (int i = 0; i < m_peak.size(); i++) m_flag_change[i] = false;
		for (int i = 0; i < m_num_changing_peak; i++) m_flag_change[a[i]] = true;
		
	}

	void UncertaintyContinuous::setHeightSeverity(const Real rS) {
		for (int i = 0; i < m_peak.size(); i++) 	m_height_severity[i] = rS;
	}

	void UncertaintyContinuous::setWidthSeverity(const Real rS) {
		for (int i = 0; i < m_peak.size(); i++) 	m_width_severity[i] = rS;
	}

	void UncertaintyContinuous::setType(ChangeType rT) {
		m_change.type = rT;
		m_params["Change type"] = ms_type[static_cast<int>(m_change.type)];
	}

	//int UncertaintyContinuous::get_num_peak()const {
	//	return m_num_peaks;
	//}

	int UncertaintyContinuous::updateEvaluationTag(SolutionBase& s, Algorithm *alg) {
		return Continuous::updateEvaluationTag(s, alg) | Dynamic::updateEvaluationTag(s, alg);
	}

	void UncertaintyContinuous::updateCandidates(const SolutionBase& sol, std::list<std::unique_ptr<SolutionBase>>& candidates) const {
		if (candidates.empty()) {
			candidates.emplace_back(new Solution<>(sol));
		}
		if (!candidates.front() || sol.dominate(*candidates.front(), m_optimize_mode))
			candidates.front().reset(new Solution<>(sol));
	}


	//void UncertaintyContinuous::setType(ChangeType rT) {
	//	m_change.type = rT;
	//	m_params["Change type"] = ms_type[static_cast<int>(m_change.type)];
	//}

	//void UncertaintyContinuous::set_dimension_change(bool rFlag) {
	//	m_flag_dimension_change = rFlag;
	//	m_params["Flag dimensional change"] = m_flag_dimension_change;
	//}
	//void UncertaintyContinuous::set_noise_severity(Real value) {
	//	m_noise_severity = value;
	//	m_params["Noise severity"] = m_noise_severity;
	//}
	//void UncertaintyContinuous::set_time_linkage_severity(Real value) {
	//	m_time_linkage_severity = value;
	//	m_params["Time-linkage severity"] = m_time_linkage_severity;
	//}

	//void UncertaintyContinuous::set_change_dirction(bool rFlag) {
	//	m_direction_dimension_change = rFlag;
	//}
	//void UncertaintyContinuous::set_flag_synchronize(bool rFlag) {
	//	m_synchronize = rFlag;
	//}

	//void UncertaintyContinuous::set_recurrent_noisy_severity(Real rSeverity) {
	//	m_recurrent_noisy_severity = rSeverity;
	//	m_params["recurrent noisy severity"] = m_recurrent_noisy_severity;
	//}



	void UncertaintyContinuous::change() {
		Dynamic::change();
		switch (m_change.type) {
		case CT_Random:
			changeRandom();
			break;
		case CT_Recurrent:
			changeRecurrent();
			break;
		case CT_RecurrentNoisy:
			changeRecurrentNoisy();
			break;
		case CT_SmallStep:
			changeSmallStep();
			break;
		case CT_LargeStep:
			changeLargeStep();
			break;
		case CT_Chaotic:
			changeChaotic();
			break;
		default:
			break;
		}

		if (m_flag_variable_memory_change) {
			if (m_temp_dimension == m_min_dimension)
				m_flag_variable_memory_dir = 1;
			if (m_temp_dimension == m_max_dimension)
				m_flag_variable_memory_dir = -1;
			m_temp_dimension += m_flag_variable_memory_dir;
			changeVarMemory();
		}

		if (m_flag_num_peak_change) {
			if (m_mode == 1 || m_mode == 2) {
				if ((unsigned int)m_num_peaks >= m_max_peaks - 1) m_dir_num_peak_change = false;
				if ((unsigned int)m_num_peaks <= m_min_peaks + 1) m_dir_num_peak_change = true;
				int step = 2;

				if (m_mode == 2) {
					step = m_random->uniform.nextNonStd(step / 2, 5 * step / 2);
				}

				if (m_dir_num_peak_change == true) {
					if (m_num_peaks + step <= m_max_peaks)		m_temp_num_peak = m_num_peaks + step;
					else m_temp_num_peak = m_max_peaks;
				}
				else {
					if (m_num_peaks - step >= m_min_peaks)		m_temp_num_peak = m_num_peaks - step;
					else m_temp_num_peak = m_min_peaks;
				}
			}
			else {
				//random change
				m_temp_num_peak = m_random->uniform.nextNonStd(m_min_peaks, m_max_peaks);
			}
			changeNumPeak();
		}
	}


	void UncertaintyContinuous::copy(const Problem& rP) {
		Problem::copy(rP);
		Dynamic::copy(rP);
		Noisy::copy(rP);
		Continuous::copy(rP);

		auto& dcp = dynamic_cast<const UncertaintyContinuous &>(rP);
		m_change = dcp.m_change;
		m_recurrent_noisy_severity = dcp.m_recurrent_noisy_severity;
		m_alpha = dcp.m_alpha;
		m_max_alpha = dcp.m_max_alpha;
		m_chaotic_constant = dcp.m_chaotic_constant;
	//	m_noise_severity = dcp.m_noise_severity;
	//	m_time_linkage_severity = dcp.m_time_linkage_severity;


		//m_flag_dimension_change = dcp.m_flag_dimension_change;
		//m_direction_dimension_change = dcp.m_direction_dimension_change;
		m_synchronize = dcp.m_synchronize;


		m_flag_num_peak_change = dcp.m_flag_num_peak_change;
		m_dir_num_peak_change = dcp.m_dir_num_peak_change;
		m_mode = dcp.m_mode;



		int dim = m_temp_dimension < dcp.numberVariables() ? m_temp_dimension : dcp.numberVariables();
		int peaks = m_peak.size() < dcp.m_peak.size() ? m_peak.size() : dcp.m_peak.size();

		for (int i = 0; i < peaks; i++) {
			std::copy(dcp.m_peak[i].begin(), dcp.m_peak[i].begin()+ dim, m_peak[i].begin());
			std::copy(dcp.m_pre_peak[i].begin(), dcp.m_pre_peak[i].begin() + dim, m_pre_peak[i].begin());
			std::copy(dcp.m_ini_peak[i].begin(), dcp.m_ini_peak[i].begin() + dim, m_ini_peak[i].begin());
		}
		std::copy(dcp.m_height.begin(), dcp.m_height.begin()+ peaks, m_height.begin());
		std::copy(dcp.m_width.begin(), dcp.m_width.begin() + peaks, m_width.begin());
		std::copy(dcp.m_pre_height.begin(), dcp.m_pre_height.begin() + peaks, m_pre_height.begin());
		std::copy(dcp.m_pre_width.begin(), dcp.m_pre_width.begin() + peaks, m_pre_width.begin());
		std::copy(dcp.m_fitness.begin(), dcp.m_fitness.begin() + peaks, m_fitness.begin());
		std::copy(dcp.m_flag_change.begin(), dcp.m_flag_change.begin() + peaks, m_flag_change.begin());

		m_min_height = dcp.m_min_height;
		m_max_height = dcp.m_max_height;

		m_min_width = dcp.m_min_width;
		m_max_width = dcp.m_max_width;

		std::copy(dcp.m_height_severity.begin(), dcp.m_height_severity.begin() + peaks, m_height_severity.begin());
		std::copy(dcp.m_width_severity.begin(), dcp.m_width_severity.begin() + peaks, m_width_severity.begin());
		std::copy(dcp.m_flag_global_optima.begin(), dcp.m_flag_global_optima.begin() + peaks, m_flag_global_optima.begin());

		m_current_peak = dcp.m_current_peak;

		m_ratio_changing_peak = dcp.m_ratio_changing_peak;
		m_num_changing_peak = (int)(m_ratio_changing_peak*peaks) > 1 ? (int)(m_ratio_changing_peak*peaks) : 1;//dcp.m_num_changing_peak;

		std::copy(dcp.m_num_tracking.begin(), dcp.m_num_tracking.begin() + peaks, m_num_tracking.begin());
		std::copy(dcp.m_height_order.begin(), dcp.m_height_order.begin() + peaks, m_height_order.begin());
		std::copy(dcp.m_tracked.begin(), dcp.m_tracked.begin() + peaks, m_tracked.begin());
		m_num_peak_tracked = dcp.m_num_peak_tracked;
		std::copy(dcp.m_time_linkage.begin(), dcp.m_time_linkage.begin() + peaks, m_time_linkage.begin());

	}

	Real UncertaintyContinuous::getRecurrentNoise(int x, Real min, Real max, Real amplitude, Real angle, Real noisy_severity) {
		// return a value in recurrent with noisy dynamism environment
		Real y;
		Real noisy, t;
		y = min + amplitude * (sin(2 * OFEC_PI * (x + angle) / m_period) + 1) / 2.;
		noisy = noisy_severity * m_random->normal.next();
		t = y + noisy;
		if (t > min && t < max) y = t;
		else y = t - noisy;
		return y;
	}

	Real UncertaintyContinuous::chaoticStep(const Real x, const Real min, const Real max, const Real scale) {
		if (min > max) return -1.;
		Real chaotic_value;
		chaotic_value = (x - min) / (max - min);
		chaotic_value = m_chaotic_constant * chaotic_value * (1 - chaotic_value);
		return chaotic_value * scale;
	}

	bool UncertaintyContinuous::predictChange(int eff_evals, int evalsMore) {
		int fre = getFrequency();
		int evals = eff_evals % fre;
		if (evals + evalsMore >= fre) return true;
		else return false;
	}

	void UncertaintyContinuous::setNumPeakChangeMode(const int mode) {
		m_mode = mode;
	}

	int UncertaintyContinuous::getNumPeakChangeMode() {
		return m_mode;
	}

	//void UncertaintyContinuous::setFlagNoise(const bool flag) {
	//	m_flag_noise = flag;
	//	m_params["Flag of noise"] = flag;
	//}

	//void UncertaintyContinuous::setFlagTimeLinkage(const bool flag) {
	//	m_flag_time_linkage = flag;
	//	m_params["Flag time-linkage"] = m_flag_time_linkage;
	//}

	void UncertaintyContinuous::initialize_() {
		//Problem::initialize_();
		Dynamic::initialize_();
		Noisy::initialize_();
		Continuous::initialize_();

		auto& v = *m_param;

		m_num_peaks = v.get<int>("numPeak", 1);
		m_number_variables = v.get<int>("number of variables");
		m_frequency = v.get<int>("changeFre",1000);
		m_ratio_changing_peak = v.get<Real>("changePeakRatio",0.1);
		//if (v.has("changePeakRatio"))

		//else
		//	throw MyExcept("changePeakRatio is not given@UncertaintyContinuous::initialize_()");

		m_flag_variable_memory_change = v.has("flagNumDimChange") ? v.get<bool>("flagNumDimChange") : false;
		
		m_mode = v.has("change mode") ? v.get<int>("change mode") : 1;
		
		m_flag_noisy_from_variable = v.has("flagNoise") ? v.get<bool>("flagNoise") : false;
		if (m_flag_noisy_from_variable)
			m_var_noisy_severity = v.has("noiseSeverity") ? v.get<Real>("noiseSeverity") : 0.01;
		
		m_flag_time_linkage = v.has("flagTimeLinkage") ? v.get<bool>("flagTimeLinkage") : false;
		if (m_flag_time_linkage)
			m_time_linkage_severity = v.has("timelinkageSeverity") ? v.get<Real>("timelinkageSeverity") : 0.1;

		resizeVariable(m_number_variables);
		resizeObjective(1);

		m_max_dimension = 15;
		m_min_dimension = 2;     //should be greater than 1
		m_max_peaks = 100;
		m_min_peaks = 10;


		m_change.type = CT_Random;
		m_change.counter = 0;
		m_period = 0;
		m_flag_variable_memory_dir = 1;
		m_synchronize = true;
		m_recurrent_noisy_severity = 0.8;

		m_alpha = 0.04;
		m_max_alpha = 0.1;
		m_chaotic_constant = 3.67; //in [3.57,4]

		m_flag_num_peak_change = false;
		m_dir_num_peak_change = true;

		m_init_peaks = m_num_peaks;
		m_init_dimensions = m_number_variables;

		addTag(ProblemTag::kDOP);

		m_objective_accuracy = 0.01;
		m_peak.resize(m_num_peaks);
		m_pre_peak.resize(m_num_peaks);
		m_ini_peak.resize(m_num_peaks);
		for (int i = 0; i < m_num_peaks; i++) {
			m_peak[i].resize(m_number_variables);
			m_pre_peak[i].resize(m_number_variables);
			m_ini_peak[i].resize(m_number_variables);
		}

		m_width.resize(m_num_peaks);
		m_height.resize(m_num_peaks);
		m_pre_height.resize(m_num_peaks);
		m_pre_width.resize(m_num_peaks);
		m_fitness.resize(m_num_peaks);
		m_width_severity.resize(m_num_peaks);
		m_height_severity.resize(m_num_peaks);

		m_flag_change.resize(m_num_peaks);
		m_flag_global_optima.resize(m_num_peaks);
		m_num_tracking.resize(m_num_peaks);
		m_height_order.resize(m_num_peaks);
		m_tracked.resize(m_num_peaks);

		m_time_linkage.resize(m_num_peaks);
		for (int i = 0; i < m_num_peaks; ++i) {
			m_tracked[i] = false;
			m_flag_global_optima[i] = false;
			m_flag_change[i] = true;
			m_height_order[i] = -1;
		}

		setNumChange(m_ratio_changing_peak);

		Dynamic::m_flag_change = true;

	}

	void UncertaintyContinuous::calculateGlobalOptima() {
		Real gobj;
		if (m_optimize_mode[0] == OptimizeMode::kMaximize) gobj = *std::max_element(m_height.begin(), m_height.end());
		else gobj = *min_element(m_height.begin(), m_height.end());

		m_optima.reset(new Optima<>());;
		Real mindis = std::numeric_limits<Real>::max();
		for (int i = 0; i < m_peak.size(); i++) {
			m_flag_global_optima[i] = false;
			if (m_height[i] == gobj) {
				for (int j = 0; j < m_peak.size(); ++j) {
					if (j == i) continue;
					Vector s1(m_peak[i].vect()), s2(m_peak[j].vect());
					Real dis = s1.distance(s2);
					if (mindis > dis) {
						mindis = dis;
					}
				}
				m_flag_global_optima[i] = true;
				dynamic_cast<Optima<>&>(*m_optima).appendVar(VariableVector<Real>(m_peak[i]));
				m_optima->appendObj(m_height[i]);
			}
		}
		m_optima->setObjectiveGiven(true);
		m_optima->setVariableGiven(true);

		if (m_name=="DYN_CONT_RotationDBG" || m_name == "Moving-Peaks") {
			if (mindis / 2 < m_variable_accuracy)		
				m_variable_accuracy=(mindis / 2);
		}
		updateNumVisablePeak();
		m_num_peak_tracked = 0;
		if (m_flag_time_linkage) updateTimeLinkage();
		for (int i = 0; i < m_peak.size(); i++) {
			m_height_order[i] = i;
			m_tracked[i] = false;
		}
		mergeSort(m_height, m_peak.size(), m_height_order);		
		m_height_order=amend_order(m_height, m_height_order);
	}

	void UncertaintyContinuous::setHeight(const Real *h) {
		std::copy(h, h + m_peak.size(), m_height.begin());
	}


	//const Real **p

	void UncertaintyContinuous::setLocation(const std::vector<std::vector<Real>> &p) {
		for (int i = 0; i < m_peak.size(); i++) {
			m_peak[i].vect() = p[i];
		}
	}

	void UncertaintyContinuous::setInitialLocation(const std::vector<std::vector<Real>> &p) {
		for (int i = 0; i < m_peak.size(); i++) {
			m_ini_peak[i].vect() = m_peak[i].vect() = p[i];
		}
	}

	void UncertaintyContinuous::setWidth(const Real w) {
		for (int i = 0; i < m_peak.size(); i++)
			m_width[i] = w;
	}

	int UncertaintyContinuous::getNumVisiblePeak() {
		return m_num_visable_peak;

	}

	void UncertaintyContinuous::updateNumVisablePeak() {
		m_num_visable_peak = m_peak.size();
		for (int i = 0; i < m_peak.size(); i++) {
			if (!isVisible(i)) --m_num_visable_peak;
		}

	}

	bool UncertaintyContinuous::isVisible(const int rIdx) {
        Solution<VariableVector<Real>> s(m_number_objectives, m_number_constraints, m_number_variables);
        s.variable() = m_peak[rIdx];
		evaluate_(s);
  //      s.evaluate(false,caller::Problem);
		Real height = s.objective(0);
		switch (m_optimize_mode[0]) {
		case OptimizeMode::kMinimize:
			if (height < m_height[rIdx]) return false;
			break;
		case OptimizeMode::kMaximize:
			if (height > m_height[rIdx]) return false;
			break;
		}
		return true;

	}

	void UncertaintyContinuous::addNoise(Real *x_) {
		for (int d = 0; d < m_number_variables; d++) {
			Real x = x_[d];
			x += m_var_noisy_severity * m_random->normal.next();
			if (m_domain[d].limit.second < x) x = m_domain[d].limit.second;
			if (m_domain[d].limit.first > x)  x = m_domain[d].limit.first;
			x_[d] = x;
		}
	}

	bool UncertaintyContinuous::isTracked(Real *gen, Real obj) {
		bool flag = false, movepeaks = false;
		for (int i = 0; i < m_peak.size(); i++) {
			Real dis = 0, dis1 = fabs(obj - m_height[i]);
			for (int j = 0; j < m_number_variables; j++) dis += (gen[j] - m_peak[i][j])*(gen[j] - m_peak[i][j]);
			dis = sqrt(dis);
			if (dis <= m_variable_accuracy&&dis1 <= m_objective_accuracy) {
				// peak[i] assumed to be found
				int j = 0;
				while (m_height_order[j++] != i&&j < m_peak.size());
				if (!m_tracked[i]) {
					m_num_tracking[j - 1]++;
					m_tracked[i] = true;
					m_num_peak_tracked++;
					flag = true;
				}
			}
			if (dis < m_time_linkage_severity) {
				// move peak[i] to a near random position when it was tracked
				if (m_flag_time_linkage) {
					movePeak(i);
					updateTimeLinkage();
					movepeaks = true;
					m_flag_trigger_time_linkage = true;
				}
			}
		}
		if (movepeaks) {
#ifdef DEMON_OFEC
			calculateSamplePoints();
#endif
		}
		return flag;
	}

	int UncertaintyContinuous::getNumPeakFound() {
		return m_num_peak_tracked;
	}

	void UncertaintyContinuous::updateTimeLinkage() {
		if (!m_flag_time_linkage) return;
		Real range;
		for (int j = 0; j < m_number_variables; j++) {
			range = fabs(Real(m_domain[j].limit.second - m_domain[j].limit.first));
			m_time_linkage[j] = 0.2 * range * (m_random->uniform.next() - 0.5);

		}
	}

	void UncertaintyContinuous::movePeak( int idx) {
		if (idx < 0 || idx >= m_peak.size()) throw MyExcept("index out of boundary @ UncertaintyContinuous::move_peak(const int idx)");
		for (int d = 0; d < m_number_variables; d++) {
			Real x = m_peak[idx][d];
			x += m_time_linkage[d];
			if (m_domain[d].limit.second < x) x = m_domain[d].limit.second;
			if (m_domain[d].limit.first > x)  x = m_domain[d].limit.first;
			m_peak[idx][d] = x;
		}
	}

	bool UncertaintyContinuous::isGloOptTracked() {
		// the global optimum is assumed to be tracked if any one of the global optima is tracked	
		for (int i = 0; i < m_peak.size(); i++) {
			if (m_flag_global_optima[i] && m_tracked[i]) return true;
		}
		return false;
	}

	const std::vector<Real>& UncertaintyContinuous::getNearestPeak(const std::vector<Real>& p) {
		int nearest = 0;
		Vector peak(m_peak[0].vect());
		Real dis = peak.distance(p);
		for (int i = 1; i < m_peak.size(); i++) {
			std::copy(m_peak[i].begin(), m_peak[i].end(), peak.begin());
			Real d = peak.distance(p);
			if (d < dis) {
				dis = d;
				nearest = i;
			}
		}
		return m_peak[nearest].vect();
	}

	void UncertaintyContinuous::updatePrePeak() {
		m_pre_peak = m_peak;
		m_pre_height = m_height;
		m_pre_width = m_width;
	}
}