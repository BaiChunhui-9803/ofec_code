#include "../../../../core/global.h"
#include "dynamic.h"

#ifdef OFEC_DEMO
#include "../../../../../ui/buffer/scene.h"
//extern std::unique_ptr<ofec_demo::scene> ofec_demo::msp_buffer;
//extern bool ofg_algTermination;
#endif

namespace ofec {


	const std::vector<std::string> Dynamic::ms_type = { "SmallStep", "LargeStep", "Random", "Recurrent", "Chaotic", "RecurrentNoisy" };

	Dynamic::Dynamic(const std::string &name ,size_t num_peaks,size_t number_objectives = 1, size_t num_cons = 0) :Problem(name, number_objectives, num_cons), m_counter(0)
		,  m_num_peaks(num_peaks), m_temp_num_peak(num_peaks), m_flag_noise(false), m_flag_time_linkage(false), m_flag_trigger_time_linkage(false) {
		//ctor
		m_max_dimension = 15;
		m_min_dimension = 2;     //should be greater than 1
		m_max_peaks = 100;
		m_min_peaks = 10;
		m_frequency = 5000;
		m_change.type = CT_Random;
		m_change.counter = 0;
		m_period = 0;
		m_flag_dimension_change = false;
		m_direction_dimension_change = true;
		m_synchronize = true;
		m_recurrent_noisy_severity = 0.8;

		m_alpha = 0.04;
		m_max_alpha = 0.1;
		m_chaotic_constant = 3.67; //in [3.57,4]

		m_flag_num_peak_change = false;
		m_dir_num_peak_change = true;
		m_mode = 1;
		m_noise_severity = 0.01;
		m_time_linkage_severity = 0.1;


		m_init_peaks = m_num_peaks;
		m_init_dimensions = 1;
		
		addTag(ProblemTag::DOP);
	}

	Dynamic::~Dynamic() {
		//dtor
	}

	void Dynamic::setFrequency(int rChangeFre) {
		m_frequency = rChangeFre;
		m_params["Change frequency"] = m_frequency;
	}

	void Dynamic::setPeriod( int rPeriod) {
		m_period = rPeriod;
		m_params["Period"] = m_period;
	}

	void Dynamic::setType( ChangeType rT) {
		m_change.type = rT;
		m_params["Change type"] = ms_type[static_cast<int>(m_change.type)];
	}


	void Dynamic::change() {
		m_counter++;
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



		if (m_flag_num_peak_change) {
			if (m_mode == 1 || m_mode == 2) {
				if ((unsigned int)m_num_peaks >= m_max_peaks - 1) m_dir_num_peak_change = false;
				if ((unsigned int)m_num_peaks <= m_min_peaks + 1) m_dir_num_peak_change = true;
				int step = 2;

				if (m_mode == 2) {
					step=m_random->uniform.nextNonStd(step / 2, 5 * step / 2);
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

#ifdef OFEC_DEMO
//		msp_buffer->updateFitnessLandsacpe_();
#endif
	}

	void Dynamic::changeDimension(int cur_dim)
	{

		if (m_flag_dimension_change) {
			if (cur_dim == m_min_dimension)
				m_direction_dimension_change = true;
			if (cur_dim == m_max_dimension)
				m_direction_dimension_change = false;

			if (m_direction_dimension_change == true) {
				m_temp_dimension += 1;
		}
			else {
				m_temp_dimension -= 1;
			}
			changeDimension();
	}
	}

	void  Dynamic::copy(const Problem &rdp) {
		Problem::copy(rdp);

		const Dynamic &rDP = dynamic_cast<const Dynamic&>(rdp);
		m_change = rDP.m_change;
		m_frequency = rDP.m_frequency;
		m_period = rDP.m_period;

		m_flag_noise = rDP.m_flag_noise;
		m_flag_time_linkage = rDP.m_flag_time_linkage;
		m_flag_trigger_time_linkage = rDP.m_flag_trigger_time_linkage;

		m_recurrent_noisy_severity = rDP.m_recurrent_noisy_severity;
		m_alpha = rDP.m_alpha;
		m_max_alpha = rDP.m_max_alpha;
		m_chaotic_constant = rDP.m_chaotic_constant;


		m_noise_severity = rDP.m_noise_severity;
		m_time_linkage_severity = rDP.m_time_linkage_severity;


		m_flag_dimension_change = rDP.m_flag_dimension_change;
		m_direction_dimension_change = rDP.m_direction_dimension_change;
		m_synchronize = rDP.m_synchronize;


		m_flag_num_peak_change = rDP.m_flag_num_peak_change;
		m_dir_num_peak_change = rDP.m_dir_num_peak_change;
		m_mode = rDP.m_mode;


	}

	Real Dynamic::getRecurrentNoise(int x, Real min, Real max, Real amplitude, Real angle, Real noisy_severity) {
		// return a value in recurrent with noisy dynamism environment
		Real y;
		Real noisy, t;
		y = min + amplitude*(sin(2 * OFEC_PI*(x + angle) / m_period) + 1) / 2.;
		noisy = noisy_severity * m_random->normal.next();
		t = y + noisy;
		if (t > min&&t < max) y = t;
		else y = t - noisy;
		return y;
	}

	Real Dynamic::chaoticStep(const Real x, const Real min, const Real max, const Real scale) {
		if (min > max) return -1.;
		Real chaotic_value;
		chaotic_value = (x - min) / (max - min);
		chaotic_value = m_chaotic_constant*chaotic_value*(1 - chaotic_value);
		return chaotic_value*scale;
	}

	bool Dynamic::predictChange(const int evalsMore) {
		int fre = getFrequency();
		int evals = numEvals() % fre;
		if (evals + evalsMore >= fre) return true;
		else return false;
	}
	void Dynamic::setNumPeakChangeMode(const int mode) {
		m_mode = mode;
	}
	int Dynamic::getNumPeakChangeMode() {
		return m_mode;
	}
	void Dynamic::setFlagNoise(const bool flag) {
		m_flag_noise = flag;
		m_params["Flag of noise"] = flag;
	}
	void Dynamic::setFlagTimeLinkage(const bool flag) {
		m_flag_time_linkage = flag;
		m_params["Flag time-linkage"] = m_flag_time_linkage;
	}
}