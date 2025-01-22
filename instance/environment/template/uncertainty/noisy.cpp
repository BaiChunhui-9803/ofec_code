#include"noisy.h"

namespace ofec {
	Noisy::Noisy() {
		m_input_parameters.add("flagNoiseObjectives", new Bool(m_flag_noisy_from_objective, false));
		m_input_parameters.add("flagNoiseVariables", new Bool(m_flag_noisy_from_variable, false));
		m_input_parameters.add("flagNoiseEnvironments", new Bool(m_flag_noisy_from_environment, false));
	}

	void Noisy::initialize_() {
		m_objective_noisy_severity = 0.0;
		m_variable_noisy_severity = 0.0;
		m_environment_noisy_severity = 0.0;
	}

	void Noisy::copy(const Noisy &noisy) {
		Problem::copy(noisy);
		m_flag_noisy_from_objective = noisy.m_flag_noisy_from_objective;
		m_objective_noisy_severity = noisy.m_objective_noisy_severity;
		m_flag_noisy_from_variable = noisy.m_flag_noisy_from_variable;
		m_variable_noisy_severity = noisy.m_variable_noisy_severity;
		m_flag_noisy_from_environment = noisy.m_flag_noisy_from_environment;
		m_environment_noisy_severity = noisy.m_environment_noisy_severity;
	}
}
