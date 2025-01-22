#ifndef DMOPs_H
#define DMOPs_H

#include "../../../../../core/problem/continuous/continuous.h"
#include <cmath>
#include "../../../../../utility/linear_algebra/matrix.h"
#include "../../../../../utility/linear_algebra/vector.h"
#include "../../../../../utility/functional.h"

namespace ofec {
	class DMOPs :virtual public Continuous {
	public:
		Real get_time();
		bool time_changed();
		size_t get_duration_gen() { return m_duration_gen; }
		size_t get_change_fre() { return m_change_fre; }
		size_t get_change_severity() { return m_severity; }
		void set_change_fre(size_t val) { m_change_fre = val; }
		void set_duration_gen(size_t val) { m_duration_gen = val; }
		void set_change_severity(Real val) { m_severity = val; }
		void set_updated_state(bool it) { updated_state = it; }
		bool get_updated_state() { return updated_state; }
		void set_effective_eval(size_t v) { m_evaluations = v; }
		size_t get_effective_eval() { return m_evaluations; }
		virtual void initialize_();
		virtual void generateAdLoadPF() {}
	private:
		size_t m_duration_gen;  //change frequency, every m_duration_gen generation
		Real m_severity;      //change severity
		size_t m_change_fre;     //change frequency, m_evaluations
		bool updated_state = false;//denote if updated problem
		size_t m_evaluations;
	};
}


#endif