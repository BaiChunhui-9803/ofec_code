/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (ofec)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of ofec. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/ofec for more information
*
*-------------------------------------------------------------------------------
* Interface of the GDBG
*
*
*********************************************************************************/

/*
1. Changhe Li and Shengxiang Yang. A Generalized Approach to Construct Benchmark Problems for Dynamic Optimization.
2008. Berlin, Heidelberg: Springer Berlin Heidelberg.
*/

// Created: 19 June 2018
// Last modified:
#ifndef GDBG_H
#define GDBG_H
#include "../../../../core/problem/uncertainty/dynamic.h"

namespace ofec {
#define CAST_DYN dynamic_cast<GDBG*>(global::ms_global->m_problem.get())

	class GDBG : virtual public Dynamic
	{
	public:
		enum ChangeType { CT_SmallStep = 0, CT_LargeStep, CT_Random, CT_Recurrent, CT_Chaotic, CT_RecurrentNoisy };
		struct Change {
			ChangeType type;
			int counter;
		};
	protected:
		const static std::vector<std::string> ms_type;
	protected:
		Change m_change;
		Real m_recurrent_noisy_severity;			// deviation servity from the trajactory of recurrent change
		bool m_synchronize;                 // default=true all dimensions change at the same time

		int m_temp_dimension;                //a temporal variable for dimensional change only

		int m_num_peaks;
		bool m_flag_num_peak_change;                  // flag of the change of the number of peaks
		bool m_dir_num_peak_change;                   // true for increasing the number of peaks, otherwise decreasing the number of peaks
		int m_temp_num_peak;                         // temporal varibal for number of peaks change only


		unsigned m_max_dimension;
		unsigned m_min_dimension;     //should be greater than 1

		int m_max_peaks;
		int m_min_peaks;

		int m_init_peaks, m_init_dimensions;
		Real m_alpha, m_max_alpha;              // to control step severity
		Real m_chaotic_constant;

		// features below added on NOV 22 2012
		int m_mode;		// for the number of peaks change; 1: periodic with fixed step, 2: periodic with random step, 3: chaotic change


		bool m_flag_trigger_time_linkage = false;


	public:
		Dynamic() = default;
		Dynamic(const std::string &name, size_t num_peaks, size_t number_objectives = 1, size_t num_cons = 0);
		virtual ~Dynamic() = 0;
		Dynamic& operator=(const Dynamic &rhs ) = default;


		void setType(ChangeType rT);
		void setFlagNumPeaksChange(bool rPC) {
			m_flag_num_peak_change = rPC;
		}
		bool getFlagNumPeaksChange() {
			return m_flag_num_peak_change;
		}
		void setFlagSynchronize(bool rFlag) {
			m_synchronize = rFlag;
		}
		void setRecurrentNoisySeverity(Real rSeverity) {
			m_recurrent_noisy_severity = rSeverity;
		}

		void setAlpha(Real rAlpha) {
			m_alpha = rAlpha;
		};
		void setMaxAlpha(Real rMaxAlpha) {
			m_max_alpha = rMaxAlpha;
		};
		void setChoaticConstant(Real rValue) {
			m_chaotic_constant = rValue;
		}

		int getFrequency()const {
			return m_frequency;
		};
		int getChangeCounter()const {
			return m_counter;
		};
		int getPeriod()const {
			return m_period;
		}
		ChangeType getType() const {
			return m_change.type;
		};
		bool getFlagDimensionChange() const {
			return m_flag_dimension_change;
		};
		bool getDirDimensionChange() const {
			return m_direction_dimension_change;
		};
		bool getFlagSynchronizeChange()const {
			return m_synchronize;
		};

		void setNumPeakChangeMode(int mode);
		int getNumPeakChangeMode();
		void setFlagNoise(bool flag);
		int getNumPeak()const {
			return m_num_peaks;
		}
		void setFlagTimeLinkage(bool flag);
		bool getFlagTimeLinkage() const {
			return m_flag_time_linkage;
		}
		virtual void change();
		Real getRecurrentNoise(int x, Real min, Real max, Real amplitude, Real angle, Real noisy_severity = 1.);
		Real chaoticStep(Real x, Real min, Real max, Real scale = 1.0);
		bool predictChange(int evalsMore);
		void setNoiseSeverity(Real value) {
			m_noise_severity = value;
		}
		void setTimeLinkageSeverity(Real value) {
			m_time_linkage_severity = value;
		}
		bool& getTriggerFlagTimeLinkage() {
			return m_flag_trigger_time_linkage;
		}

		virtual void changeDimension() {};
	protected:
		virtual void changeRandom() {};
		virtual void changeSmallStep() {};
		virtual void changeLargeStep() {};
		virtual void changeRecurrent() {};
		virtual void changeChaotic() {};
		virtual void changeRecurrentNoisy() {};

		virtual void changeDimension();
		virtual void changeNumPeak() {};


		void copy(const Problem& rP);
	
	};
}

#endif