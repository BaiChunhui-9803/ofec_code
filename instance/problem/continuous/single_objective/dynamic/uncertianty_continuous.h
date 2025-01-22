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
* 1. Changhe Li and Shengxiang Yang. A Generalized Approach to Construct Benchmark Problems for Dynamic Optimization.
* 2008. Berlin, Heidelberg: Springer Berlin Heidelberg.
* 2. Changhe Li, Shengxiang Yang, D. Pelta. Benchmark Generator for CEC'2012 Competition on Evolutionary Computation for
* Dynamic Optimization Problems,Tech. Rep., the School of Computer Science, China University of Geosciences, Wuhan, China, 2012.
*
*********************************************************************************/
// Created: 11 May 2011
// modified: 30 July 2018
// modified: 10 June 2021 by DYY

#ifndef OFEC_UNCERTAIN_DYNAMIC_CONT_H
#define OFEC_UNCERTAIN_DYNAMIC_CONT_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../core/problem/uncertainty/dynamic.h"
#include "../../../../../core/problem/uncertainty/noisy.h"

#include "../../../../../core/problem/solution.h"



namespace ofec {
#define GET_DYNCONOP(pro) dynamic_cast<UncertaintyContinuous*>(pro)

	class UncertaintyContinuous : public Dynamic, public Noisy, public Continuous {
	public:
		enum ChangeType { CT_SmallStep = 0, CT_LargeStep, CT_Random, CT_Recurrent, CT_Chaotic, CT_RecurrentNoisy };
		struct Change {
			ChangeType type;
			int counter;
		};

	protected:
		const static std::vector<std::string> ms_type;
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
		Real m_alpha, m_max_alpha;              // to control Step severity
		Real m_chaotic_constant;

		// features below added on NOV 22 2012
		int m_mode;		// for the number of peaks change; 1: periodic with fixed Step, 2: periodic with random Step, 3: chaotic change


		std::vector<VariableVector<Real>> m_peak;					    	// positions of local or global optima(local optima in Rotation_DBG,
		std::vector<VariableVector<Real>> m_pre_peak;							// global optima of basic function in Composition_DBG)
		std::vector<VariableVector<Real>> m_ini_peak;				        // save the initial positions
		std::vector<Real> m_height;							// peak height in Rotation_DBG, height of global optima in Composition_DBG
		std::vector<Real> m_width;                           // weight value of each basic function in Composition_DBG,  peak width in Rotation_DBG

		///TODO preHeight and preWidth not considered in current version
		std::vector<Real> m_pre_height;
		std::vector<Real>m_pre_width;

		Real m_min_height, m_max_height;		// minimum\maximum height of all peaks(local optima) in Rotation_DBG(Composition_DBG)
		std::vector<Real> m_height_severity;
		Real m_min_width, m_max_width;
		std::vector<Real> m_width_severity;

		std::vector<Real> m_fitness;						    	// objective value of each basic funciton in Composition_DBG, peak height in Rotation_DBG
		std::vector<bool> m_flag_global_optima;                      // the index of the global optimal peak

		int m_current_peak;                         // the peak where the best Solution is located
		std::vector<bool> m_flag_change;                      // whether peaks change or not
		int m_num_changing_peak;                       // the number of peaks that change
		Real m_ratio_changing_peak = 1.0;                    // the ratio of changing peaks

		int m_num_visable_peak;                      // number of visable peaks, a peak is visable only if no peak is on top of it
		std::vector<int> m_num_tracking;							// accumulated number of peak[i] being tracked
		std::vector<int> m_height_order;
		int m_num_peak_tracked;
		std::vector<bool> m_tracked;
		std::vector<Real> m_time_linkage;



		void set(const UncertaintyContinuous& uc) {

			m_change = uc.m_change;
			m_recurrent_noisy_severity = uc.m_recurrent_noisy_severity;
			m_synchronize = uc.m_synchronize;
			m_temp_dimension = uc.m_temp_dimension;
			m_num_peaks = uc.m_num_peaks;
			m_flag_num_peak_change = uc.m_flag_num_peak_change;


			m_dir_num_peak_change = uc.m_dir_num_peak_change;
			m_temp_num_peak = uc.m_temp_num_peak;
			m_max_dimension = uc.m_max_dimension;
			m_min_dimension = uc.m_min_dimension;
			m_max_peaks = uc.m_max_peaks;
			m_min_peaks = uc.m_min_peaks;

			m_init_peaks = uc.m_init_peaks;
			m_init_dimensions = uc.m_init_dimensions;
			m_alpha = uc.m_alpha;
			m_max_alpha = uc.m_max_alpha;
			m_chaotic_constant = uc.m_chaotic_constant;
			m_mode = uc.m_mode;


			m_peak = uc.m_peak;
			m_pre_peak = uc.m_pre_peak;
			m_ini_peak = uc.m_ini_peak;
			m_height = uc.m_height;
			m_width = uc.m_width;

			m_pre_height = uc.m_pre_height;
			m_pre_width = uc.m_pre_width;
			m_min_height = uc.m_min_height;
			m_max_height = uc.m_max_height;
			m_height_severity = uc.m_height_severity;


			m_min_width = uc.m_min_width;
			m_max_width = uc.m_max_width;
			m_width_severity = uc.m_width_severity;
			m_fitness = uc.m_fitness;
			m_flag_global_optima = uc.m_flag_global_optima;
			m_current_peak = uc.m_current_peak;

			m_flag_change = uc.m_flag_change;
			m_num_changing_peak = uc.m_num_changing_peak;
			m_ratio_changing_peak = uc.m_ratio_changing_peak;
			m_num_visable_peak = uc.m_num_visable_peak;
			m_num_tracking = uc.m_num_tracking;
			m_height_order = uc.m_height_order;

			


			m_num_peak_tracked = uc.m_num_peak_tracked;
			m_tracked = uc.m_tracked;
			m_time_linkage = uc.m_time_linkage;

			if (uc.m_cur_working_solution == nullptr) {
				m_cur_working_solution.reset(nullptr);
			}
			else {
				m_cur_working_solution.reset(new Solution<>
					(dynamic_cast<const Solution<>&>(*uc.m_cur_working_solution)));
			}
		}

	public:
		UncertaintyContinuous() = default;
		UncertaintyContinuous(const UncertaintyContinuous& uc):
			Dynamic(uc), Noisy(uc), Continuous(uc){
			set(uc);
		//	updateOptBase();
			
		}
		UncertaintyContinuous& operator=(const UncertaintyContinuous& rhs) {
			if (this == &rhs)return *this;
			Dynamic::operator=(rhs);
			Noisy::operator=(rhs);
			Continuous::operator=(rhs);
			set(rhs);
		//	updateOptBase();
			return *this;
		}



		virtual ~UncertaintyContinuous() = default;


		virtual void resizeVariable(size_t num_vars);

		void setType(ChangeType rT);
		void setFlagNumPeaksChange(bool rPC) {
			m_flag_num_peak_change = rPC;
			m_params["Flag number of peaks change"] = m_flag_num_peak_change;
		}
		bool getFlagNumPeaksChange() { return m_flag_num_peak_change; }
		void setFlagSynchronize(bool rFlag) { m_synchronize = rFlag; }
		void setRecurrentNoisySeverity(Real rSeverity) { m_recurrent_noisy_severity = rSeverity; }

		void setAlpha(Real rAlpha) { m_alpha = rAlpha; }
		void setMaxAlpha(Real rMaxAlpha) { m_max_alpha = rMaxAlpha; }
		void setChaoticConstant(Real rValue) { m_chaotic_constant = rValue; }

		ChangeType getType() const { return m_change.type; }
		bool getFlagSynchronizeChange()const { return m_synchronize; }
		void setNumPeakChangeMode(int mode);
		int getNumPeakChangeMode();
		int getNumPeak()const { return m_num_peaks; }

		virtual void change();
		Real getRecurrentNoise(int x, Real min, Real max, Real amplitude, Real angle, Real noisy_severity = 1.);
		Real chaoticStep(Real x, Real min, Real max, Real scale = 1.0);
		bool predictChange(int eff_evals, int evalsMore);

		void setNumChange(Real rRatio);
		void setHeightSeverity(const Real rS);
		void setWidthSeverity(const Real rS);
		void setHeight(const Real *h);
		void setLocation(const std::vector<std::vector<Real>> &);
		void setInitialLocation(const std::vector<std::vector<Real>> &);
		virtual void setWidth(const Real w);
		int getNumVisiblePeak();
		bool isVisible(int rIdx);
		bool isTracked(Real *gen, Real obj);// is any peak tracked for the first time
		int getNumPeakFound();

		//15-07-2013
		bool isGloOptTracked();
		const std::vector<Real> &getNearestPeak(const std::vector<Real> &);

		virtual int updateEvaluationTag(SolutionBase &s, Algorithm *alg)override;
		void updateCandidates(const SolutionBase &sol, std::list<std::unique_ptr<SolutionBase>> &candidates) const override;

	protected:
		virtual void initialize_() override;
		virtual void changeRandom() {};
		virtual void changeSmallStep() {};
		virtual void changeLargeStep() {};
		virtual void changeRecurrent() {};
		virtual void changeChaotic() {};
		virtual void changeRecurrentNoisy() {};
		virtual void changeNumPeak() {};

	protected:
		void calculateGlobalOptima();
		void updateNumChange();
		void updateNumVisablePeak();
		void addNoise(Real *x);
		void updateTimeLinkage();
		void movePeak(int idx);
		void copy(const Problem & );
		void updatePrePeak();
	};
}
#endif // DYNAMICCONTINUOUS_H
