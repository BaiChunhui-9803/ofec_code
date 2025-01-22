//Register rotation_DBG "DRPBG" ConOP,DOP,SOP

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
* 1. Changhe Li, Shengxiang Yang, D. Pelta. Benchmark Generator for CEC'2012 Competition on Evolutionary Computation for
* Dynamic Optimization Problems,Tech. Rep., the School of Computer Science, China University of Geosciences, Wuhan, China, 2012.
*
*********************************************************************************/
// Created: 11 May 2011
// Last modified: 4 August 2018

#ifndef ROTATIONDBG_H
#define ROTATIONDBG_H

#include "../real_DBG.h"

namespace ofec {
	class rotation_DBG final : public real_DBG
	{
	public:
		rotation_DBG(const std::string &name,  int rDimNumber,  int rNumPeaks, int fre,\
			 dynamic::ChangeType rT, Real  rChangingRatio,  bool rFlagDimChange,  bool rFlagNumPeakChange, \
			 int peakNumChangeMode,  bool flagNoise,  bool flagTimelinkage);
		rotation_DBG(const ParameterMap &v);
		void  set_width(Real w);
		void evaluateObjective(Real* x, std::vector<Real>& obj);
	protected:
		void change_width();
		virtual void change_random();
		virtual void change_small_step();
		virtual void change_large_step();
		virtual void change_recurrent();
		virtual void change_chaotic();
		virtual void change_recurrent_noisy();
		virtual void change_dimension();
		virtual void change_num_peak();
		void initialize();
	};
}
#endif // ROTATIONDBG_H
