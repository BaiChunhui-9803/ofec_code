/********* Begin Register Information **********
[
	{ "name":"MMOP_CEC2013_F01", "identifier":"MMOP_CEC2013_F01", "tags":[ "continuous", "single-objective" ] },
	{ "name":"MMOP_CEC2013_F02", "identifier":"MMOP_CEC2013_F02", "tags":[ "continuous", "single-objective" ] },
	{ "name":"MMOP_CEC2013_F03", "identifier":"MMOP_CEC2013_F03", "tags":[ "continuous", "single-objective" ] },
	{ "name":"MMOP_CEC2013_F04", "identifier":"MMOP_CEC2013_F04", "tags":[ "continuous", "single-objective" ] },
	{ "name":"MMOP_CEC2013_F05", "identifier":"MMOP_CEC2013_F05", "tags":[ "continuous", "single-objective" ] },
	{ "name":"MMOP_CEC2013_F06", "identifier":"MMOP_CEC2013_F06", "tags":[ "continuous", "single-objective" ] },
	{ "name":"MMOP_CEC2013_F07", "identifier":"MMOP_CEC2013_F07", "tags":[ "continuous", "single-objective" ] },
	{ "name":"MMOP_CEC2013_F08", "identifier":"MMOP_CEC2013_F08", "tags":[ "continuous", "single-objective" ] }
]
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
******************************************************************************************
*  Paper: Benchmark Functions for CEC��2013 Special Session and Competition on Niching
*  Methods for Multimodal Function Optimization.
*******************************************************************************************/

#ifndef OFEC_CEC13_F1_F8_H
#define OFEC_CEC13_F1_F8_H

#include "../classical/five_uneven_peak_trap.h"
#include "../classical/equal_maxima.h"
#include "../classical/uneven_de_maxima.h"
#include "f4_himmenblau.h"
#include "f5_six_hump_camel_back.h"
#include "f6_shubert.h"
#include "../classical/vincent.h"
#include "f8_modified_rastrigin.h"

namespace ofec {
	using MMOP_CEC2013_F01 = FiveUnevenPeakTrap;
	using MMOP_CEC2013_F02 = EqualMaxima;
	using MMOP_CEC2013_F03 = UnevenDeMaxima;
	using MMOP_CEC2013_F04 = cec2013::Himmenblau;
	using MMOP_CEC2013_F05 = cec2013::SixHumpCamelBack;
	using MMOP_CEC2013_F06 = cec2013::Shubert;
	using MMOP_CEC2013_F07 = Vincent;
	using MMOP_CEC2013_F08 = cec2013::ModifiedRastrigin;
}
#endif // !OFEC_CEC13_F1_F8_H