/********* Begin Register Information **********
[
	{ "name":"RW_CEC2011_T01_BOUNDARY", "identifier":"RW_CEC2011_T1_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
	{ "name":"RW_CEC2011_T02_BOUNDARY", "identifier":"RW_CEC2011_T2_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
    { "name":"RW_CEC2011_T05_BOUNDARY", "identifier":"RW_CEC2011_T5_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
	{ "name":"RW_CEC2011_T06_BOUNDARY", "identifier":"RW_CEC2011_T6_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
	{ "name":"RW_CEC2011_T07_BOUNDARY", "identifier":"RW_CEC2011_T7_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
	{ "name":"RW_CEC2011_T09_BOUNDARY", "identifier":"RW_CEC2011_T9_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
	{ "name":"RW_CEC2011_T10_BOUNDARY", "identifier":"RW_CEC2011_T10_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
	{ "name":"RW_CEC2011_T11_BOUNDARY", "identifier":"RW_CEC2011_T11_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
	{ "name":"RW_CEC2011_T12_BOUNDARY", "identifier":"RW_CEC2011_T12_BOUNDARY", "tags":["ConOP","GOP","SOP"] },
	{ "name":"RW_CEC2011_T13_BOUNDARY", "identifier":"RW_CEC2011_T13_BOUNDARY", "tags":["ConOP","GOP","SOP"] }
]
*********** End Register Information **********/

//	{ "name":"CEC2011_MAT_F03", "identifier":"CEC2011_MAT_F03", "tags":[ "continuous", "single-objective" ] },

#ifndef OFEC_RW_CEC2011_BOUNDARY_H
#define OFEC_RW_CEC2011_BOUNDARY_H

#include "../../boundary_map.h"
#include "t1_fm_soun_waves.h"
#include "t2_lj_potential.h"
#include "t5_t6_tersoff_potential.h"
#include "t7_sprd_spectrum_rad_pphase.h"
#include "t9_large_scale_transmission_pricing.h"
#include "t10_circular_antenna_array_design.h"
#include "t11_eld_problem.h"
#include "t12_messengar.h"
#include "t13_cassini2.h"

namespace ofec {
	using RW_CEC2011_T1_BOUNDARY = ProblemBoundary<RW_CEC2011_T1>;
	using RW_CEC2011_T2_BOUNDARY = ProblemBoundary<RW_CEC2011_T2>;
	using RW_CEC2011_T5_BOUNDARY = ProblemBoundary<RW_CEC2011_T5>;
	using RW_CEC2011_T6_BOUNDARY = ProblemBoundary<RW_CEC2011_T6>;
	using RW_CEC2011_T7_BOUNDARY = ProblemBoundary<RW_CEC2011_T7>;
	using RW_CEC2011_T9_BOUNDARY = ProblemBoundary<RW_CEC2011_T9>;
	using RW_CEC2011_T10_BOUNDARY = ProblemBoundary<RW_CEC2011_T10>;
	using RW_CEC2011_T11_BOUNDARY = ProblemBoundary<RW_CEC2011_T11>;
	using RW_CEC2011_T12_BOUNDARY = ProblemBoundary<RW_CEC2011_T12>;
	using RW_CEC2011_T13_BOUNDARY = ProblemBoundary<RW_CEC2011_T13>;
}

#endif // !OFEC_RW_CEC2011_INIT_POP_SELECTED_H
