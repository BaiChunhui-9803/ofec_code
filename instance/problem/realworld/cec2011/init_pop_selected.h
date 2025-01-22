/********* Begin Register Information **********
[
	{ "name":"RW-CEC2011-T01-init-pop-selected", "identifier":"RW_CEC2011_T1_InitPopSelected", "problem tags":["ConOP","GOP","SOP"] },
	{ "name":"RW-CEC2011-T02-init-pop-selected", "identifier":"RW_CEC2011_T2_InitPopSelected", "problem tags":["ConOP","GOP","SOP"] },
	{ "name":"RW-CEC2011-T05-init-pop-selected", "identifier":"RW_CEC2011_T5_InitPopSelected", "problem tags":["ConOP","GOP","SOP"] },
	{ "name":"RW-CEC2011-T10-init-pop-selected", "identifier":"RW_CEC2011_T10_InitPopSelected", "problem tags":["ConOP","GOP","SOP"] }
]
*********** End Register Information **********/

#ifndef OFEC_RW_CEC2011_INIT_POP_SELECTED_H
#define OFEC_RW_CEC2011_INIT_POP_SELECTED_H

#include "../../init_pop_selected.h"
#include "t1_fm_soun_waves.h"
#include "t2_lj_potential.h"
#include "t5_tersoff_potential.h"
#include "t10_circular_antenna_array_design.h"

namespace ofec {
	using RW_CEC2011_T1_InitPopSelected = InitPopSelected<RW_CEC2011_T1>;
	using RW_CEC2011_T2_InitPopSelected = InitPopSelected<RW_CEC2011_T2>;
	using RW_CEC2011_T5_InitPopSelected = InitPopSelected<RW_CEC2011_T5>;
	using RW_CEC2011_T10_InitPopSelected = InitPopSelected<RW_CEC2011_T10>;
}

#endif // !OFEC_RW_CEC2011_INIT_POP_SELECTED_H
