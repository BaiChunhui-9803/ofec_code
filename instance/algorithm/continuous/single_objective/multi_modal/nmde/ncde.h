/********* Begin Register Information **********
{
	"name": "NCDE",
	"identifier": "NCDE",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li and Li Zhou
* Email: changhe.lw@gmail.com, 441837060@qq.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*********************************************************************************/
// updated Mar 28, 2018 by Li Zhou
// Updated on 15th Aug, 2019 by Junchen Wang
// Updated on 3rd Sep, 2024 by Junchen Wang

/******************************************************************************************
@article{qu2012differential,
  volume = {16},
  number = {5},
  journal = {IEEE Transactions on Evolutionary Computation},
  pages = {601--614},
  year = {2012},
  author = {Boyang Qu and Ponnuthurai Nagaratnam Suganthan and Jing Liang},
  title = {Differential evolution with neighborhood mutation for multimodal optimization}
}
******************************************************************************************/

#ifndef OFEC_NCDE_H
#define OFEC_NCDE_H

#include "../../../../../../core/algorithm/algorithm.h"
#include <list>

namespace ofec {
	class NCDE  : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(NCDE)
	protected:
		size_t m_pop_size;
	
		void addInputParameters();
		void initialize_(Environment* env) override;
		void run_(Environment *env) override;
	};
}


#endif // OFEC_NCDE_H

