/********* Begin Register Information **********
{
	"name": "JADE",
	"identifier": "JADE",
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
// Fixed on Feb, 2024 by Junchen Wang (wangjunchen.chris@gmail.com)

/*
@article{zhang2009jade,
  title={JADE: adaptive differential evolution with optional external archive},
  author={Zhang, Jingqiao and Sanderson, Arthur C},
  journal={IEEE Transactions on Evolutionary Computation},
  volume={13},
  number={5},
  pages={945--958},
  year={2009},
  publisher={IEEE}
}
*/

#ifndef OFEC_JADE_H
#define OFEC_JADE_H

#include "jade_pop.h"
#include "../../../../../../core/algorithm/algorithm.h"

namespace ofec {
	class JADE : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(JADE)
	protected:
		std::unique_ptr<PopJADE<>> m_pop;
		size_t m_pop_size;

		void addInputParameters();
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	};
}
#endif // OFEC_JADE_H
