/********* Begin Register Information **********
{
	"name": "MOP_ZDT3",
	"identifier": "ZDT3",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 31 December 2014
// Modified: 29 Mar 2018 by Junchen Wang (wangjunchen@cug.edu.cn)
//f2 is less than 0; thus, add 1 to f2 by tanqingshan

#ifndef OFEC_ZDT3_H
#define OFEC_ZDT3_H

#include "zdt.h"

namespace ofec {
	class ZDT3 : public ZDT {
	public:
		void updateOptima() override;
		void sampleParetoSols(size_t sample_num);

		/*void sampleParetoFront(size_t sample_num);
		void loadParetoFront(size_t sample_num);*/
	protected:
		void evaluateObjective(Real *x, std::vector<Real> &obj) override;
	};
}

#endif //OFEC_ZDT3_H