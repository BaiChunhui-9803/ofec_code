/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* class Algorithm is an abstract for all algorithms.
*
*********************************************************************************/
#ifndef OFEC_DISTANCE_CALCULATOR_GENERATOR_H
#define OFEC_DISTANCE_CALCULATOR_GENERATOR_H


#include<memory>
#include "../../template/framework/uncertianty/distance_calculator.h"
#include "distance_calculator_weight.h"


namespace ofec {
	template<typename TSolution>
	void generateDistanceCalculator(
		std::unique_ptr<>& calculator, const std::string& name) {

	}
	
}

#endif