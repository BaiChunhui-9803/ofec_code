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

*-------------------------------------------------------------------------------



*********************************************************************************/

#ifndef OFEC_MEMERY_RECORD_H
#define OFEC_MEMERY_RECORD_H

#include <vector>
#include <functional>
#include <fstream>

namespace ofec {
	class MemeryRecord {
	protected:
		int m_problem;
		double m_eps_distance = 0.3;
		std::vector<std::unique_ptr<SolBase>> m_memory;
		std::function<double(Problem *pro, const SolBase& sol)> m_fitness_fun;
		std::vector<int> m_com_relation;

	protected:
	public:

		//virtual void initialize(Problem *pro, double eps_distance=0.3);
		std::function<double(Problem *pro, const SolBase& sol)>& getEvalFun() {
			return m_fitness_fun;
		}
		void inputFile(const std::string& filename);
		void outputFile(const std::string& filename);
		void append(std::unique_ptr<SolBase>& p);
	};



}

#endif // !OFEC_NBC_H
