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

#ifndef OFEC_DYNDESUBPOP_H
#define OFEC_DYNDESUBPOP_H

#include "dynde_Solution.h"
#include "../../../../template/classic/de/population.h"

namespace ofec {
	//class DynDE;
	class SubPopDynDE final : public PopDE<IndDynDE> {
	public:
		//friend class DynDE;
		explicit SubPopDynDE(size_t size, Problem *pro);
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
		bool getFlag() const { return m_flag; }
		void setFlag(bool flag) { m_flag = flag; }
	protected:
		void assign_type();
	protected:
		bool m_flag; // wheter need to be re-initialize, marked in the exclution_check()
		int m_num_normal;    // the number of normal Solutions of each swarm
		int m_num_brownian;       // the number of brownian Solutions of each swarm
		int m_num_quantum;    // the number of quantum Solutions of each swarm
		Real m_r_cloud; // radius of quantum swarms
		Real m_sigma;	// deviation for generation of Brownian Solutions
	};
}
#endif //OFEC_DYNDESUBPOP_H