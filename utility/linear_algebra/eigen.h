/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li and Junchen Wang
* Email: changhe.lw@gmail.com and wangjunchen.chris@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*-------------------------------------------------------------------------------
* 
* Some Eigen-based operations of matrices and vectors,
* which enable users to write Matlab-like codes in OFEC.
* 
*********************************************************************************/

#ifndef OFEC_EIGEN_H
#define OFEC_EIGEN_H

#include <Eigen/Dense>
#include "../../core/random/newran.h"

namespace Eigen {
	using ArrayXb = Array<bool, Dynamic, 1>;
}

namespace ofec {
	Eigen::MatrixXd randMatXd(int rows, int cols, Uniform *rnd);
	
	Eigen::VectorXd randVecXd(int size, Uniform *unif);

	Eigen::VectorXd randnVecXd(int size, Normal *norm);

	Eigen::VectorXi join(const Eigen::Ref<const Eigen::VectorXi> &vec, int val);
	
	Eigen::VectorXi join(const Eigen::Ref<const Eigen::VectorXi> &vec1, const Eigen::Ref<const Eigen::VectorXi> &vec2);
	
	Eigen::VectorXd join(const Eigen::Ref<const Eigen::VectorXd> &vec, double val);
	
	Eigen::VectorXd join(const Eigen::Ref<const Eigen::VectorXd> &vec1, const Eigen::Ref<const Eigen::VectorXd> &vec2);

	Eigen::VectorXi setdiff(const Eigen::Ref<const Eigen::VectorXi> &vec1, const Eigen::Ref<const Eigen::VectorXi> &vec2);

	double prctile(const Eigen::Ref<const Eigen::VectorXd> &vec, size_t p);

	double median(const Eigen::Ref<const Eigen::VectorXd> &vec);

	Eigen::VectorXi sort(const Eigen::Ref<const Eigen::VectorXd> &vec, bool ascend = true);

	Eigen::VectorXi find(const Eigen::Ref<const Eigen::ArrayXb> &vec);

	Eigen::VectorXd normcdf(const Eigen::Ref<const Eigen::VectorXd> &vec);
	
	Eigen::MatrixXd joinRow(const Eigen::Ref<const Eigen::MatrixXd> &mat1, const Eigen::Ref<const Eigen::MatrixXd> &mat2);

	Eigen::MatrixXd joinCol(const Eigen::Ref<const Eigen::MatrixXd> &mat1, const Eigen::Ref<const Eigen::MatrixXd> &mat2);
}

#endif // !OFEC_EIGEN_H
