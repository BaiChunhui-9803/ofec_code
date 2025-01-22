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
*************************************************************************/

#ifndef OFEC_CEC2013_FUNCTION_H
#define OFEC_CEC2013_FUNCTION_H

#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../../utility/matrix.h"

namespace ofec {
	struct index_map {
		unsigned arr_index1;
		unsigned arr_index2;
	};
	class function_CEC2013 : public continuous {
	protected:
		unsigned ID;
		void create_shifted_vector(std::vector<Real> &vec);
		void create_random_permutation(std::vector<size_t> &vec);
		void create_rotated_matrix(size_t dim, matrix & mat);
		void create_rotated_matrix_1D(matrix & mat, std::vector<Real> & vec);
		void create_multi_rotated_matrix_1D(size_t dim, size_t num, std::vector<std::vector<Real>> &vvec);

		
		// Basic mathematical functions' declaration
		Real* multiply(Real * vector, Real * matrix, size_t dim);
		Real* multiply(Real * vector, Real ** matrix, size_t dim);
		Real elliptic(Real * x, size_t dim);
		
		Real rastrigin(Real * x, size_t dim);
		Real rastrigin(Real * x, size_t dim, size_t k);
		Real ackley(Real * x, size_t dim);
		Real ackley(Real * x, size_t dim, size_t k);
		
		Real rot_rastrigin(Real * x, size_t dim);
		Real rot_rastrigin(Real * x, size_t dim, size_t k);
		Real rot_ackley(Real * x, size_t dim);
		Real rot_ackley(Real * x, size_t dim, size_t k);
		Real schwefel(Real * x, size_t dim);
		Real schwefel(Real * x, size_t dim, size_t k);
		Real sphere(Real * x, size_t dim);
		Real sphere(Real * x, size_t dim, size_t k);
		Real rosenbrock(Real * x, size_t dim);
		Real rosenbrock(Real * x, size_t dim, size_t k);
		unsigned convertMatrixToArrayIndex(unsigned i, unsigned j);
		void create_index_mapping();
	
		Real* rotate_vector(size_t i, size_t &c);
		Real* rotate_vector_conform(size_t i, size_t &c);
		Real* rotate_vector_conflict(size_t i, size_t &c, Real* x);

		void setGlobalOpt(Real *tran = 0);
		void setOriginalGlobalOpt(Real *opt = 0);

		Real* mp_Ovector;
		size_t* mp_Pvector;
		Real* mp_anotherz;
		Real* mp_anotherz1;
		Real* mp_anotherz2;
		Real* mp_rot_matrix;
		size_t* mp_s;
		Real* mp_w;

		Real** mpp_OvectorVec;
		Real** mpp_MultiRotmatrix1D;
		Real** mpp_r25;
		Real** mpp_r50;
		Real** mpp_r100;

		bool mb_setOvectorToZero;

		size_t m_nonSeparableGroupNumber;
		size_t m_nonSeparableGroupSize;
		size_t m_overlap;

		std::vector<bool> mbv_interArray;
		struct index_map* mp_index_map;
		unsigned m_arrSize;

	public:
		//CEC2013(const ParameterMap &v);
		function_CEC2013(const std::string &name, size_t size_var, size_t size_obj);
		virtual ~function_CEC2013();

		std::vector<bool> getInterArray();
		void ArrToMat(unsigned I1, unsigned I2, unsigned &matIndex);
		void MatToArr(unsigned &I1, unsigned &I2, unsigned matIndex);

		/* for CEC2013SS */
		Real* readOvector();
		Real** readOvectorVec();
		size_t* readPermVector();
		Real** readR(size_t sub_dim);
		size_t* readS(size_t num);
		Real* readW(size_t num);

		size_t* getS();

		void transform_osz(Real* z, size_t dim);
		void transform_asy(Real* z, Real beta, size_t dim);
		void lambda(Real* z, Real alpha, size_t dim);
		int sign(Real x);
		Real hat(Real x);
		Real c1(Real x);
		Real c2(Real x);

		optima<VariableVector<Real>, Real> m_original_global_opt;
	};
}
#endif  // !OFEC_CEC2013_FUNCTION_H
