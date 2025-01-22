/********* Begin Register Information **********
[
	{ "name":"CEC2022_F01", "identifier":"CEC2022_F01", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F02", "identifier":"CEC2022_F02", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F03", "identifier":"CEC2022_F03", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F04", "identifier":"CEC2022_F04", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F05", "identifier":"CEC2022_F05", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F06", "identifier":"CEC2022_F06", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F07", "identifier":"CEC2022_F07", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F08", "identifier":"CEC2022_F08", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F09", "identifier":"CEC2022_F09", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F10", "identifier":"CEC2022_F10", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F11", "identifier":"CEC2022_F11", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2022_F12", "identifier":"CEC2022_F12", "tags":[ "continuous", "single-objective" ] }
]
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao & Changhe Li 
* Email: diaoyiyacug@163.com&changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

/*
* source code refered to CEC2022.v15.02, http://coco.gforge.inria.fr/doku.php?id=CEC2022-2015-downloads
* @techreport{finck2010real,
  institution = {Citeseer},
  year = {2010},
  author = {Steffen Finck and Nikolaus Hansen and Raymond Ros and Anne Auger},
  title = {Real-parameter black-box optimization benchmarking 2009: Presentation of the noiseless functions}
}
***********************************************************************/

#ifndef OFEC_CEC2022_COMPOSITION_H
#define OFEC_CEC2022_COMPOSITION_H

#include "../../../../../../core/instance.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../utility/linear_algebra/matrix.h"

namespace ofec {
		class CEC2022Function : virtual public Continuous {
			OFEC_ABSTRACT_INSTANCE(CEC2022Function);


		protected:
			//void ackley_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Ackley's */
			//void bent_cigar_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Discus */
			//void discus_func(double*, double* y, double* z, double*, int, double*, double*, int, int);  /* Bent_Cigar */
			//void ellips_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Ellipsoidal */
			//void escaffer6_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Expanded Scaffer¡¯s F6  */
			//void griewank_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Griewank's  */
			//void grie_rosen_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Griewank-Rosenbrock  */
			//void happycat_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* HappyCat */
			//void hgbat_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* HGBat  */
			//void rosenbrock_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Rosenbrock's */
			//void rastrigin_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Rastrigin's  */
			//void schwef 
			// 
			// el_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Schwefel's */
			//void schaffer_F7_func(double*, double* y, double* z, double*, int, double*, double*, int, int); /* Schwefel's F7 */
			//void step_rastrigin_func(double*, double*, int, double*, double*, int, int); /* Noncontinuous Rastrigin's  */
			//void levy_func(double*, double*, int, double*, double*, int, int); /* Levy */
			//void zakharov_func(double*, double*, int, double*, double*, int, int); /* ZAKHAROV */
			//void katsuura_func(double*, double*, int, double*, double*, int, int); /* Katsuura */



			void hf02(double*, double*, double*, double*, int, double*, double*, int*, int, int); /* Hybrid Function 2 */
			void hf06(double*, double*, double*, double*, int, double*, double*, int*, int, int); /* Hybrid Function 6 */
			void hf10(double*, double*, double*, double*, int, double*, double*, int*, int, int); /* Hybrid Function 10 */


			void cf01(double*, double*, double*, double*, int, double*, double*, int); /* Composition Function 1 */
			void cf02(double*, double*, double*, double*, int, double*, double*, int); /* Composition Function 2 */
			void cf06(double*, double*, double*, double*, int, double*, double*, int); /* Composition Function 6 */
			void cf07(double*, double*, double*, double*, int, double*, double*, int); /* Composition Function 7 */



			//void shiftfunc(double*, double*, int, double*);
			//void rotatefunc(double*, double*, int, double*);
			//void sr_func(double*, double*, int, double*, double*, double, int, int); /* shift and rotate */
			//void asyfunc(double*, double* x, int, double);
			//void oszfunc(double*, double*, int);
			//void cf_cal(double*, double*, int, double*, double*, double*, double*, int);


			template<typename T>
			void readVector(const std::string& filename, std::vector<T>& data) {
				std::ifstream in(filename);
				for (auto& it : data) {
					in >> it;
				}
				in.close();
			}

			void addInputParameters();
			void initialize_(Environment* env) override;
			void evaluateObjective(Real* x, std::vector<Real>& obj) override;
			void updateOptima(Environment* env) override;

			

		protected:
			std::vector<double> m_M;
			std::vector<double> m_OShift;
			std::vector<int> m_SS;

			double m_offsetY = 0;



			//extern int ini_flag, n_flag, func_flag, * SS;


			int m_index_number = 1;
		};

		class CEC2022_F01 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F01) protected: void addInputParameters() { m_index_number = 1; } };
		class CEC2022_F02 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F02) protected: void addInputParameters() { m_index_number = 2; } };
		class CEC2022_F03 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F03) protected: void addInputParameters() { m_index_number = 3; } };
		class CEC2022_F04 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F04) protected: void addInputParameters() { m_index_number = 4; } };
		class CEC2022_F05 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F05) protected: void addInputParameters() { m_index_number = 5; } };
		class CEC2022_F06 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F06) protected: void addInputParameters() { m_index_number = 6; } };
		class CEC2022_F07 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F07) protected: void addInputParameters() { m_index_number = 7; } };
		class CEC2022_F08 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F08) protected: void addInputParameters() { m_index_number = 8; } };
		class CEC2022_F09 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F09) protected: void addInputParameters() { m_index_number = 9; } };
		class CEC2022_F10 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F10) protected: void addInputParameters() { m_index_number = 10; } };
		class CEC2022_F11 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F11) protected: void addInputParameters() { m_index_number = 11; } };
		class CEC2022_F12 : public CEC2022Function { OFEC_CONCRETE_INSTANCE(CEC2022_F12) protected: void addInputParameters() { m_index_number = 12; } };
		
	
}
#endif // ! OFEC_CEC2015_COMPOSITION_H
