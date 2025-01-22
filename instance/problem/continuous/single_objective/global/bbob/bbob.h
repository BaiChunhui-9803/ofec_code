/********* Begin Register Information **********
[
	{ "name":"BBOB_F01", "identifier":"BBOB_F01", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F02", "identifier":"BBOB_F02", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F03", "identifier":"BBOB_F03", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F04", "identifier":"BBOB_F04", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F05", "identifier":"BBOB_F05", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F06", "identifier":"BBOB_F06", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F07", "identifier":"BBOB_F07", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F08", "identifier":"BBOB_F08", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F09", "identifier":"BBOB_F09", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F10", "identifier":"BBOB_F10", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F11", "identifier":"BBOB_F11", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F12", "identifier":"BBOB_F12", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F13", "identifier":"BBOB_F13", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F14", "identifier":"BBOB_F14", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F15", "identifier":"BBOB_F15", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F16", "identifier":"BBOB_F16", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F17", "identifier":"BBOB_F17", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F18", "identifier":"BBOB_F18", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F19", "identifier":"BBOB_F19", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F20", "identifier":"BBOB_F20", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F21", "identifier":"BBOB_F21", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F22", "identifier":"BBOB_F22", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F23", "identifier":"BBOB_F23", "tags":[ "continuous", "single-objective" ] },
	{ "name":"BBOB_F24", "identifier":"BBOB_F24", "tags":[ "continuous", "single-objective" ] }
]
*********** End Register Information **********/

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

/*
* source code refered to bbob.v15.02, http://coco.gforge.inria.fr/doku.php?id=bbob-2015-downloads
* @techreport{finck2010real,
  institution = {Citeseer},
  year = {2010},
  author = {Steffen Finck and Nikolaus Hansen and Raymond Ros and Anne Auger},
  title = {Real-parameter black-box optimization benchmarking 2009: Presentation of the noiseless functions}
}
*/

#ifndef OFEC_BBOB_H
#define OFEC_BBOB_H

#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../utility/linear_algebra/matrix.h"
#include "../../../../single_objective/metrics_gop.h"
#include <string>
#include <algorithm>

namespace ofec {
#define GET_BBOB(pro) dynamic_cast<BBOB*>(pro)
	class BBOB : public Continuous, public MetricsGOP {
		OFEC_ABSTRACT_INSTANCE(BBOB)
	public:
		void setAN(double alpha);

	protected:
		void addInputParameters();
		void initialize_(Environment *env) override;
		void updateOptima(Environment *env) override;
		std::vector<Real> computeXopt();
		Real computeFopt();
		Real penalize(Real *x);
		bool loadRotation(Real base_);
		void computeRotation(Matrix& rot, size_t Dim);
		void reshape(Matrix &B, std::vector<Real>& vector, size_t m, size_t n);
		void irregularize(Real *x);
		void asyemmetricalize(Real *x, Real belta);
		void evaluateObjective(Real *x, std::vector<Real>& obj) override;

	private:
		void f1_sphere(Real *x, std::vector<Real>& obj);
		void f2_ellipsoidale(Real *x, std::vector<Real>& obj);
		void f3_rastrigin(Real *x, std::vector<Real>& obj);
		void f4_buche_rastrigin(Real *x, std::vector<Real>& obj);
		void f5_slope(Real *x, std::vector<Real>& obj);
		void f6_sector(Real *x, std::vector<Real>& obj);
		void f7_step_ellipsoid(Real *x, std::vector<Real>& obj);
		void f8_original_rosenbrock(Real *x, std::vector<Real>& obj);
		void f9_rotated_rosenbrock(Real *x, std::vector<Real>& obj);
		void f10_nonseparable_ellipsoid(Real *x, std::vector<Real>& obj);
		void f11_discus(Real *x, std::vector<Real>& obj);
		void f12_bent_cigar(Real *x, std::vector<Real>& obj);
		void f13_sharp_ridge(Real *x, std::vector<Real>& obj);
		void f14_different_powers(Real *x, std::vector<Real>& obj);
		void f15_nonseparable_rastrigin(Real *x, std::vector<Real>& obj);
		void f16_weierstrass(Real *x, std::vector<Real>& obj);
		void f17_schaffers_F7(Real *x, std::vector<Real>& obj);
		void f18_illconditioned_schaffers_F7(Real *x, std::vector<Real>& obj);
		void f19_composite_griewank_rosenbrock(Real *x, std::vector<Real>& obj);
		void f20_schwefel(Real *x, std::vector<Real>& obj);
		void f21_gallagher_gaussian101me_peaks(Real *x, std::vector<Real>& obj);
		void f22_gallagher_gaussian21hi_peaks(Real *x, std::vector<Real>& obj);
		void f23_katsuura(Real *x, std::vector<Real>& obj);
		void f24_lunacekbi_rastrigin(Real *x, std::vector<Real>& obj);

	protected:
		typedef void (BBOB::*function_ptr)(Real *x, std::vector<Real>& obj);
		Real m_beta, m_alpha;
		function_ptr m_fun = nullptr;
		Matrix m_rot, m_rot2, m_linearTF;
		std::vector<Real> m_norRand;
		Real m_mu, m_scales, m_bias;
		Real m_condition_number;
	
		//for f21 and f22
		std::vector<Real> m_peak_values;
		std::vector < std::vector < Real >> m_arrScales;
		std::vector < std::vector < Real >> m_Xlocal;
		//for F16
		std::vector<Real> m_aK, m_bK;
		Real m_F0;

		size_t m_index_number;		// the index number in the benchmark suite
	};

	namespace bbob {
		int compareDoubles(const void *a, const void *b);
	}

	class BBOB_F01 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F01) protected: void addInputParameters() { m_index_number = 1; } };
	class BBOB_F02 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F02) protected: void addInputParameters() { m_index_number = 2; } };
	class BBOB_F03 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F03) protected: void addInputParameters() { m_index_number = 3; } };
	class BBOB_F04 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F04) protected: void addInputParameters() { m_index_number = 4; } };
	class BBOB_F05 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F05) protected: void addInputParameters() { m_index_number = 5; } };
	class BBOB_F06 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F06) protected: void addInputParameters() { m_index_number = 6; } };
	class BBOB_F07 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F07) protected: void addInputParameters() { m_index_number = 7; } };
	class BBOB_F08 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F08) protected: void addInputParameters() { m_index_number = 8; } };
	class BBOB_F09 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F09) protected: void addInputParameters() { m_index_number = 9; } };
	class BBOB_F10 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F10) protected: void addInputParameters() { m_index_number = 10; } };
	class BBOB_F11 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F11) protected: void addInputParameters() { m_index_number = 11; } };
	class BBOB_F12 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F12) protected: void addInputParameters() { m_index_number = 12; } };
	class BBOB_F13 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F13) protected: void addInputParameters() { m_index_number = 13; } };
	class BBOB_F14 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F14) protected: void addInputParameters() { m_index_number = 14; } };
	class BBOB_F15 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F15) protected: void addInputParameters() { m_index_number = 15; } };
	class BBOB_F16 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F16) protected: void addInputParameters() { m_index_number = 16; } };
	class BBOB_F17 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F17) protected: void addInputParameters() { m_index_number = 17; } };
	class BBOB_F18 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F18) protected: void addInputParameters() { m_index_number = 18; } };
	class BBOB_F19 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F19) protected: void addInputParameters() { m_index_number = 19; } };
	class BBOB_F20 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F20) protected: void addInputParameters() { m_index_number = 20; } };
	class BBOB_F21 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F21) protected: void addInputParameters() { m_index_number = 21; } };
	class BBOB_F22 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F22) protected: void addInputParameters() { m_index_number = 22; } };
	class BBOB_F23 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F23) protected: void addInputParameters() { m_index_number = 23; } };
	class BBOB_F24 : public BBOB { OFEC_CONCRETE_INSTANCE(BBOB_F24) protected: void addInputParameters() { m_index_number = 24; } };
}
#endif // !OFEC_BBOB_H
