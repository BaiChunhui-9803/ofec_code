#include "bbob.h"
#include "../../../../../../utility/functional.h"

namespace ofec {
	constexpr size_t NHIGHPEAKS21 = 101;
	constexpr size_t NHIGHPEAKS22 = 21;
	
	void BBOB::setAN(double alpha) {
		for (size_t i = 0; i < 12; ++i) { /// number of summands, 20 in CEC2005, 10/12 saves 30% of time
			m_aK[i] = pow(alpha, (Real)i);
			m_F0 += m_aK[i] * cos(2 * OFEC_PI * m_bK[i] * 0.5);
		}
	}

	void BBOB::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
	}

	void BBOB::initialize_(Environment *env) {
		Continuous::initialize_(env);
		resizeVariable(m_number_variables);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		setDomain(-5, 5);
		m_bias = computeFopt();
		if (m_index_number == 1) {
			m_fun = &BBOB::f1_sphere;
		}
		else if (m_index_number == 2) {
			m_condition_number = 1e6;
			m_fun = &BBOB::f2_ellipsoidale;
		}
		else if (m_index_number == 3) {
			m_condition_number = 10;
			m_beta = 0.2;
			m_fun = &BBOB::f3_rastrigin;
		}
		else if (m_index_number == 4) {
			m_condition_number = 10.;
			m_alpha = 100.;
			m_fun = &BBOB::f4_buche_rastrigin;
		}
		else if (m_index_number == 5) {
			m_alpha = 100.;
			m_fun = &BBOB::f5_slope;
			for (size_t i = 0; i < m_number_variables; ++i) {
				Real tmp = pow(sqrt(m_alpha), ((Real)i) / ((Real)(m_number_variables - 1)));
				m_bias += 5. * tmp;
			}
		}
		else if (m_index_number == 6) {
			m_fun = &BBOB::f6_sector;
			m_alpha = 100.;
			m_condition_number = 10.;
			loadRotation(sqrt(m_condition_number));
		}
		else if (m_index_number == 7) {
			m_fun = &BBOB::f7_step_ellipsoid;
			m_condition_number = 100.;
			m_alpha = 10.;
			loadRotation(sqrt(m_condition_number));
		}
		else if (m_index_number == 8) {
			m_fun = &BBOB::f8_original_rosenbrock;
			m_scales = fmax(1., sqrt((Real)m_number_variables) / 8.);
		}
		else if (m_index_number == 9) {
			m_fun = &BBOB::f9_rotated_rosenbrock;
			m_scales = fmax(1., sqrt((Real)m_number_variables) / 8.);
			computeRotation(m_rot, m_number_variables);
			m_linearTF.resize(m_number_variables, m_number_variables);
			for (size_t i = 0; i < m_number_variables; ++i){
				for (size_t j = 0; j < m_number_variables; ++j)
					m_linearTF[i][j] = m_scales * m_rot[i][j];
			}
		}
		else if (m_index_number == 10) {
			m_fun = &BBOB::f10_nonseparable_ellipsoid;
			computeRotation(m_rot, m_number_variables);
			m_condition_number = 1e6;
		}
		else if (m_index_number == 11) {
			m_fun = &BBOB::f11_discus;
			computeRotation(m_rot, m_number_variables);
			m_condition_number = 1e6;
		}
		else if (m_index_number == 12) {
			m_fun = &BBOB::f12_bent_cigar;
			computeRotation(m_rot, m_number_variables);
			m_condition_number = 1e6;
			m_beta = 0.5;
		}
		else if (m_index_number == 13) {
			m_fun = &BBOB::f13_sharp_ridge;
			m_condition_number = 10.;
			m_alpha = 100.;
			loadRotation(sqrt(m_condition_number));
		}
		else if (m_index_number == 14) {
			m_fun = &BBOB::f14_different_powers;
			m_alpha = 4.;
			computeRotation(m_rot, m_number_variables);
		}
		else if (m_index_number == 15) {
			m_fun = &BBOB::f15_nonseparable_rastrigin;
			m_condition_number = 10.;
			m_beta = 0.2;
			loadRotation(sqrt(m_condition_number));
		}
		else if (m_index_number == 16) {
			m_fun = &BBOB::f16_weierstrass;
			m_condition_number = 100.;
			m_F0 = 0.;
			m_aK.resize(12);
			m_bK.resize(12);
			loadRotation(1. / sqrt(m_condition_number));
			for (size_t i = 0; i < 12; ++i) /// number of summands, 20 in CEC2005, 10/12 saves 30% of time
			{
				m_aK[i] = pow(0.5, (Real)i);
				m_bK[i] = pow(3., (Real)i);
				m_F0 += m_aK[i] * cos(2 * OFEC_PI * m_bK[i] * 0.5);
			}
		}
		else if (m_index_number == 17) {
			m_fun = &BBOB::f17_schaffers_F7;
			computeRotation(m_rot, m_number_variables);
			computeRotation(m_rot2, m_number_variables);
			m_condition_number = 10.;
			m_beta = 0.5;
		}
		else if (m_index_number == 18) {
			m_fun = &BBOB::f18_illconditioned_schaffers_F7;
			computeRotation(m_rot, m_number_variables);
			computeRotation(m_rot2, m_number_variables);
			m_condition_number = 1e3;
			m_beta = 0.5;
		}
		else if (m_index_number == 19) {
			m_fun = &BBOB::f19_composite_griewank_rosenbrock;
			m_scales = fmax(1., sqrt((Real)m_number_variables) / 8.);
			computeRotation(m_rot, m_number_variables);
			m_linearTF.resize(m_number_variables, m_number_variables);
			for (size_t i = 0; i < m_number_variables; ++i) {
				for (size_t j = 0; j < m_number_variables; ++j) {
					m_linearTF[i][j] = m_scales * m_rot[i][j];
				}
			}
		}
		else if (m_index_number == 20) {
			m_fun = &BBOB::f20_schwefel;
			m_condition_number = 10.;
		}
		else if (m_index_number == 21) {
			m_fun = &BBOB::f21_gallagher_gaussian101me_peaks;
			computeRotation(m_rot, m_number_variables);
			m_condition_number = 1000.;

			m_peak_values.resize(NHIGHPEAKS21);
			std::vector<Real> fitvalues = { (Real)1.1, (Real)9.1 };
			std::vector<size_t> rperm(std::max(m_number_variables, NHIGHPEAKS21 - 1));
			std::vector<Real> arrCondition(NHIGHPEAKS21);
			for (size_t i = 0; i < NHIGHPEAKS21 - 1; ++i)
				rperm[i] = i;
			qsort(rperm.data(), NHIGHPEAKS21 - 1, sizeof(int), bbob::compareDoubles);
			arrCondition[0] = sqrt(m_condition_number);
			m_peak_values[0] = 10;
			for (size_t i = 1; i < NHIGHPEAKS21; ++i) {
				arrCondition[i] = pow(m_condition_number, (Real)(rperm[i - 1]) / ((Real)(NHIGHPEAKS21 - 2)));
				m_peak_values[i] = (Real)(i - 1) / (Real)(NHIGHPEAKS21 - 2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
			}

			for (size_t i = 0; i < NHIGHPEAKS21; ++i) {
				m_arrScales.push_back(std::vector<Real>(m_number_variables));
			}
			for (size_t i = 0; i < NHIGHPEAKS21; ++i) {
				for (size_t j = 0; j < m_number_variables; ++j)
					rperm[j] = j;
				qsort(rperm.data(), m_number_variables, sizeof(int), bbob::compareDoubles);
				for (size_t j = 0; j < m_number_variables; ++j)	{
					m_arrScales[i][j] = pow(arrCondition[i], ((Real)rperm[j]) / ((Real)(m_number_variables - 1)) - 0.5);
				}
			}
		}
		else if (m_index_number == 22) {
			m_fun = &BBOB::f22_gallagher_gaussian21hi_peaks;
			computeRotation(m_rot, m_number_variables);
			m_condition_number = 1000.;

			m_peak_values.resize(NHIGHPEAKS22);
			std::vector<Real> fitvalues = { (Real)1.1, (Real)9.1 };
			std::vector<size_t> rperm(std::max(m_number_variables, NHIGHPEAKS22 - 1));
			std::vector<Real> arrCondition(NHIGHPEAKS22);
			for (size_t i = 0; i < NHIGHPEAKS22 - 1; ++i)
				rperm[i] = i;
			qsort(rperm.data(), NHIGHPEAKS22 - 1, sizeof(int), bbob::compareDoubles);
			arrCondition[0] = sqrt(m_condition_number);
			m_peak_values[0] = 10;
			for (size_t i = 1; i < NHIGHPEAKS22; ++i) {
				arrCondition[i] = pow(m_condition_number, (Real)(rperm[i - 1]) / ((Real)(NHIGHPEAKS22 - 2)));
				m_peak_values[i] = (Real)(i - 1) / (Real)(NHIGHPEAKS22 - 2) * (fitvalues[1] - fitvalues[0]) + fitvalues[0];
			}

			for (size_t i = 0; i < NHIGHPEAKS22; ++i) {
				m_arrScales.push_back(std::vector<Real>(m_number_variables));
			}
			for (size_t i = 0; i < NHIGHPEAKS22; ++i) {
				for (size_t j = 0; j < m_number_variables; ++j)
					rperm[j] = j;
				qsort(rperm.data(), m_number_variables, sizeof(int), bbob::compareDoubles);
				for (size_t j = 0; j < m_number_variables; ++j) {
					m_arrScales[i][j] = pow(arrCondition[i], ((Real)rperm[j]) / ((Real)(m_number_variables - 1)) - 0.5);
				}
			}
		}
		else if (m_index_number == 23) {
			m_fun = &BBOB::f23_katsuura;
			m_condition_number = 100.;
			loadRotation(sqrt(m_condition_number));
		}
		else if (m_index_number == 24) {
			m_fun = &BBOB::f24_lunacekbi_rastrigin;
			m_condition_number = 100.;
			loadRotation(sqrt(m_condition_number));
			m_mu = 2.5;
		}
	}

	void BBOB::updateOptima(Environment *env) {
		auto temp_var = computeXopt();
		if (m_index_number == 4) {
			for (size_t i = 0; i < m_number_variables; i += 2) {
				temp_var[i] = fabs(temp_var[i]); /*Skew*/
			}
		}
		else if (m_index_number == 5) {
			for (size_t i = 0; i < m_number_variables; ++i) {
				if (temp_var[i] > 0) {
					temp_var[i] = 5.;
				}
				else if (temp_var[i] < 0) {
					temp_var[i] = -5.;
				}
			}
		}
		else if (m_index_number == 8) {
			for (size_t i = 0; i < m_number_variables; ++i)
				temp_var[i] *= 0.75;
		}
		else if (m_index_number == 9) {
			for (size_t i = 0; i < m_number_variables; ++i) {
				temp_var[i] = 0.;
				for (size_t j = 0; j < m_number_variables; ++j) {
					temp_var[i] += m_linearTF[j][i] * 0.5 / m_scales / m_scales;
				}
			}
		}
		else if (m_index_number == 19) {
			for (size_t i = 0; i < m_number_variables; ++i) {
				temp_var[i] = 0.;
				for (size_t j = 0; j < m_number_variables; ++j) {
					temp_var[i] += m_linearTF[j][i] * 0.5 / m_scales / m_scales;
				}
			}
		}
		else if (m_index_number == 20) {
			std::vector<Real> tmpvect(m_number_variables);
			for (auto &i : tmpvect) i = m_random->uniform.next();
			for (size_t i = 0; i < m_number_variables; ++i) {
				temp_var[i] = 0.5 * 4.2096874633;
				if (tmpvect[i] - 0.5 < 0)
					temp_var[i] *= -1.;
			}
		}
		else if (m_index_number == 21) {
			m_Xlocal.resize(m_number_variables, std::vector<Real>(NHIGHPEAKS21));
			std::vector<Real> peaks(m_number_variables * NHIGHPEAKS21);
			for (auto &j : peaks)
				j = m_random->uniform.next();
			for (size_t i = 0; i < m_number_variables; ++i) {
				temp_var[i] = 0.8 * (10. * peaks[i] - 5.);
				for (size_t j = 0; j < NHIGHPEAKS21; ++j) {
					m_Xlocal[i][j] = 0.;
					for (size_t k = 0; k < m_number_variables; ++k) {
						m_Xlocal[i][j] += m_rot[i][k] * (10. * peaks[j * m_number_variables + k] - 5.);
					}
					if (j == 0)
						m_Xlocal[i][j] *= 0.8;
				}
			}
		}
		else if (m_index_number == 22) {
			m_Xlocal.resize(m_number_variables, std::vector<Real>(NHIGHPEAKS22));
			std::vector<Real> peaks(m_number_variables * NHIGHPEAKS22);
			for (auto &j : peaks) 
				j = m_random->uniform.next();
			for (size_t i = 0; i < m_number_variables; ++i) {
				temp_var[i] = 0.8 * (10. * peaks[i] - 5.);
				for (size_t j = 0; j < NHIGHPEAKS22; ++j) {
					m_Xlocal[i][j] = 0.;
					for (size_t k = 0; k < m_number_variables; ++k) {
						m_Xlocal[i][j] += m_rot[i][k] * (10. * peaks[j * m_number_variables + k] - 5.);
					}
					if (j == 0)
						m_Xlocal[i][j] *= 0.8;
				}
			}
		}
		else if (m_index_number == 24) {
			std::vector<Real> tmpvect(m_number_variables);
			for (size_t i = 0; i < m_number_variables; ++i) {
				tmpvect[i] = m_random->uniform.next();
			}
			for (size_t i = 0; i < m_number_variables; ++i) {
				temp_var[i] = 0.5 * m_mu;
				if (tmpvect[i] < 0.)
					temp_var[i] *= -1.;
			}
		}
		m_optima.reset(new Optima<>());
		Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
		temp_sol.objective(0) = m_bias;
		temp_sol.variable().vector() = temp_var;
		dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);
		m_objective_accuracy = 1e-8;
	}

	std::vector<Real> BBOB::computeXopt() {
		std::vector<Real> temp(m_number_variables);
		for (size_t i = 0; i < m_number_variables; ++i) {
			temp[i] = 8 * floor(1e4 * m_random->uniform.next()) / 1e4 - 4;
			if (temp[i] == 0.0) {
				temp[i] = -1e-5;
			}
		}
		return temp;
	}

	Real BBOB::computeFopt() {
		Real gval, gval2;
		gval = m_random->normal.next();
		gval2 = m_random->normal.next();
		return fmin(1000., fmax(-1000., (round(100.*100.*gval / gval2) / 100.)));
	}

	Real BBOB::penalize(Real *x) {
		Real tmp, fPen = 0;
		for (size_t i = 0; i < m_number_variables; ++i) {
			tmp = fabs(x[i]) - 5.;
			if (tmp > 0.)
				fPen += tmp * tmp;
		}
		return fPen;
	}

	void BBOB::reshape(Matrix &B, std::vector<Real>& vector, size_t m, size_t n) {
		B.resize(m, n);
		for (size_t i = 0; i < m; ++i) {
			for (size_t j = 0; j < n; ++j) {
				B[i][j] = vector[j * m + i];
			}
		}
	}

	void BBOB::computeRotation(Matrix& rot, size_t Dim) {
		m_norRand.resize(Dim * Dim);
		for (auto &i : m_norRand) i = m_random->normal.next();
		reshape(rot, m_norRand, Dim, Dim);
		/*1st coordinate is row, 2nd is column.*/
		size_t i, j, k;
		Real prod;
		for (i = 0; i < Dim; ++i)
		{
			for (j = 0; j < i; ++j)
			{
				prod = 0;
				for (k = 0; k < Dim; ++k)
				{
					prod += rot[k][i] * rot[k][j];
				}
				for (k = 0; k < Dim; ++k)
				{
					rot[k][i] = rot[k][i] - prod * rot[k][j];
				}
			}
			prod = 0;
			for (k = 0; k < Dim; ++k)
			{
				prod += rot[k][i] * rot[k][i];
			}
			for (k = 0; k < Dim; ++k)
			{
				//rot[k][i] /= sqrt(prod);
				rot[k][i] = rot[k][i] / sqrt(prod);
			}
		}
	}
	bool BBOB::loadRotation(Real base_) {
		computeRotation(m_rot, m_number_variables);
		computeRotation(m_rot2, m_number_variables);
		/* decouple scaling from function definition*/
		size_t i, j, k;
		m_linearTF.resize(m_number_variables, m_number_variables);
		for (i = 0; i < m_number_variables; ++i) {
			for (j = 0; j < m_number_variables; ++j) {
				m_linearTF[i][j] = 0.;
				for (k = 0; k < m_number_variables; ++k) {
					m_linearTF[i][j] += m_rot[i][k] * pow(base_, ((Real)k) / ((Real)(m_number_variables - 1))) * m_rot2[k][j];
				}
			}
		}
		return true;
	}


	void BBOB::irregularize(Real *x) {
		// this method from BBOB
		Real c1, c2, x_;
		for (size_t i = 0; i < m_number_variables; ++i) {
			if (x[i]>0) {
				c1 = 10;	c2 = 7.9;
			}
			else {
				c1 = 5.5;	c2 = 3.1;
			}
			if (x[i] != 0) {
				x_ = log(fabs(x[i]));
			}
			else x_ = 0;
			x[i] = sign(x[i])*exp(x_ + 0.049*(sin(c1*x_) + sin(c2*x_)));
		}
	}

	void BBOB::asyemmetricalize(Real *x, Real belta) {
		// this method from BBOB
		if (m_number_variables == 1) return;
		for (size_t i = 0; i < m_number_variables; ++i) {
			if (x[i]>0) {
				x[i] = pow(x[i], 1 + belta*i*sqrt(x[i]) / (m_number_variables - 1));
			}
		}
	}



	void BBOB::evaluateObjective(Real* x, std::vector<Real>& obj) {
		(this->*m_fun)(x, obj);
	}

	void BBOB::f1_sphere(Real *x, std::vector<Real>& obj) {

		obj[0] = 0;
		for (size_t i = 0; i < m_number_variables; ++i) {
			x[i] -= optima()->solution(0).variable()[i];
		}
		for (size_t i = 0; i < m_number_variables; ++i) {
			obj[0] += x[i] * x[i];
		}
		obj[0] += m_bias;
	}

	void BBOB::f2_ellipsoidale(Real *x, std::vector<Real>& obj) {

		obj[0] = 0;
		for (size_t i = 0; i < m_number_variables; ++i) {
			x[i] -= optima()->solution(0).variable()[i];
		}
		irregularize(x);

		for (size_t i = 0; i < m_number_variables; ++i) {
			obj[0] += pow(m_condition_number, ((Real)i) / ((Real)(m_number_variables - 1))) * x[i] * x[i];
		}
		obj[0] += m_bias;
	}

	void BBOB::f3_rastrigin(Real *x, std::vector<Real>& obj) {

		for (size_t i = 0; i < m_number_variables; ++i) {
			x[i] -= optima()->solution(0).variable()[i];
		}
		irregularize(x);
		asyemmetricalize(x, m_beta);
		Real tmp, tmp2;
		for (size_t i = 0; i < m_number_variables; ++i) {
			tmp = ((Real)i) / ((Real)(m_number_variables - 1));
			x[i] = pow(sqrt(m_condition_number), tmp) * x[i];
		}
		/* COMPUTATION core*/
		tmp = 0.;
		tmp2 = 0.;
		for (size_t i = 0; i < m_number_variables; ++i) {
			tmp += cos(2 * OFEC_PI*x[i]);
			tmp2 += x[i] * x[i];
		}
		obj[0] = 10 * (m_number_variables - tmp) + tmp2;
		obj[0] += m_bias;
	}

	void BBOB::f4_buche_rastrigin(Real *x, std::vector<Real>& obj) {

		Real fPen = penalize(x);
		fPen *= 1e2;

		Real fAdd = optima()->solution(0).objective()[0] + fPen;

		for (size_t i = 0; i < m_number_variables; ++i) {
			x[i] -= optima()->solution(0).variable()[i];
		}

		irregularize(x);
		for (size_t i = 0; i < m_number_variables; ++i)
		{
			if (i % 2 == 0 && x[i] > 0)
				x[i] = sqrt(m_alpha) * x[i];
			x[i] = pow(sqrt(m_condition_number), ((Real)i) / ((Real)(m_number_variables - 1))) * x[i];
		}
		/* COMPUTATION core*/
		Real tmp = 0.;
		Real tmp2 = 0.;
		for (size_t i = 0; i < m_number_variables; ++i)
		{
			tmp += cos(2 * OFEC_PI*x[i]);
			tmp2 += x[i] * x[i];
		}
		obj[0] = 10 * (m_number_variables - tmp) + tmp2;
		obj[0] += fAdd;
	}

	void BBOB::f5_slope(Real *x, std::vector<Real>& obj) {

		/* BOUNDARY HANDLING*/
		/* move "too" good coordinates back into domain*/
		for (size_t i = 0; i < m_number_variables; ++i) {
			if ((optima()->solution(0).variable()[i] == 5.) && (x[i] > 5))
				x[i] = 5.;
			else if ((optima()->solution(0).variable()[i] == -5.) && (x[i] < -5))
				x[i] = -5.;
		}

		obj[0] = 0;
		/* COMPUTATION core*/
		for (size_t i = 0; i < m_number_variables; ++i)
		{
			if (optima()->solution(0).variable()[i] > 0) {
				obj[0] -= pow(sqrt(m_alpha), ((Real)i) / ((Real)(m_number_variables - 1))) * x[i];
			}
			else {
				obj[0] += pow(sqrt(m_alpha), ((Real)i) / ((Real)(m_number_variables - 1))) * x[i];
			}
		}
		obj[0] += m_bias;

	}

	void BBOB::f6_sector(Real *x, std::vector<Real>& obj) {
		/* attractive sector function*/
		size_t i, j; /*Loop over dim*/

					 /* BOUNDARY HANDLING*/
					 /* TRANSFORMATION IN SEARCH SPACE*/
		Real tmp;
		std::vector<Real> trasX(x, x + m_number_variables);
		for (i = 0; i < m_number_variables; ++i) {
			tmp = 0.;
			for (j = 0; j < m_number_variables; ++j) {
				tmp += m_linearTF[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}
			trasX[i] = tmp;
		}

		obj[0] = 0.;
		/* COMPUTATION core*/
		for (i = 0; i < m_number_variables; ++i)
		{
			if (trasX[i] * optima()->solution(0).variable()[i] > 0)
				trasX[i] *= m_alpha;
			obj[0] += trasX[i] * trasX[i];
		}

		/*MonotoneTFosc...*/
		if (obj[0] > 0)
		{
			obj[0] = pow(exp(log(obj[0]) / 0.1 + 0.49*(sin(log(obj[0]) / 0.1) + sin(0.79*log(obj[0]) / 0.1))), 0.1);
		}
		else if (obj[0] < 0)
		{
			obj[0] = -pow(exp(log(-obj[0]) / 0.1 + 0.49*(sin(0.55 * log(-obj[0]) / 0.1) + sin(0.31*log(-obj[0]) / 0.1))), 0.1);
		}
		obj[0] = pow(obj[0], 0.9);
		obj[0] += m_bias;

	}

	void BBOB::f7_step_ellipsoid(Real *x, std::vector<Real>& obj) {

		size_t i, j;
		Real x1, tmp;

		/* BOUNDARY HANDLING*/
		Real fPen = penalize(x);

		Real fAdd = optima()->solution(0).objective()[0] + fPen;

		std::vector<Real> tmpVect(m_number_variables), trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i) {
			tmpVect[i] = 0.;
			tmp = sqrt(pow(m_condition_number / 10., ((Real)i) / ((Real)(m_number_variables - 1))));
			for (j = 0; j < m_number_variables; ++j) {
				tmpVect[i] += tmp * m_rot2[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}

		}
		x1 = tmpVect[0];

		for (i = 0; i < m_number_variables; ++i) {
			if (fabs(tmpVect[i]) > 0.5)
				tmpVect[i] = round(tmpVect[i]);
			else
				tmpVect[i] = round(m_alpha * tmpVect[i]) / m_alpha;
		}

		for (i = 0; i < m_number_variables; ++i) {
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_rot[i][j] * tmpVect[j];
			}
		}

		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables; ++i)
		{
			obj[0] += pow(m_condition_number, ((Real)i) / ((Real)(m_number_variables - 1))) * trasX[i] * trasX[i];
		}
		obj[0] = 0.1 * fmax(1e-4 * fabs(x1), obj[0]);
		obj[0] += fAdd;
	}

	void BBOB::f8_original_rosenbrock(Real *x, std::vector<Real>& obj) {

		size_t i;
		Real tmp;

		/* TRANSFORMATION IN SEARCH SPACE*/
		std::vector<Real> trasX(m_number_variables);
		for (i = 0; i < m_number_variables; ++i) {
			trasX[i] = m_scales * (x[i] - optima()->solution(0).variable()[i]) + 1;
		}

		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables - 1; ++i)
		{
			tmp = (trasX[i] * trasX[i] - trasX[i + 1]);
			obj[0] += tmp * tmp;
		}
		obj[0] *= 1e2;
		for (i = 0; i < m_number_variables - 1; ++i)
		{
			tmp = (trasX[i] - 1.);
			obj[0] += tmp * tmp;
		}
		obj[0] += optima()->solution(0).objective()[0];
	}

	void BBOB::f9_rotated_rosenbrock(Real *x, std::vector<Real>& obj) {

		size_t i, j;
		Real  tmp;

		std::vector<Real> trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i <m_number_variables; ++i) {
			trasX[i] = 0.5;
			for (j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_linearTF[i][j] * x[j];
			}
		}

		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables - 1; ++i)
		{
			tmp = (trasX[i] * trasX[i] - trasX[i + 1]);
			obj[0] += tmp * tmp;
		}
		obj[0] *= 1e2;
		for (i = 0; i < m_number_variables - 1; ++i)
		{
			tmp = (trasX[i] - 1.);
			obj[0] += tmp * tmp;
		}
		obj[0] += m_bias;
	}

	void BBOB::f10_nonseparable_ellipsoid(Real *x, std::vector<Real>& obj) {

		size_t i, j;
		std::vector<Real> trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_rot[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}
		}

		irregularize(trasX.data());
		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables; ++i)
		{
			obj[0] += pow(m_condition_number, ((Real)i) / ((Real)(m_number_variables - 1))) * trasX[i] * trasX[i];
		}
		obj[0] += m_bias;
	}

	void BBOB::f11_discus(Real *x, std::vector<Real>& obj) {

		size_t i, j;
		std::vector<Real> trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_rot[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}
		}

		irregularize(trasX.data());

		/* COMPUTATION core*/
		obj[0] = m_condition_number * trasX[0] * trasX[0];
		for (i = 1; i < m_number_variables; ++i)
		{
			obj[0] += trasX[i] * trasX[i];
		}
		obj[0] += m_bias;
	}

	void BBOB::f12_bent_cigar(Real *x, std::vector<Real>& obj) {

		size_t i, j;

		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmpvect[i] = 0.;
			for (j = 0; j < m_number_variables; ++j) {
				tmpvect[i] += m_rot[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}
			if (tmpvect[i] > 0)
			{
				tmpvect[i] = pow(tmpvect[i], 1 + m_beta * ((Real)i) / ((Real)(m_number_variables - 1)) * sqrt(tmpvect[i]));
			}
		}

		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_rot[i][j] * tmpvect[j];
			}
		}

		/* COMPUTATION core*/
		obj[0] = trasX[0] * trasX[0];
		for (i = 1; i < m_number_variables; ++i)
		{
			obj[0] += m_condition_number * trasX[i] * trasX[i];
		}
		obj[0] += m_bias;
	}


	void BBOB::f13_sharp_ridge(Real *x, std::vector<Real>& obj) {

		size_t i, j/*, k*/; // keep an eye on warnings 

							/* BOUNDARY HANDLING*/

		std::vector<Real> trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_linearTF[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}
		}

		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 1; i < m_number_variables; ++i)
		{
			obj[0] += trasX[i] * trasX[i];
		}
		obj[0] = m_alpha * sqrt(obj[0]);
		obj[0] += trasX[0] * trasX[0];
		obj[0] += m_bias;
	}

	void BBOB::f14_different_powers(Real *x, std::vector<Real>& obj) {

		size_t i, j;

		std::vector<Real> trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_rot[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}
		}

		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables; ++i)
		{
			obj[0] += pow(fabs(trasX[i]), 2. + m_alpha * ((Real)i) / ((Real)(m_number_variables - 1)));
		}
		obj[0] = sqrt(obj[0]);
		obj[0] += m_bias;
	}

	void BBOB::f15_nonseparable_rastrigin(Real *x, std::vector<Real>& obj) {

		size_t i, j;

		Real tmp = 0., tmp2 = 0.;
		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmpvect[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				tmpvect[i] += m_rot[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}
		}
		irregularize(tmpvect.data());
		asyemmetricalize(tmpvect.data(), m_beta);
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				trasX[i] += m_linearTF[i][j] * tmpvect[j];
			}
		}
		/* COMPUTATION core*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmp += cos(2. * OFEC_PI * trasX[i]);
			tmp2 += trasX[i] * trasX[i];
		}
		obj[0] = 10. * ((Real)m_number_variables - tmp) + tmp2;
		obj[0] += m_bias;
	}

	void BBOB::f16_weierstrass(Real *x, std::vector<Real>& obj) {
		size_t i, j;
		Real tmp;

		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables);
		/* BOUNDARY HANDLING*/
		Real fPen = penalize(x), Fopt = m_bias;
		Fopt += 10. / (Real)m_number_variables * fPen;

		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmpvect[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				tmpvect[i] += m_rot[i][j] * (x[j] - optima()->solution(0).variable()[j]);
			}
		}

		irregularize(tmpvect.data());
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				trasX[i] += m_linearTF[i][j] * tmpvect[j];
			}
		}
		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables; ++i)
		{
			tmp = 0.;
			for (j = 0; j < 12; ++j)
			{
				tmp += cos(2 * OFEC_PI * (trasX[i] + 0.5) * m_bK[j]) * m_aK[j];
			}
			obj[0] += tmp;
		}
		obj[0] = 10. * pow(obj[0] / (Real)m_number_variables - m_F0, 3.);
		obj[0] += Fopt;
	}

	void BBOB::f17_schaffers_F7(Real *x, std::vector<Real>& obj) {

		size_t i, j;

		Real tmp;

		/* BOUNDARY HANDLING*/
		Real fPen = penalize(x), fOpt = m_bias;
		fOpt += 10. * fPen;

		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmpvect[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
				tmpvect[i] += m_rot[i][j] * (x[j] - optima()->solution(0).variable()[j]);
		}
		asyemmetricalize(tmpvect.data(), m_beta);

		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			tmp = pow(sqrt(m_condition_number), ((Real)i) / ((Real)(m_number_variables - 1)));
			for (j = 0; j < m_number_variables; ++j)
			{
				trasX[i] += tmp * m_rot2[i][j] * tmpvect[j];
			}
		}

		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables - 1; ++i)
		{
			tmp = trasX[i] * trasX[i] + trasX[i + 1] * trasX[i + 1];
			obj[0] += pow(tmp, 0.25) * (pow(sin(50 * pow(tmp, 0.1)), 2.) + 1.);
		}
		obj[0] = pow(obj[0] / (Real)(m_number_variables - 1), 2.);
		obj[0] += fOpt;
	}

	void BBOB::f18_illconditioned_schaffers_F7(Real *x, std::vector<Real>& obj) {

		size_t i, j;

		Real tmp;

		/* BOUNDARY HANDLING*/
		Real fPen = penalize(x), fOpt = m_bias;
		fOpt += 10. * fPen;

		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmpvect[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
				tmpvect[i] += m_rot[i][j] * (x[j] - optima()->solution(0).variable()[j]);
		}
		asyemmetricalize(tmpvect.data(), m_beta);

		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			tmp = pow(sqrt(m_condition_number), ((Real)i) / ((Real)(m_number_variables - 1)));
			for (j = 0; j < m_number_variables; ++j)
			{
				trasX[i] += tmp * m_rot2[i][j] * tmpvect[j];
			}
		}

		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables - 1; ++i)
		{
			tmp = trasX[i] * trasX[i] + trasX[i + 1] * trasX[i + 1];
			obj[0] += pow(tmp, 0.25) * (pow(sin(50. * pow(tmp, 0.1)), 2.) + 1.);
		}
		obj[0] = pow(obj[0] / (Real)(m_number_variables - 1), 2.);
		obj[0] += fOpt;
	}

	void BBOB::f19_composite_griewank_rosenbrock(Real *x, std::vector<Real>& obj) {

		size_t i, j;

		Real F2, tmp = 0., tmp2;

		std::vector<Real> trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i) {
			trasX[i] = 0.5;
			for (j = 0; j < m_number_variables; ++j) {
				trasX[i] += m_linearTF[i][j] * x[j];
			}
		}
		/* COMPUTATION core*/
		for (i = 0; i < m_number_variables - 1; ++i)
		{
			tmp2 = trasX[i] * trasX[i] - trasX[i + 1];
			F2 = 100. * tmp2 * tmp2;
			tmp2 = 1 - trasX[i];
			F2 += tmp2 * tmp2;
			tmp += F2 / 4000. - cos(F2);
		}
		obj[0] = 10. + 10. * tmp / (Real)(m_number_variables - 1);
		obj[0] += m_bias;
	}

	void BBOB::f20_schwefel(Real *x, std::vector<Real>& obj) {

		size_t i;
		Real tmp;
		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables);

		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmpvect[i] = 2. * x[i];
			if (optima()->solution(0).variable()[i] < 0.)
				tmpvect[i] *= -1.;
		}

		trasX[0] = tmpvect[0];
		for (i = 1; i < m_number_variables; ++i)
		{
			trasX[i] = tmpvect[i] + 0.25 * (tmpvect[i - 1] - 2. * fabs(optima()->solution(0).variable()[i - 1]));
		}

		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] -= 2 * fabs(optima()->solution(0).variable()[i]);
			trasX[i] *= pow(sqrt(m_condition_number), ((Real)i) / ((Real)(m_number_variables - 1)));
			trasX[i] = 100. * (trasX[i] + 2 * fabs(optima()->solution(0).variable()[i]));
		}

		Real fPen = 0, fOpt = m_bias;
		/* BOUNDARY HANDLING*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmp = fabs(trasX[i]) - 500.;
			if (tmp > 0.)
			{
				fPen += tmp * tmp;
			}
		}

		fOpt += 0.01 * fPen;
		/* COMPUTATION core*/
		obj[0] = 0.;
		for (i = 0; i < m_number_variables; ++i)
		{
			obj[0] += trasX[i] * sin(sqrt(fabs(trasX[i])));
		}
		obj[0] = 0.01 * ((418.9828872724339) - obj[0] / (Real)m_number_variables);
		obj[0] += fOpt;
	}

	void BBOB::f21_gallagher_gaussian101me_peaks(Real *x, std::vector<Real>& obj) {
		size_t i, j;
		Real a = 0.1;
		Real tmp2, f = 0., tmp;
		Real fac = -0.5 / (Real)m_number_variables;

		Real fAdd = penalize(x);

		/* BOUNDARY HANDLING*/

		m_bias = 1e4 * fAdd + optima()->solution(0).objective()[0];

		std::vector<Real> trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				trasX[i] += m_rot[i][j] * x[j];
			}
		}

		/* COMPUTATION core*/
		for (i = 0; i <NHIGHPEAKS21; ++i)
		{
			tmp2 = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				tmp = (trasX[j] - m_Xlocal[j][i]);
				tmp2 += m_arrScales[i][j] * tmp * tmp;
			}
			tmp2 = m_peak_values[i] * exp(fac * tmp2);
			f = fmax(f, tmp2);
		}

		f = 10. - f;
		if (f > 0)
		{
			obj[0] = log(f) / a;
			obj[0] = pow(exp(obj[0] + 0.49*(sin(obj[0]) + sin(0.79*obj[0]))), a);
		}
		else if (f < 0)
		{
			obj[0] = log(-f) / a;
			obj[0] = -pow(exp(obj[0] + 0.49*(sin(0.55 * obj[0]) + sin(0.31*obj[0]))), a);
		}
		else
			obj[0] = f;

		obj[0] *= obj[0];
		obj[0] += m_bias;
	}

	void BBOB::f22_gallagher_gaussian21hi_peaks(Real *x, std::vector<Real>& obj) {

		size_t i, j;
		Real a = 0.1;
		Real tmp2, f = 0., tmp;
		Real fac = -0.5 / (Real)m_number_variables;

		/* BOUNDARY HANDLING*/
		Real fPen = penalize(x), fOpt = m_bias;
		fOpt += fPen;

		std::vector<Real> trasX(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				trasX[i] += m_rot[i][j] * x[j];
			}
		}

		/* COMPUTATION core*/
		for (i = 0; i < NHIGHPEAKS22; ++i)
		{
			tmp2 = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				tmp = (trasX[j] - m_Xlocal[j][i]);
				tmp2 += m_arrScales[i][j] * tmp * tmp;
			}
			tmp2 = m_peak_values[i] * exp(fac * tmp2);
			f = fmax(f, tmp2);
		}

		f = 10. - f;
		if (f > 0)
		{
			obj[0] = log(f) / a;
			obj[0] = pow(exp(obj[0] + 0.49*(sin(obj[0]) + sin(0.79*obj[0]))), a);
		}
		else if (f < 0)
		{
			obj[0] = log(-f) / a;
			obj[0] = -pow(exp(obj[0] + 0.49*(sin(0.55 * obj[0]) + sin(0.31*obj[0]))), a);
		}
		else
			obj[0] = f;

		obj[0] *= obj[0];
		obj[0] += fOpt;
	}

	void BBOB::f23_katsuura(Real *x, std::vector<Real>& obj) {

		size_t i, j;

		Real tmp, arr, prod = 1., tmp2;

		/* BOUNDARY HANDLING*/
		Real fPen = penalize(x), fAdd = m_bias;
		fAdd += fPen;

		std::vector<Real> trasX(m_number_variables), tmpvect(m_number_variables);
		/* TRANSFORMATION IN SEARCH SPACE*/
		/* write rotated difference std::vector into tmx*/
		for (j = 0; j < m_number_variables; ++j)  /* store difference std::vector*/
			tmpvect[j] = x[j] - optima()->solution(0).variable()[j];
		Real *ptmx, *plinTF, *ptmp;
		for (i = 0; i < m_number_variables; ++i) {
			trasX[i] = 0.;
			ptmx = &trasX[i];
			plinTF = m_linearTF[i].vect().data();
			ptmp = tmpvect.data();
			for (j = 0; j < m_number_variables; ++j) {
				*ptmx += *plinTF++ * *ptmp++;
			}
		}

		/* COMPUTATION core*/
		for (i = 0; i < m_number_variables; ++i)
		{
			tmp = 0.;
			for (j = 1; j < 33; ++j)
			{
				tmp2 = pow(2., (Real)j);
				arr = trasX[i] * tmp2;
				tmp += fabs(arr - round(arr)) / tmp2;
			}
			tmp = 1. + tmp * (Real)(i + 1);
			prod *= tmp;
		}
		obj[0] = 10. / (Real)m_number_variables / (Real)m_number_variables * (-1. + pow(prod, 10. / pow((Real)m_number_variables, 1.2)));
		obj[0] += fAdd;
	}

	void BBOB::f24_lunacekbi_rastrigin(Real *x, std::vector<Real>& obj) {
		size_t i, j;

		Real tmp, tmp2 = 0., tmp3 = 0., tmp4 = 0.;
		Real s = 1. - 0.5 / (sqrt((Real)(m_number_variables + 20)) - 4.1);
		Real d = 1.;
		Real mu2 = -sqrt((m_mu * m_mu - d) / s);


		/* BOUNDARY HANDLING*/
		Real fPen = penalize(x);
		Real fAdd = optima()->solution(0).objective()[0] + 1e4 * fPen;
		std::vector<Real> trasX(m_number_variables);

		/* TRANSFORMATION IN SEARCH SPACE*/
		for (i = 0; i < m_number_variables; ++i)
		{
			trasX[i] = 2. * x[i];
			if (optima()->solution(0).variable()[i] < 0.)
				trasX[i] *= -1.;
		}

		/* COMPUTATION core*/
		tmp = 0.;
		for (i = 0; i < m_number_variables; ++i)
		{
			tmp2 += (trasX[i] - m_mu) * (trasX[i] - m_mu);
			tmp3 += (trasX[i] - mu2) * (trasX[i] - mu2);
			tmp4 = 0.;
			for (j = 0; j < m_number_variables; ++j)
			{
				tmp4 += m_linearTF[i][j] * (trasX[j] - m_mu);
			}
			tmp += cos(2 * OFEC_PI * tmp4);
		}
		obj[0] = fmin(tmp2, d * (Real)m_number_variables + s * tmp3) + 10. * ((Real)m_number_variables - tmp);
		obj[0] += fAdd;
	}

	int bbob::compareDoubles(const void *a, const void *b) {
		Real temp = (*(Real *)a) - (*(Real *)b);
		if (temp > 0)
			return 1;
		else if (temp < 0)
			return -1;
		else
			return 0;
	}
}