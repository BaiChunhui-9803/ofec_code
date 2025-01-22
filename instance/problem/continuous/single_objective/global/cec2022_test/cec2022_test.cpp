#include "cec2022_test.h"

#include <filesystem>
#include <fstream>
#include "../../../../../../core/global.h"

namespace ofec {


#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029






	void CEC2022TestFunction::addInputParameters() {
		m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 1, 1000, 2));
		m_input_parameters.add("generate by random", new Bool(m_randGen, false));
	}

	void CEC2022TestFunction::initialize_(Environment* env) {
		Continuous::initialize_(env);
		m_number_objectives = 1;
		resizeVariable(m_number_variables);
		resizeObjective(1);
		m_optimize_mode[0] = OptimizeMode::kMinimize;
		m_domain.resize(m_number_variables);
		for (int idx(0); idx < m_number_variables; ++idx) {
			m_domain.setRange(-100, 100, idx);
		}

		int cf_num = 12;
		int nx = m_number_variables;
		auto func_num = m_index_number;


		Real condition_number = 1.0;

	//	std::string dataDir = "";

		std::string dataDir = "";
		if (!m_randGen) {
			dataDir = g_working_directory + "/instance/problem/continuous/single_objective/global/cec2022/";

			if (!std::filesystem::exists(dataDir)) {
				throw Exception("CEC2022Function data dir not exist");
			}
		}
		


		{
			std::string filename = "M_" + std::to_string(func_num) + "_D" + std::to_string(nx) + ".txt";
			if (func_num < 9)
			{
				//M = (double*)malloc(nx * nx * sizeof(double));
				m_M.resize(nx * nx);
			}
			else {
				m_M.resize(cf_num * nx * nx);
			}

			if (std::filesystem::exists(dataDir + filename)) {
				readVector(dataDir + filename, m_M);
			}
			else {

				if (func_num < 9)
				{
					//M = (double*)malloc(nx * nx * sizeof(double));
					m_M.resize(nx * nx);

					Matrix m_rotation(nx, nx);
					//			Matrix m_rotation(nx, nx);
					for (int icf(0); icf < 1; ++icf) {
						m_rotation.generateRotationClassical(&(m_random->normal), condition_number);
						for (int ix(0); ix < nx; ++ix) {
							for (int iy(0); iy < nx; ++iy) {
								m_M[icf * nx * nx + ix * nx + iy] = m_rotation[iy][ix];
							}
						}
					}
				}
				else {
					Matrix m_rotation(nx, nx);
					for (int icf(0); icf < cf_num; ++icf) {
						m_rotation.generateRotationClassical(&(m_random->normal), condition_number);
						for (int ix(0); ix < nx; ++ix) {
							for (int iy(0); iy < nx; ++iy) {
								m_M[icf * nx * nx + ix * nx + iy] = m_rotation[iy][ix];
							}
						}
					}
				}

			}
		}

		{
			FILE* fpt;
			//char FileName[256];
			auto filepath = dataDir + "shift_data_" + std::to_string(func_num) + ".txt";
			//	sprintf(FileName, dataDir+ "input_data/shift_data_%d.txt", func_num);
			fpt = fopen(filepath.c_str(), "r");
			if (fpt == NULL)
			{
				if (func_num < 9)
				{
					m_OShift.resize(nx);
				}
				else
				{
					m_OShift.resize(nx * cf_num);

				}
				for (auto& it : m_OShift) {
					it = m_random->uniform.nextNonStd<double>(-80, 80);
				}
			}

			else {
				if (func_num < 9)
				{
					m_OShift.resize(nx);
					for (int i = 0; i < nx; i++)
					{
						fscanf(fpt, "%lf", &m_OShift[i]);
					}
				}
				else
				{
					//OShift=(double *)malloc(nx*sizeof(double));


					m_OShift.resize(nx * cf_num);
					//OShift = (double*)malloc(nx * cf_num * sizeof(double));

					for (int i = 0; i < cf_num - 1; i++)
					{
						for (int j = 0; j < nx; j++)
						{
							fscanf(fpt, "%lf", &m_OShift[i * nx + j]);

						}
						fscanf(fpt, "%*[^\n]%*c");
					}
					for (int j = 0; j < nx; j++)
					{
						fscanf(fpt, "%lf", &m_OShift[nx * (cf_num - 1) + j]);

					}

				}
				fclose(fpt);
			}
		}

		//{
		//	if (func_num < 9) {

		//		std::string filename = "shift_data_" + std::to_string(func_num) + ".txt";
		//		m_OShift.resize(nx);
		//		std::ifstream in(dataDir + filename);
		//		for (int i = 0; i < nx; i++)
		//		{
		//			in >> m_OShift[i];
		//			//fscanf(fpt, "%lf", &OShift[i]);
		//		}
		//		in.close();
		//	}
		//	else {
		//		m_OShift.resize(nx * cf_num);

		//		std::string filename = "shift_data_" + std::to_string(func_num) + ".txt";
		//		{
		//			std::ifstream in(dataDir + filename);
		//			for (int i = 0; i < cf_num - 1; i++)
		//			{
		//				for (int j = 0; j < nx; j++)
		//				{
		//					in >> m_OShift[i * nx + j];
		//					//	fscanf(fpt, "%lf", &OShift[i * nx + j]);

		//				}
		//				//fscanf(fpt, "%*[^\n]%*c");
		//			}
		//			for (int j = 0; j < nx; j++)
		//			{
		//				in >> m_OShift[nx * (cf_num - 1) + j];
		//				//fscanf(fpt, "%lf", &OShift[nx * (cf_num - 1) + j]);
		//			}
		//			in.close();
		//		}
		//	}



		//	//readVector(dataDir + filename, m_OShift);
		//}
		//
		if (func_num >= 6 && func_num <= 8) {
			m_SS.resize(nx);
			std::string filename = "shuffle_data_" + std::to_string(func_num) + "_D" + std::to_string(nx) + ".txt";
			if (std::filesystem::exists(dataDir + filename)) {
				readVector(dataDir + filename, m_SS);
			}
			else {
				for (int idx(0); idx < m_SS.size(); ++idx) {
					m_SS[idx] = idx + 1;
				}
				m_random->uniform.shuffle(m_SS.begin(), m_SS.end());
			}
		}


		switch (m_index_number)
		{
		case 1:
			m_offsetY = 300.0;
			break;
		case 2:
			m_offsetY = 400.0;
			break;
		case 3:
			m_offsetY = 600.0;
			break;
		case 4:
			m_offsetY = 800.0;
			break;
		case 5:
			m_offsetY = 900.0;
			break;
		case 6:
			m_offsetY = 1800.0;
			break;
		case 7:
			m_offsetY = 2000.0;
			break;
		case 8:
			m_offsetY = 2200.0;
			break;
		case 9:
			m_offsetY = 2300.0;
			break;
		case 10:
			m_offsetY = 2400.0;
			break;
		case 11:
			m_offsetY = 2600.0;
			break;
		case 12:
			m_offsetY = 2700.0;
			break;
		default:
			throw Exception("Error: There are only 10 test functions in CEC2022 single objective test suite!\n");
			//printf("\nError: There are only 10 test functions in this test suite!\n");
			//f[i] = 0.0;
			break;
		}



	}


	void CEC2022TestFunction::evaluateObjective(Real* x, std::vector<Real>& obj)const {
		int nx = m_number_variables;
		std::vector<double> m_z, m_y;
		m_z.resize(nx + 10, 0);
		std::fill(m_z.begin(), m_z.end(), 0);
		m_y.resize(nx + 10, 0);
		std::fill(m_y.begin(), m_y.end(), 0);


		auto OShift= m_OShift;
		auto SS = m_SS;
		auto M = m_M;
		//std::vector<int> m_SS;
		EvaluateFun evalfun;
		evalfun.z = m_z.data();
		evalfun.y = m_y.data();
		evalfun.SS = SS.data();

		evalfun.OShift = OShift.data();
		evalfun.M = M.data();


		switch (m_index_number)
		{
		case 1:
			evalfun.zakharov_func(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1, 1);
			//f[i] += 300.0;
			break;
		case 2:
			evalfun.rosenbrock_func(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1, 1);
			//f[i] += 400.0;
			break;
		case 3:
			evalfun.schaffer_F7_func(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1, 1);
		//	f[i] += 600.0;
			break;
		case 4:
			evalfun.step_rastrigin_func(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1, 1);
		//	f[i] += 800.0;
			break;
		case 5:
			evalfun.levy_func(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1, 1);
		//	f[i] += 900.0;
			break;
		case 6:
			evalfun.hf02(x, obj.data(), nx, evalfun.OShift, evalfun.M, evalfun.SS, 1, 1);
		//	f[i] += 1800.0;
			break;
		case 7:
			evalfun.hf10(x, obj.data(), nx, evalfun.OShift, evalfun.M, evalfun.SS, 1, 1);
//			f[i] += 2000.0;
			break;
		case 8:
			evalfun.hf06(x, obj.data(), nx, evalfun.OShift, evalfun.M, evalfun.SS, 1, 1);
		//	f[i] += 2200.0;
			break;
		case 9:
			evalfun.cf01(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1);
		//	f[i] += 2300.0;
			break;
		case 10:
			evalfun.cf02(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1);
		//	f[i] += 2400.0;
			break;
		case 11:
			evalfun.cf06(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1);
		//	f[i] += 2600.0;
			break;
		case 12:
			evalfun.cf07(x, obj.data(), nx, evalfun.OShift, evalfun.M, 1);
		//	f[i] += 2700.0;
			break;
		default:

			obj[0] = 0;
		//	printf("\nError: There are only 10 test functions in this test suite!\n");
		//	f[i] = 0.0;
			break;
		}


		obj[0] += m_offsetY;

	}


	void CEC2022TestFunction::updateOptima(Environment* env) {

		m_optima.reset(new Optima<>());
		int cf = m_OShift.size() / m_number_variables;
		std::vector<Real> temp(m_number_variables);
		for (int idx(0); idx < cf; ++idx) {
			for (int j = 0; j < m_number_variables; ++j) {
				temp[j] = m_OShift[idx * m_number_variables + j];
			}
			Solution<> temp_sol(m_number_objectives, m_number_constraints, m_number_variables);
			//temp_sol.objective(0) = m_bias;
			temp_sol.variable().vector() = temp;
			evaluateObjective(temp_sol.variable().data(), temp_sol.objective());

			dynamic_cast<Optima<>&>(*m_optima).appendSolution(temp_sol);

		}

		double minvalue = std::numeric_limits<double>::max();
		for (int idx(0); idx < optima()->numberSolutions(); ++idx) {
			minvalue = std::min(minvalue, m_optima->solutionBase(idx).objective()[0]);
		}

		//if (minvalue != m_offsetY) {
		//	throw ofec::Exception("error optima position");
		//}



	}




	void  CEC2022TestFunction::EvaluateFun::ellips_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Ellipsoidal */
	{
		int i;
		f[0] = 0.0;
		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
		for (i = 0; i < nx; i++)
		{
			f[0] += pow(10.0, 6.0 * i / (nx - 1)) * z[i] * z[i];
		}
	}



	void  CEC2022TestFunction::EvaluateFun::bent_cigar_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Bent_Cigar */
	{
		int i;
		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
		f[0] = z[0] * z[0];
		for (i = 1; i < nx; i++)
		{
			f[0] += pow(10.0, 6.0) * z[i] * z[i];
		}

	}

	void  CEC2022TestFunction::EvaluateFun::discus_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Discus */
	{
		int i;
		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
		f[0] = pow(10.0, 6.0) * z[0] * z[0];
		for (i = 1; i < nx; i++)
		{
			f[0] += z[i] * z[i];
		}
	}



	void  CEC2022TestFunction::EvaluateFun::rosenbrock_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Rosenbrock's */
	{
		int i;
		double tmp1, tmp2;
		f[0] = 0.0;
		sr_func(x, z, nx, Os, Mr, 2.048 / 100.0, s_flag, r_flag); /* shift and rotate */
		z[0] += 1.0;//shift to orgin
		for (i = 0; i < nx - 1; i++)
		{
			z[i + 1] += 1.0;//shift to orgin
			tmp1 = z[i] * z[i] - z[i + 1];
			tmp2 = z[i] - 1.0;
			f[0] += 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
		}
	}



	void  CEC2022TestFunction::EvaluateFun::ackley_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Ackley's  */
	{
		int i;
		double sum1, sum2;
		sum1 = 0.0;
		sum2 = 0.0;

		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			sum1 += z[i] * z[i];
			sum2 += cos(2.0 * PI * z[i]);
		}
		sum1 = -0.2 * sqrt(sum1 / nx);
		sum2 /= nx;
		f[0] = E - 20.0 * exp(sum1) - exp(sum2) + 20.0;
	}




	void  CEC2022TestFunction::EvaluateFun::griewank_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Griewank's  */
	{
		int i;
		double s, p;
		s = 0.0;
		p = 1.0;

		sr_func(x, z, nx, Os, Mr, 600.0 / 100.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			s += z[i] * z[i];
			p *= cos(z[i] / sqrt(1.0 + i));
		}
		f[0] = 1.0 + s / 4000.0 - p;
	}

	void  CEC2022TestFunction::EvaluateFun::rastrigin_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Rastrigin's  */
	{
		int i;
		f[0] = 0.0;

		sr_func(x, z, nx, Os, Mr, 5.12 / 100.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * PI * z[i]) + 10.0);
		}
	}



	void  CEC2022TestFunction::EvaluateFun::schwefel_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Schwefel's  */
	{
		int i;
		double tmp;
		f[0] = 0.0;

		sr_func(x, z, nx, Os, Mr, 1000.0 / 100.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			z[i] += 4.209687462275036e+002;
			if (z[i] > 500)
			{
				f[0] -= (500.0 - fmod(z[i], 500)) * sin(pow(500.0 - fmod(z[i], 500), 0.5));
				tmp = (z[i] - 500.0) / 100;
				f[0] += tmp * tmp / nx;
			}
			else if (z[i] < -500)
			{
				f[0] -= (-500.0 + fmod(fabs(z[i]), 500)) * sin(pow(500.0 - fmod(fabs(z[i]), 500), 0.5));
				tmp = (z[i] + 500.0) / 100;
				f[0] += tmp * tmp / nx;
			}
			else
				f[0] -= z[i] * sin(pow(fabs(z[i]), 0.5));
		}
		f[0] += 4.189828872724338e+002 * nx;

	}



	void  CEC2022TestFunction::EvaluateFun::grie_rosen_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Griewank-Rosenbrock  */
	{
		int i;
		double temp, tmp1, tmp2;
		f[0] = 0.0;

		sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

		z[0] += 1.0;//shift to orgin
		for (i = 0; i < nx - 1; i++)
		{
			z[i + 1] += 1.0;//shift to orgin
			tmp1 = z[i] * z[i] - z[i + 1];
			tmp2 = z[i] - 1.0;
			temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;
			f[0] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
		}
		tmp1 = z[nx - 1] * z[nx - 1] - z[0];
		tmp2 = z[nx - 1] - 1.0;
		temp = 100.0 * tmp1 * tmp1 + tmp2 * tmp2;;
		f[0] += (temp * temp) / 4000.0 - cos(temp) + 1.0;
	}


	void  CEC2022TestFunction::EvaluateFun::escaffer6_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Expanded Scaffer��s F6  */
	{
		int i;
		double temp1, temp2;

		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

		f[0] = 0.0;
		for (i = 0; i < nx - 1; i++)
		{
			temp1 = sin(sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));
			temp1 = temp1 * temp1;
			temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]);
			f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
		}
		temp1 = sin(sqrt(z[nx - 1] * z[nx - 1] + z[0] * z[0]));
		temp1 = temp1 * temp1;
		temp2 = 1.0 + 0.001 * (z[nx - 1] * z[nx - 1] + z[0] * z[0]);
		f[0] += 0.5 + (temp1 - 0.5) / (temp2 * temp2);
	}

	void  CEC2022TestFunction::EvaluateFun::happycat_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* HappyCat, provdided by Hans-Georg Beyer (HGB) */
		/* original global optimum: [-1,-1,...,-1] */
	{
		int i;
		double alpha, r2, sum_z;
		alpha = 1.0 / 8.0;

		sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

		r2 = 0.0;
		sum_z = 0.0;
		for (i = 0; i < nx; i++)
		{
			z[i] = z[i] - 1.0;//shift to orgin
			r2 += z[i] * z[i];
			sum_z += z[i];
		}
		f[0] = pow(fabs(r2 - nx), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
	}

	void  CEC2022TestFunction::EvaluateFun::hgbat_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* HGBat, provdided by Hans-Georg Beyer (HGB)*/
		/* original global optimum: [-1,-1,...,-1] */
	{
		int i;
		double alpha, r2, sum_z;
		alpha = 1.0 / 4.0;

		sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

		r2 = 0.0;
		sum_z = 0.0;
		for (i = 0; i < nx; i++)
		{
			z[i] = z[i] - 1.0;//shift to orgin
			r2 += z[i] * z[i];
			sum_z += z[i];
		}
		f[0] = pow(fabs(pow(r2, 2.0) - pow(sum_z, 2.0)), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;
	}


	void  CEC2022TestFunction::EvaluateFun::schaffer_F7_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Schwefel's 1.2  */
	{
		int i;
		double tmp;
		f[0] = 0.0;
		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */
		for (i = 0; i < nx - 1; i++)
		{
			z[i] = pow(y[i] * y[i] + y[i + 1] * y[i + 1], 0.5);
			tmp = sin(50.0 * pow(z[i], 0.2));
			f[0] += pow(z[i], 0.5) + pow(z[i], 0.5) * tmp * tmp;
		}
		f[0] = f[0] * f[0] / (nx - 1) / (nx - 1);
	}

	void  CEC2022TestFunction::EvaluateFun::step_rastrigin_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Noncontinuous Rastrigin's  */
	{
		int i;
		f[0] = 0.0;
		for (i = 0; i < nx; i++)
		{
			if (fabs(y[i] - Os[i]) > 0.5)
				y[i] = Os[i] + floor(2 * (y[i] - Os[i]) + 0.5) / 2;
		}

		sr_func(x, z, nx, Os, Mr, 5.12 / 100.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			f[0] += (z[i] * z[i] - 10.0 * cos(2.0 * PI * z[i]) + 10.0);
		}
	}

	void  CEC2022TestFunction::EvaluateFun::levy_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Levy */
	{
		int i;
		f[0] = 0.0;
		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

		double* w;
		w = (double*)malloc(sizeof(double) * nx);

		double sum1 = 0.0;
		for (i = 0; i < nx; i++)
		{
			w[i] = 1.0 + (z[i] - 0.0) / 4.0;
		}

		double term1 = pow((sin(PI * w[0])), 2);
		double term3 = pow((w[nx - 1] - 1), 2) * (1 + pow((sin(2 * PI * w[nx - 1])), 2));

		double sum = 0.0;

		for (i = 0; i < nx - 1; i++)
		{
			double wi = w[i];
			double newv = pow((wi - 1), 2) * (1 + 10 * pow((sin(PI * wi + 1)), 2));
			sum = sum + newv;
		}

		f[0] = term1 + sum + term3;// - 1.442600987052770; // - 1.442600987052770
		free(w);   // ADD THIS LINE to free memory! Thanks for Dr. Janez
	}

	void  CEC2022TestFunction::EvaluateFun::zakharov_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* zakharov */
	{
		int i;
		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); // shift and rotate
		f[0] = 0.0;
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (i = 0; i < nx; i++)
		{
			double xi = z[i];
			sum1 = sum1 + pow(xi, 2);
			sum2 = sum2 + 0.5 * (i + 1) * xi;
		}

		f[0] = sum1 + pow(sum2, 2) + pow(sum2, 4);
	}

	void  CEC2022TestFunction::EvaluateFun::katsuura_func(double* x, double* f, int nx, double* Os, double* Mr, int s_flag, int r_flag) /* Katsuura  */
	{
		int i, j;
		double temp, tmp1, tmp2, tmp3;
		f[0] = 1.0;
		tmp3 = pow(1.0 * nx, 1.2);

		sr_func(x, z, nx, Os, Mr, 5.0 / 100.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			temp = 0.0;
			for (j = 1; j <= 32; j++)
			{
				tmp1 = pow(2.0, j);
				tmp2 = tmp1 * z[i];
				temp += fabs(tmp2 - floor(tmp2 + 0.5)) / tmp1;
			}
			f[0] *= pow(1.0 + (i + 1) * temp, 10.0 / tmp3);
		}
		tmp1 = 10.0 / nx / nx;
		f[0] = f[0] * tmp1 - tmp1;

	}


	void  CEC2022TestFunction::EvaluateFun::hf02(double* x, double* f, int nx, double* Os, double* Mr, int* S, int s_flag, int r_flag) /* Hybrid Function 2 */
	{
		int i, tmp, cf_num = 3;
		double fit[3];
		int G[3], G_nx[3];
		double Gp[3] = { 0.4,0.4,0.2 };

		tmp = 0;
		for (i = 0; i < cf_num - 1; i++)
		{
			G_nx[i] = ceil(Gp[i] * nx);
			tmp += G_nx[i];
		}
		G_nx[cf_num - 1] = nx - tmp;

		G[0] = 0;
		for (i = 1; i < cf_num; i++)
		{
			G[i] = G[i - 1] + G_nx[i - 1];
		}

		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			y[i] = z[S[i] - 1];
		}
		i = 0;
		bent_cigar_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 1;
		hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 2;
		rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

		f[0] = 0.0;
		for (i = 0; i < cf_num; i++)
		{
			f[0] += fit[i];
		}
	}

	void  CEC2022TestFunction::EvaluateFun::hf10(double* x, double* f, int nx, double* Os, double* Mr, int* S, int s_flag, int r_flag) /* Hybrid Function 6 */
	{
		int i, tmp, cf_num = 6;
		double fit[6];
		int G[6], G_nx[6];
		double Gp[6] = { 0.1,0.2,0.2,0.2,0.1,0.2 };

		tmp = 0;
		for (i = 0; i < cf_num - 1; i++)
		{
			G_nx[i] = ceil(Gp[i] * nx);
			tmp += G_nx[i];
		}
		G_nx[cf_num - 1] = nx - tmp;

		G[0] = 0;
		for (i = 1; i < cf_num; i++)
		{
			G[i] = G[i - 1] + G_nx[i - 1];
		}

		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			y[i] = z[S[i] - 1];
		}

		i = 0;
		hgbat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 1;
		katsuura_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 2;
		ackley_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 3;
		rastrigin_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 4;
		schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 5;
		schaffer_F7_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);

		f[0] = 0.0;
		for (i = 0; i < cf_num; i++)
		{
			f[0] += fit[i];
		}
	}

	void  CEC2022TestFunction::EvaluateFun::hf06(double* x, double* f, int nx, double* Os, double* Mr, int* S, int s_flag, int r_flag) /* Hybrid Function 6 */
	{
		int i, tmp, cf_num = 5;
		double fit[5];
		int G[5], G_nx[5];
		double Gp[5] = { 0.3,0.2,0.2,0.1,0.2 };

		tmp = 0;
		for (i = 0; i < cf_num - 1; i++)
		{
			G_nx[i] = ceil(Gp[i] * nx);
			tmp += G_nx[i];
		}
		G_nx[cf_num - 1] = nx - tmp;

		G[0] = 0;
		for (i = 1; i < cf_num; i++)
		{
			G[i] = G[i - 1] + G_nx[i - 1];
		}

		sr_func(x, z, nx, Os, Mr, 1.0, s_flag, r_flag); /* shift and rotate */

		for (i = 0; i < nx; i++)
		{
			y[i] = z[S[i] - 1];
		}
		i = 0;
		katsuura_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 1;
		happycat_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 2;
		grie_rosen_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 3;
		schwefel_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		i = 4;
		ackley_func(&y[G[i]], &fit[i], G_nx[i], Os, Mr, 0, 0);
		f[0] = 0.0;
		for (i = 0; i < cf_num; i++)
		{
			f[0] += fit[i];
		}
	}

	void  CEC2022TestFunction::EvaluateFun::cf01(double* x, double* f, int nx, double* Os, double* Mr, int r_flag) /* Composition Function 1 */
	{
		int i, cf_num = 5;
		double fit[5];
		double delta[5] = { 10, 20, 30, 40, 50 };
		double bias[5] = { 0, 200, 300, 100, 400 };

		i = 0;
		rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1e+4;
		i = 1;
		ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1e+10;
		i = 2;
		bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1e+30;
		i = 3;
		discus_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1e+10;
		i = 4;
		ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, 0);
		fit[i] = 10000 * fit[i] / 1e+10;
		cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
	}

	void  CEC2022TestFunction::EvaluateFun::cf02(double* x, double* f, int nx, double* Os, double* Mr, int r_flag) /* Composition Function 2 */
	{
		int i, cf_num = 3;
		double fit[3];
		double delta[3] = { 20,10,10 };
		double bias[3] = { 0, 200, 100 };

		i = 0;
		schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, 0);
		i = 1;
		rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		i = 2;
		hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
	}


	void  CEC2022TestFunction::EvaluateFun::cf06(double* x, double* f, int nx, double* Os, double* Mr, int r_flag) /* Composition Function 4 */
	{
		int i, cf_num = 5;
		double fit[5];
		double delta[5] = { 20,20,30,30,20 };
		double bias[5] = { 0, 200, 300, 400, 200 };
		i = 0;
		escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 2e+7;
		i = 1;
		schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		i = 2;
		griewank_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 1000 * fit[i] / 100;
		i = 3;
		rosenbrock_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		i = 4;
		rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1e+3;
		cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
	}

	void  CEC2022TestFunction::EvaluateFun::cf07(double* x, double* f, int nx, double* Os, double* Mr, int r_flag) /* Composition Function 4 */
	{
		int i, cf_num = 6;
		double fit[6];
		double delta[6] = { 10,20,30,40,50,60 };
		double bias[6] = { 0, 300, 500, 100, 400, 200 };
		i = 0;
		hgbat_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1000;
		i = 1;
		rastrigin_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1e+3;
		i = 2;
		schwefel_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 4e+3;
		i = 3;
		bent_cigar_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1e+30;
		i = 4;
		ellips_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 1e+10;
		i = 5;
		escaffer6_func(x, &fit[i], nx, &Os[i * nx], &Mr[i * nx * nx], 1, r_flag);
		fit[i] = 10000 * fit[i] / 2e+7;
		cf_cal(x, f, nx, Os, delta, bias, fit, cf_num);
	}




	void  CEC2022TestFunction::EvaluateFun::shiftfunc(double* x, double* xshift, int nx, double* Os)
	{
		int i;
		for (i = 0; i < nx; i++)
		{
			xshift[i] = x[i] - Os[i];
		}
	}

	void  CEC2022TestFunction::EvaluateFun::rotatefunc(double* x, double* xrot, int nx, double* Mr)
	{
		int i, j;
		for (i = 0; i < nx; i++)
		{
			xrot[i] = 0;
			for (j = 0; j < nx; j++)
			{
				xrot[i] = xrot[i] + x[j] * Mr[i * nx + j];
			}
		}
	}

	void  CEC2022TestFunction::EvaluateFun::sr_func(double* x, double* sr_x, int nx, double* Os, double* Mr, double sh_rate, int s_flag, int r_flag) /* shift and rotate */
	{
		int i;
		if (s_flag == 1)
		{
			if (r_flag == 1)
			{
				shiftfunc(x, y, nx, Os);
				for (i = 0; i < nx; i++)//shrink to the orginal search range
				{
					y[i] = y[i] * sh_rate;
				}
				rotatefunc(y, sr_x, nx, Mr);
			}
			else
			{
				shiftfunc(x, sr_x, nx, Os);
				for (i = 0; i < nx; i++)//shrink to the orginal search range
				{
					sr_x[i] = sr_x[i] * sh_rate;
				}
			}
		}
		else
		{

			if (r_flag == 1)
			{
				for (i = 0; i < nx; i++)//shrink to the orginal search range
				{
					y[i] = x[i] * sh_rate;
				}
				rotatefunc(y, sr_x, nx, Mr);
			}
			else
				for (i = 0; i < nx; i++)//shrink to the orginal search range
				{
					sr_x[i] = x[i] * sh_rate;
				}
		}
	}

	void  CEC2022TestFunction::EvaluateFun::asyfunc(double* x, double* xasy, int nx, double beta)
	{
		int i;
		for (i = 0; i < nx; i++)
		{
			if (x[i] > 0)
				xasy[i] = pow(x[i], 1.0 + beta * i / (nx - 1) * pow(x[i], 0.5));
		}
	}

	void  CEC2022TestFunction::EvaluateFun::oszfunc(double* x, double* xosz, int nx)
	{
		int i, sx;
		double c1, c2, xx;
		for (i = 0; i < nx; i++)
		{
			if (i == 0 || i == nx - 1)
			{
				if (x[i] != 0)
					xx = log(fabs(x[i]));
				if (x[i] > 0)
				{
					c1 = 10;
					c2 = 7.9;
				}
				else
				{
					c1 = 5.5;
					c2 = 3.1;
				}
				if (x[i] > 0)
					sx = 1;
				else if (x[i] == 0)
					sx = 0;
				else
					sx = -1;
				xosz[i] = sx * exp(xx + 0.049 * (sin(c1 * xx) + sin(c2 * xx)));
			}
			else
				xosz[i] = x[i];
		}
	}


	void  CEC2022TestFunction::EvaluateFun::cf_cal(double* x, double* f, int nx, double* Os, double* delta, double* bias, double* fit, int cf_num)
	{
		int i, j;
		double* w;
		double w_max = 0, w_sum = 0;
		w = (double*)malloc(cf_num * sizeof(double));
		for (i = 0; i < cf_num; i++)
		{
			fit[i] += bias[i];
			w[i] = 0;
			for (j = 0; j < nx; j++)
			{
				w[i] += pow(x[j] - Os[i * nx + j], 2.0);
			}
			if (w[i] != 0)
				w[i] = pow(1.0 / w[i], 0.5) * exp(-w[i] / 2.0 / nx / pow(delta[i], 2.0));
			else
				w[i] = INF;
			if (w[i] > w_max)
				w_max = w[i];
		}

		for (i = 0; i < cf_num; i++)
		{
			w_sum = w_sum + w[i];
		}
		if (w_max == 0)
		{
			for (i = 0; i < cf_num; i++)
				w[i] = 1;
			w_sum = cf_num;
		}
		f[0] = 0.0;
		for (i = 0; i < cf_num; i++)
		{
			f[0] = f[0] + w[i] / w_sum * fit[i];
		}
		free(w);
	}

}