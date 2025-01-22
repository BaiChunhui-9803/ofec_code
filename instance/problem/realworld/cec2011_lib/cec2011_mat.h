/********* Begin Register Information **********
[
	{ "name":"CEC2011_MAT_F01", "identifier":"CEC2011_MAT_F01", "tags":[ "continuous", "single-objective" ], "dependency on libraries": [ "matlab" ] },
	{ "name":"CEC2011_MAT_F02", "identifier":"CEC2011_MAT_F02", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F03", "identifier":"CEC2011_MAT_F03", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F04", "identifier":"CEC2011_MAT_F04", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F05", "identifier":"CEC2011_MAT_F05", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F06", "identifier":"CEC2011_MAT_F06", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F07", "identifier":"CEC2011_MAT_F07", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F08", "identifier":"CEC2011_MAT_F08", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F09", "identifier":"CEC2011_MAT_F09", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F10", "identifier":"CEC2011_MAT_F10", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F11", "identifier":"CEC2011_MAT_F11", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F12", "identifier":"CEC2011_MAT_F12", "tags":[ "continuous", "single-objective" ] },
	{ "name":"CEC2011_MAT_F13", "identifier":"CEC2011_MAT_F13", "tags":[ "continuous", "single-objective" ] }

]
*********** End Register Information **********/

/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao & Changhe Li 
* Email: diaoyiyacug@gmail.com & changhe.lw@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/

/*
*  to do
	//{ "name":"CEC2011_MAT_F11", "identifier":"CEC2011_MAT_F11", "tags":[ "continuous", "single-objective" ] },
	
*/

#ifndef OFEC_CEC2011_MAT_H
#define OFEC_CEC2011_MAT_H


#include "../../../../core/problem/continuous/continuous.h"
#include "../../../../core/global.h"
#include "../../../../utility/functional.h"
#include <string>
#include <algorithm>
#include "CEC2011_MATLAB.h"
//#pragma comment(lib, "CEC2011_MATLAB.lib")
//#pragma comment(lib, "mclmcrrt.lib")

namespace ofec {
#define GET_CEC2011_MAT(pro) dynamic_cast<CEC2011_MAT*>(pro)


	template <int TfunId>
	class CEC2011_MAT : public Continuous {
		OFEC_ABSTRACT_INSTANCE(CEC2011_MAT)
	protected:

		void callBenchFun(Real* x, std::vector<Real>& obj) {
			mwArray xIn(1, m_number_variables, mxDOUBLE_CLASS);
			for (size_t i = 0; i < m_number_variables; ++i) {
				xIn(1, i + 1) = x[i];
			}
			mxDouble fun_numInData = TfunId;
			mwArray fun_numIn(fun_numInData);
			
			mwArray yOut;

			bench_func(1, yOut, xIn, fun_numIn);

			
			int numElements = yOut.NumberOfElements();
			obj.resize(numElements);
			for (int i = 0; i < numElements; ++i) {
				obj[i] = yOut(i + 1, 1);
			}
		}
		void callCostFun(Real* x, std::vector<Real>& obj) {

			mwArray xIn(m_number_variables, 1, mxDOUBLE_CLASS);
			for (size_t i = 0; i < m_number_variables; ++i) {
				xIn(i + 1, 1) = x[i];
			}

			mwArray yOut;
			mwArray PENALTYOut;
			mwArray rate_dOut;
			cost_fn(3, yOut, PENALTYOut, rate_dOut, xIn);


			int numElements = yOut.NumberOfElements();
			obj.resize(numElements);
			for (int i = 0; i < numElements; ++i) {
				obj[i] = yOut(i + 1, 1);
			}
		}


		void callMY_FUNCTION11_5(Real* x, std::vector<Real>& obj) {
			mwArray xIn(1,m_number_variables, mxDOUBLE_CLASS);
			for (size_t i = 0; i < m_number_variables; ++i) {
				xIn(1, i+1) = x[i];
			}

			mwArray yOut;
			mwArray CountOut;

			MY_FUNCTION11_5(2, yOut, CountOut, xIn);
		}
		void callAntennafunccircular(Real* x, std::vector<Real>& obj) {


			mwArray xIn(m_number_variables, 1, mxDOUBLE_CLASS);
			for (size_t i = 0; i < m_number_variables; ++i) {
				xIn(i + 1, 1) = x[i];
			}


			mxDouble nullInData[] = { 50.0, 120.0 };
			mwArray nullIn(1, 2, mxDOUBLE_CLASS);
			nullIn.SetData(nullInData, 2);

			mxDouble phi_desiredInData = 180.0;
			mwArray phi_desiredIn(phi_desiredInData);
			mxDouble distanceInData = 0.5;
			mwArray distanceIn(distanceInData);
			mwArray yOut;
			mwArray sllreturnOut;
			mwArray bwfnOut;


			antennafunccircular(3, yOut, sllreturnOut, bwfnOut, xIn, nullIn, phi_desiredIn, distanceIn);

			int numElements = yOut.NumberOfElements();
			obj.resize(numElements);
			for (int i = 0; i < numElements; ++i) {
				obj[i] = yOut(i + 1, 1);
			}

		}
		

		void callFUNCTION12_6(Real* x, std::vector<Real>& obj) {

			mwArray xIn(m_number_variables, 1, mxDOUBLE_CLASS);
			for (size_t i = 0; i < m_number_variables; ++i) {
				xIn(i + 1, 1) = x[i];
			}

			mwArray yOut;
			mwArray CountOut;
			MY_FUNCTION12_6(2, yOut, CountOut, xIn);

			int numElements = yOut.NumberOfElements();
			obj.resize(numElements);
			for (int i = 0; i < numElements; ++i) {
				obj[i] = yOut(i + 1, 1);
			}
		}


		void callFUNCTION13(Real* x, std::vector<Real>& obj) {

			mwArray xIn(m_number_variables, 1, mxDOUBLE_CLASS);
			for (size_t i = 0; i < m_number_variables; ++i) {
				xIn(i + 1, 1) = x[i];
			}

			mwArray yOut;
			mwArray CountOut;
			MY_FUNCTION13_1(2, yOut, CountOut, xIn);

			int numElements = yOut.NumberOfElements();
			obj.resize(numElements);
			for (int i = 0; i < numElements; ++i) {
				obj[i] = yOut(i + 1, 1);
			}
		}

		void addInputParameters() {
			switch (TfunId)
			{
			case 2: case 5: case 6: {
				m_input_parameters.add("number of variables", new RangedSizeT(m_number_variables, 3, 444, 30));
				break;
			}
			default:
				break;
			}
			
		}
		void initialize_(Environment* env) override {
			Continuous::initialize_(env);
			m_number_objectives = 1;
			m_optimize_mode.resize(m_number_objectives);
			m_optimize_mode[0] = OptimizeMode::kMinimize;
			
			if (TfunId == 1) {
				m_number_variables = 6;
				m_domain.resize(m_number_variables);
				setDomain(-6.4, 6.35);
			}
			else if(TfunId==2|| TfunId == 5 || TfunId == 6){
				if (m_number_variables < 3 || m_number_variables % 3 != 0) {
					throw Exception("The number of varialbes must be a multiple of 3.");
				}
				resizeVariable(m_number_variables);
				m_domain.setRange(0, 4, 0);
				m_domain.setRange(0, 4, 1);
				m_domain.setRange(0, OFEC_PI, 2);
				for (size_t i = 3; i < m_number_variables; ++i) {
					m_domain.setRange(-4 - floor((i - 4) / 3.0) / 4, 4 + floor((i - 4) / 3.0)/4, i);
				}
			}
			else if (TfunId == 3) {
				m_number_variables = 1;
				m_domain.resize(m_number_variables);
				setDomain(0.6,0.9);
			}
			else if (TfunId == 4) {
				m_number_variables = 1;
				m_domain.resize(m_number_variables);
				setDomain(0, 5);
				m_domain[0].limited = false;
				
			}
			else if (TfunId == 7) {
				resizeVariable(20);
				for (int idx(0); idx < m_number_variables; ++idx) {
					m_domain.setRange(0, 2 * OFEC_PI, idx);
				}
			}
			else if (TfunId == 8) {
				resizeVariable(7);
				for (int idx(0); idx < m_number_variables; ++idx) {
					m_domain.setRange(0, 15,idx);
				}
			}
			else if (TfunId == 9) {

				double basemva = 100; 
				double accuracy = 0.0001; 
				double maxiter = 10;
				double  tc = 283.4 * 4.52 * 2;
				std::vector<std::vector<double>> ebeData;
				readEBEdata(ebeData);
				std::vector<double> Pg;
				std::vector<double> Pd;

				for (const auto& bus : ebeData) {
					Pg.push_back(bus[6] / 100.0);
					Pd.push_back(bus[5] / 100.0);
				}
				

				std::vector<int> g;
				std::vector<int> d;

				for (int i = 0; i < Pg.size(); ++i) {
					if (Pg[i] > 0) {
						g.push_back(i);
					}
					if (Pd[i] > 0) {
						d.push_back(i);
					}
				}

				std::vector<std::vector<double>> BT(g.size(), std::vector<double>(d.size(), 0.0));

				BT[0][3] = 5;
				BT[0][4] = 10;
				BT[0][5] = 5;
				BT[1][2] = 5;
				BT[2][20] = 2.5;
				BT[3][20] = 2.5;
				BT[3][15] = 15;
				BT[4][11] = 2.5;
				BT[5][7] = 2.5;

				for (auto& row : BT) {
					for (auto& value : row) {
						value /= 100.0;
					}
				}

				m_GD_max.resize(g.size(), std::vector<double>(d.size(), 0.0));

				for (int i = 0; i < g.size(); ++i) {
					for (int j = 0; j < d.size(); ++j) {
						double temp1 = Pg[g[i]] - BT[i][j];
						double temp2 = Pd[d[j]] - BT[i][j];
						m_GD_max[i][j] = std::min(temp1, temp2);
					}
				}


				m_number_variables = g.size() * d.size();
				resizeVariable(m_number_variables);
				for (int i = 0; i < g.size(); ++i) {
					for (int j = 0; j < d.size(); ++j) {
						m_domain.setRange(0, m_GD_max[i][j], (i* d.size() + j));
					}
				}
				
			}
			else if (TfunId == 10) {
				resizeVariable(12);
				for (size_t i = 0; i < 6; ++i)
					m_domain.setRange(0.2, 1, i);
				for (size_t i = 6; i < 12; ++i)
					m_domain.setRange(-180, 180, i);
			}

			else if (TfunId == 11) {
				resizeVariable(120);

				const int No_of_Load_Hours = 24;
				const int No_of_Units = 5;
				Eigen::Matrix<double, 5, 24> Lower_Limit = Eigen::Matrix<double, 5, 24>::Constant(0);
				Eigen::Matrix<double, 5, 24> Upper_Limit = Eigen::Matrix<double, 5, 24>::Constant(0);

				Eigen::Matrix<double, 5, 1> Pmin;
				Pmin << 10, 20, 30, 40, 50;

				Eigen::Matrix<double, 5, 1> Pmax;
				Pmax << 75, 125, 175, 250, 300;

				for (int i = 0; i < 24; ++i) {
					Lower_Limit.col(i) = Pmin;
					Upper_Limit.col(i) = Pmax;
				}


				for (int i = 0; i < No_of_Load_Hours; ++i) {
					for (int j = 0; j < No_of_Units; ++j) {
						m_domain.setRange(Lower_Limit(j, i), Upper_Limit(j, i), j * No_of_Load_Hours + i);
						//Input_Generations(i, j) = input_vector(j * No_of_Load_Hours + i);
					}
				}

			}
			else if (TfunId==12){
				resizeVariable(26);
				m_domain.setRange(1900, 2300, 0);
				m_domain.setRange(2.5, 4.05, 1);
				m_domain.setRange(0, 1, 2);
				m_domain.setRange(0, 1, 3);
				for (size_t i = 4; i < 9; ++i)
					m_domain.setRange(100, 500, i);
				m_domain.setRange(100, 600, 9);
				for (size_t i = 10; i < 16; ++i)
					m_domain.setRange(0.01, 0.99, i);
				m_domain.setRange(1.1, 6, 16);
				m_domain.setRange(1.1, 6, 17);
				for (size_t i = 18; i < 21; ++i)
					m_domain.setRange(1.05, 6, i);
				for (size_t i = 21; i < 26; ++i)
					m_domain.setRange(-OFEC_PI, OFEC_PI, i);
			}
			else if (TfunId == 13) {
				resizeVariable(22);
				m_domain.setRange(-1000, 0, 0);
				m_domain.setRange(3, 5, 1);
				m_domain.setRange(0, 1, 2);
				m_domain.setRange(0, 1, 3);
				m_domain.setRange(100, 400, 4);
				m_domain.setRange(100, 500, 5);
				m_domain.setRange(30, 300, 6);
				m_domain.setRange(400, 1600, 7);
				m_domain.setRange(800, 2200, 8);
				for (size_t i = 9; i < 14; ++i)
					m_domain.setRange(0.01, 0.9, i);
				for (size_t i = 14; i < 18; ++i)
					m_domain.setRange(1.05, 6, i);
				m_domain.setRange(1.15, 6.5, 16);
				m_domain.setRange(1.7, 291, 17);
				for (size_t i = 18; i < 22; ++i)
					m_domain.setRange(-OFEC_PI, OFEC_PI, i);
			}
			
		}

		void evaluateObjective(Real* x, std::vector<Real>& obj) override {

			if (TfunId >= 1 && TfunId <= 8) {
				callBenchFun(x, obj);
			}
			else if (TfunId == 9) {
				callCostFun(x, obj);
			}
			else if (TfunId == 10) {
				callAntennafunccircular(x, obj);
			}
			else if (TfunId == 11) {
				callMY_FUNCTION11_5(x, obj);
			}
			
		}

		void readEBEdata(std::vector<std::vector<double>>& ebeData) {
			std::ifstream inputFile(std::string(OFEC_DIR) + "/instance/problem/realworld/cec2011_lib/data/EBEinputfile.txt");

			std::vector<double> curValue(11);
			std::string line;
			while (std::getline(inputFile, line)) {
				std::istringstream iss(line);
				for (auto& it : curValue) {
					iss >> it;
				}
				ebeData.push_back(curValue);
			}
			inputFile.close();
		}
		


	protected:
		// for problem9
		
		std::vector<std::vector<double>> m_GD_max;
		
	};


	class CEC2011_MAT_F01 : public CEC2011_MAT<1> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F01) protected: void addInputParameters() { } };
	class CEC2011_MAT_F02 : public CEC2011_MAT<2> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F02) protected: void addInputParameters() { } };
	class CEC2011_MAT_F03 : public CEC2011_MAT<3> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F03) protected: void addInputParameters() { } };
	class CEC2011_MAT_F04 : public CEC2011_MAT<4> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F04) protected: void addInputParameters() { } };
	class CEC2011_MAT_F05 : public CEC2011_MAT<5> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F05) protected: void addInputParameters() { } };
	class CEC2011_MAT_F06 : public CEC2011_MAT<6> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F06) protected: void addInputParameters() { } };
	class CEC2011_MAT_F07 : public CEC2011_MAT<7> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F07) protected: void addInputParameters() { } };
	class CEC2011_MAT_F08 : public CEC2011_MAT<8> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F08) protected: void addInputParameters() { } };
	class CEC2011_MAT_F09 : public CEC2011_MAT<9> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F09) protected: void addInputParameters() { } };
	class CEC2011_MAT_F10 : public CEC2011_MAT<10> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F10) protected: void addInputParameters() { } };
	class CEC2011_MAT_F11 : public CEC2011_MAT<11> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F11) protected: void addInputParameters() { } };
	class CEC2011_MAT_F12 : public CEC2011_MAT<12> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F12) protected: void addInputParameters() { } };
	class CEC2011_MAT_F13 : public CEC2011_MAT<13> { OFEC_CONCRETE_INSTANCE(CEC2011_MAT_F13) protected: void addInputParameters() { } };
}
#endif // !OFEC_CEC2011_MAT_H
