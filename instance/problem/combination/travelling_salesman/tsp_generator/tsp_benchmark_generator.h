/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: YiyaDiao
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// link: http://archive.dimacs.rutgers.edu/Challenges/TSP/download.html
// Created: 27 DEC 2022

#ifndef OFEC_TSP_BENCHMARK_GENERATOR_H
#define OFEC_TSP_BENCHMARK_GENERATOR_H

#include <string>
#include <vector>
#include <sstream>

namespace ofec {

	struct PortGen {

		const int MAXCOORD = 1000000;
		const int PRANDMAX = 1000000000;

		int a, b;
		std::vector<int> arr;

		void init() {
			arr.resize(55);
		}

		void sprand(int);
		int lprand(void);
		void generate(int n,int seed, std::ostream & out);
	};


	struct PortCGen {

		const int  MAXN = 1000000;
		const int  MAXCOORD = 1000000;
		const int  PRANDMAX = 1000000000;
		const int  CLUSTERFACTOR = 100;
		const int  SCALEFACTOR = 1.0;

		std::vector<std::vector<int>> center;// [MAXN + 1] [2] ;
		int a, b;
		std::vector<int> arr;

		void init() {
			arr.resize(55);
			center.resize(MAXN + 1);
			for (auto& it : center) {
				it.resize(2);
			}
		}

		void sprand(int);
		double lprand(void);
		double normal(void);

		void generate(int n, int seed, std::ostream& out);
	};

	


	//void portgen(int n, int seed, const std::string& tspFile);
	//void portcgen(int n, int seed, const std::string& tspFile);
	//void portmgen(int n, int seed, const std::string& tspFile);
}

#endif