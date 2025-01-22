#ifndef OFEC_GL_CALCULATOR_H
#define OFEC_GL_CALCULATOR_H


#include <vector>
#include "../../../../core/problem/problem.h"


namespace ofec {
	struct GL_calculator {
		std::vector<std::vector<ofec::Real>> mvv_pro_mat;

		void initialize(ofec::Problem* pro);

		void calWeight(std::vector<double>& weight, const std::vector<double>& obj, ofec::Problem* pro);
		void udpateProMat(
			const std::vector<std::vector<std::vector<int>>*>& edges,
			std::vector<double>& weight);
		
		void generateTSPsolution(
			std::vector<int>& tour, 
			ofec::Problem* pro, ofec::Random* rnd);
	



		void generateTSPsolution(
			std::vector<int>& tour,
			ofec::Problem* pro, ofec::Random* rnd,
			const std::vector<std::vector<int>>& centerSol, double maxRadius
			);
	};
}


#endif 