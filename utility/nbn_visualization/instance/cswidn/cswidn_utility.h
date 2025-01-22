#ifndef OFEC_NBN_UTILITY_CSWIDN_H
#define OFEC_NBN_UTILITY_CSWIDN_H


#include <vector>
#include <set>
#include "../../../../core/random/newran.h"
#include "../../../../core/problem/solution.h"
#include "../../../../core/environment/environment.h"




namespace ofec {
	namespace nbn {
		namespace cswidn {

			void generateSolRasiu(std::shared_ptr<ofec::SolutionBase>& initSol,
				const ofec::SolutionBase& centerSol,
				double sampleRadius,
				ofec::Environment* env,
				ofec::Random* rnd);



			void generateSolIntRasiu(std::shared_ptr<ofec::SolutionBase>& initSol,
				const ofec::SolutionBase& centerSol,
				int sampleRadius,
				ofec::Environment* env,
				ofec::Random* rnd);
			

			void generateSolutionsInts(
				const ofec::SolutionBase& centerSol,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				ofec::Environment* env,
				ofec::Random* rnd
			);
			void generateSolutionsConRandom(
				const ofec::SolutionBase& centerSol,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				ofec::Environment* env,
				ofec::Random* rnd);

			void generateSolutionsConRadius(
				const ofec::SolutionBase& centerSol,
				double sampleRadius,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				ofec::Environment* env,
				ofec::Random* rnd
			);


			

		}
	}
}

#endif