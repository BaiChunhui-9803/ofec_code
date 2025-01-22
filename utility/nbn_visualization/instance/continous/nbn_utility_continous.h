#ifndef NBN_MODIFY_SOLUTIONS_CONTINOUS_H
#define NBN_MODIFY_SOLUTIONS_CONTINOUS_H


#include <vector>
#include <set>
#include "../../../../core/random/newran.h"
#include "../../../../core/problem/solution.h"
#include "../../../../core/environment/environment.h"
#include "../../../../core/problem/continuous/continuous.h"


namespace ofec {
	namespace nbn {
		namespace continous {

			void fitlerSameSolutions(std::vector<ofec::SolutionBase*>& sols,
				std::vector<ofec::SolutionBase*>& filterSols, std::vector<int>& mapSolId, ofec::Environment* env) ;


			void generateSolutionsConRadius(
				const ofec::SolutionBase& centerSol,
				double sampleRadius,
				std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
				ofec::Environment* env,
				ofec::Random* rnd
			);
			
			void genenrateProblem(
				const std::string& proname,
				const ofec::ParameterMap& params,
				std::shared_ptr<ofec::Environment>& env);

		}
	}
}




#endif