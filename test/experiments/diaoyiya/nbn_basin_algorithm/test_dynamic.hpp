#ifndef OFEC_CUSTOM_METHOD_HPP
#define OFEC_CUSTOM_METHOD_HPP

#include "../core/global.h"
#include "interface.h"
#include "../instance/problem/continuous/single_objective/dynamic/shifting_sphere/evaluation_shifting_sphere.h"
#include "../instance/algorithm/continuous/single_objective/global/canonical_de/canonical_de.h"

namespace ofec {
	void run() {
		// the path to data
		ofec::g_working_directory = "E:/Diao_Yiya/code/OFEC_v2_data";

		/* generate instances by identifier */
		{
			std::shared_ptr<Environment> env(EvaluationDynamicEnvironment::create());
			env->inputParameters().at("changeFre")->setValue("500");
			env->recordInputParameters();
			env->initialize(0.5);
			;
			env->setProblem(ShiftingSphere::create());
			//env->setProblem(ShfitingSphere::create());
			env->problem()->inputParameters().at("number of variables")->setValue("2");
			env->problem()->recordInputParameters();
			env->initializeProblem();

			env->setAlgorithm(Canonical_DE::create());
			env->algorithm()->inputParameters().at("population size")->setValue("10");
			env->algorithm()->inputParameters().at("maximum evaluations")->setValue("10000");
			env->algorithm()->recordInputParameters();
			env->initializeAlgorithm(0.5);
			env->runAlgorithm();
		}

		/* generate instances by name */
		{
			registerInstance();

			std::string environment_problem_name = "Evaluation_shifting_sphere";

			std::shared_ptr<Environment> env(generateEnvironmentByFactory(environment_problem_name));
			env->inputParameters().at("changeFre")->setValue("500");
			env->recordInputParameters();
			env->initialize(0.5);

			env->setProblem(generateProblemByFactory(environment_problem_name));
			env->problem()->inputParameters().at("number of variables")->setValue("2");
			env->problem()->recordInputParameters();
			env->initializeProblem();





			std::string algorithm_name = "Canonical-DE";
			if (checkValidation(environment_problem_name, algorithm_name)) {
				env->setAlgorithm(generateAlgorithmByFactory(algorithm_name));
				env->algorithm()->inputParameters().at("population size")->setValue("10");
				env->algorithm()->inputParameters().at("maximum evaluations")->setValue("10000");
				env->algorithm()->recordInputParameters();
				env->initializeAlgorithm(0.5);
				env->runAlgorithm();
			}
		}
	}
}

#endif // !OFEC_CUSTOM_METHOD_HPP