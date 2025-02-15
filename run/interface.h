#ifndef OFEC_INTERFACE_H
#define OFEC_INTERFACE_H

#include <string>
#include <memory>
#include <set>
#include <map>
#include "../core/environment/environment.h"
#include "factory.h"

namespace ofec {
	// Factory class for creating Environment objects
	extern Factory<Environment> g_factory_environment;
	// Factory class for creating Problem objects
	extern Factory<Problem> g_factory_problem;
	// Factory class for creating Algorithm objects
	extern Factory<Algorithm> g_factory_algorithm;
	// Mapping table storing the relationship between environment problems and supported algorithms
	extern std::map<std::string, std::set<std::string>> g_algorithm_for_environment_problem;
	void registerInstance();

	/**
 	 * @brief_byBCH Check if a specific environment problem supports a given algorithm.
 	 *
 	 * @param environment_problem_name The name of the environment problem.
 	 * @param algorithm_name The name of the algorithm.
 	 * @return true If the environment problem supports the algorithm.
 	 * @return false If the environment problem does not support the algorithm.
 	 */
	bool checkValidation(
		const std::string &environment_problem_name, 
		const std::string &algorithm_name
	);

	/**
 	 * @brief_byBCH Generate an Environment object based on the environment problem name.
 	 *
 	 * @param environment_problem_name The name of the environment problem.
 	 * @return Environment* Pointer to the created Environment object.
 	 */
	Environment* generateEnvironmentByFactory(const std::string &environment_problem_name);

	/**
  	 * @brief_byBCH Generate a Problem object based on the environment problem name.
  	 *
  	 * @param environment_problem_name The name of the environment problem.
  	 * @return Problem* Pointer to the created Problem object.
  	 */
	Problem* generateProblemByFactory(const std::string &environment_problem_name);

	/**
  	 * @brief_byBCH Generate an Algorithm object based on the algorithm name.
  	 *
  	 * @param algorithm_name The name of the algorithm.
  	 * @return Algorithm* Pointer to the created Algorithm object.
  	 */
	Algorithm* generateAlgorithmByFactory(const std::string &algorithm_name);
}

#endif // !OFEC_INTERFACE_H
