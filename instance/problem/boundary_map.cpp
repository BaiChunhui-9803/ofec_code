#include "boundary_map.h"
#include "../../core/global.h"
#include "../../core/parameter/parameter_variant.h"
#include "../../core/parameter/variants_to_stream.h"

namespace ofec {


	//void ProblemBoundaryBase::addInputParameters() {
	//	m_input_parameters.add("tabu setup", new Bool(m_tabuSetup, false));
	//}

	void ProblemBoundaryBase::inputFile() {
		if (m_tabuSetup) {
			std::string dir = g_working_directory + "/instance/problem/continuous/tabu_setup/";
			std::string filename = name();

		//	std::string absFilename = dir + name() + "_bestSolution.txt";

			m_radius.clear();
		    m_radius.push_back(0.01);
			m_tabuSolutions.clear();
			m_tabuSolutions.emplace_back(createSolution());

			{
				ofec::ParameterVariantStream paramsStream;
			//	std::cout << "input solutions" << std::endl;
				ofec::variants_stream::inputFromFile(paramsStream, dir + name() + "_bestSolution.txt");

				parameterVariantsToSolution(paramsStream ,*m_tabuSolutions.front());



				//ParameterVariantStream paramstream;
			}
		}
	}
}