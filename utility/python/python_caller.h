#ifndef PYTHON_CALLER_H
#define PYTHON_CALLER_H

#include <vector>
#include <string>


namespace ofec {
	struct PythonCaller {

		static void initPython();
		static void finishPython();
		static void initPythonFunction(const std::string& pythonFilePath);
		static void callFunction(
			const std::string& pythonFileName,
			const std::string& pythonFunName,
			const std::vector<std::string>& pars);
	};
}


#endif