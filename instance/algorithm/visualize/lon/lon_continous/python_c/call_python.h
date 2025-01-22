#ifndef LON_CALL_PYTHON_H
#define LON_CALL_PYTHON_H

#include <vector>
#include <string>
#include <mutex>


namespace lon {
	struct PythonCaller {


		struct PythonCallerPar {
			std::string m_pythonFileName;
			std::string m_pythonFunName;
			std::vector<std::string>  m_pars;
		};


		static bool m_python_flag;
		static std::mutex m_finish_flag_mtx;
		static bool m_finish_flag;
		static std::vector<PythonCallerPar> m_pars;





		static bool initPython(std::vector<std::string>& path);
		static void finishPython();

		static void callFunction(const std::string& pythonFileName,
			const std::string& pythonFunName,
			const std::vector<std::string>& pars);


	};
}


#endif