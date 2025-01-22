#include "python_caller.h"
#ifdef  OFEC_PYTHON
#include "Python.h"
#endif //  OFEC_PYTHON
#include "../../core/exception.h"


namespace ofec {

	void PythonCaller::initPython() {
#ifdef  OFEC_PYTHON
		Py_Initialize();// before python£¬call Py_Initialize();
		if (!Py_IsInitialized())
		{
			throw ofec::Exception("Python can not be initialized");
		}
		PyRun_SimpleString("import sys");
#endif
	}
	void PythonCaller::finishPython() {
#ifdef  OFEC_PYTHON
		Py_Finalize();
#endif
	}

	void PythonCaller::initPythonFunction(const std::string& pythonFilePath) {
#ifdef  OFEC_PYTHON

		std::string path_prefix = "sys.path.append('";
		std::string path_subfix = "')";
		std::string cur_path;
		cur_path = path_prefix + pythonFilePath + path_subfix;
		PyRun_SimpleString(cur_path.c_str());
#endif
	}

	void PythonCaller::callFunction(const std::string& pythonFileName,
		const std::string& pythonFunName,
		const std::vector<std::string>& pars) {

#ifdef  OFEC_PYTHON
		PyObject* pModule = nullptr;
		PyObject* pFunc = nullptr;
		pModule = PyImport_ImportModule(pythonFileName.c_str());//  pythonFileName.py
		if (pModule == NULL)
		{
			std::string error_info = "Python file not found: " + pythonFileName;
			throw ofec::MyExcept(error_info);
		}
		else {
			pFunc = PyObject_GetAttrString(pModule, pythonFunName.c_str());//name of function
			PyObject* pArgs = PyTuple_New(pars.size());
			for (int idx(0); idx < pars.size(); ++idx) {
				PyTuple_SetItem(pArgs, idx, Py_BuildValue("s", pars[idx].c_str()));
			}
			PyObject* pRet = PyObject_CallObject(pFunc, pArgs);//call function
		}
#endif
	}

}