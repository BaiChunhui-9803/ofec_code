#include "call_python.h"
#include "Python.h"
#include <iostream>

namespace lon {
	bool PythonCaller::m_python_flag = false;

	std::mutex PythonCaller::m_finish_flag_mtx;
	bool PythonCaller::m_finish_flag = true;
	std::vector<PythonCaller::PythonCallerPar> PythonCaller::m_pars;

	bool PythonCaller::initPython(std::vector<std::string>& path) {
		m_finish_flag = true;
		Py_Initialize();//使用python之前，要调用Py_Initialize();这个函数进行初始化
		if (!Py_IsInitialized())
		{
			printf("初始化失败！");
			m_python_flag = false;
			return false;
		}
		else {
			PyRun_SimpleString("import sys");
			std::string path_prefix = "sys.path.append('";
			std::string path_subfix = "')";
			std::string cur_path;
			for (auto& it : path) {
				cur_path = path_prefix + it + path_subfix;

				PyRun_SimpleString(cur_path.c_str());
			}

			m_python_flag = true;

			//PyRun_SimpleString("sys.path.append('E:/Diao_Yiya/paper/fitness_landscape/compare_alg/cmake_python/')");//这一步很重要，修改Python路径
			//PyRun_SimpleString("sys.path.append('E:/Diao_Yiya/paper/fitness_landscape/compare_alg/cmake_python/BasinHopping/')");//这一步很重要，修改Python路径
			return true;
		}
	}
	void PythonCaller::finishPython() {
		if (m_python_flag) {
			Py_Finalize();
		}
	}

	void PythonCaller::callFunction(const std::string& pythonFileName,
		const std::string& pythonFunName,
		const std::vector<std::string>& pars) {
		PyObject* pModule = nullptr;
		PyObject* pFunc = nullptr;


		pModule = PyImport_ImportModule(pythonFileName.c_str());//这里是要调用的文件名  pythonFileName.py
		if (pModule == NULL)
		{
			std::cout << "没找到该Python文件:\t" << pythonFileName << std::endl;
		}
		else {


			pFunc = PyObject_GetAttrString(pModule, pythonFunName.c_str());//这里是要调用的函数名

			PyObject* pArgs = PyTuple_New(pars.size());
			for (int idx(0); idx < pars.size(); ++idx) {
				PyTuple_SetItem(pArgs, idx, Py_BuildValue("s", pars[idx].c_str()));
			}

			PyObject* pRet = PyObject_CallObject(pFunc, pArgs);//调用函数

		}


	}

}