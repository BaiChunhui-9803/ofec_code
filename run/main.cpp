#include <iostream>
#ifdef OFEC_UNIT_TEST
#include "../test/catch_amalgamated.hpp"
#else
#include "custom_method.hpp"
#endif
#ifdef OFEC_MATLAB
#include "../utility/matlab/matlab_define.h"
#endif // OFEC_MATLAB
#ifdef OFEC_PYTHON
#include "../utility/python/python_caller.h"
#endif // OFEC_PYTHON
#include "../test/example/example_1.h"
#include "../test/example/example_2.h"

#define OFEC_EXAMPLE_1

int main(int argc, char* argv[]) {

#ifdef OFEC_MATLAB
	ofec::initMatlabFunction();
#endif // OFEC_MATLAB

#ifdef OFEC_PYTHON
	ofec::PythonCaller::initPython();
#endif 
	time_t timer_start, timer_end;
	time(&timer_start);
#ifdef OFEC_UNIT_TEST
	int result = Catch::Session().run(argc, argv);
#elifdef OFEC_EXAMPLE_1
	ofec::run_example_1(argc, argv);
#elifdef OFEC_EXAMPLE_2
	ofec::run_example_2(argc, argv);
#else
	ofec::run(argc, argv);
#endif
	time(&timer_end);
	std::cout << "Time used: " << difftime(timer_end, timer_start) << " seconds" << std::endl;

#ifdef OFEC_MATLAB
	ofec::terminateMatlabFunction();
#endif // OFEC_MATLAB

#ifdef OFEC_PYTHON
	ofec::PythonCaller::finishPython();
#endif // OFEC_MATLAB
	return 0;
}

