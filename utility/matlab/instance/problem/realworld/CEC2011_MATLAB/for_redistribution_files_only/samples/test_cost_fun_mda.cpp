/*=================================================================
 *
 * TEST_COST_FUN
 * Sample driver code that uses the generic interface and
 * MATLAB Data API to call a C++ shared library created using 
 * MATLAB Compiler SDK.
 * Refer to the MATLAB Compiler SDK documentation for more
 * information.
 *
 *=================================================================*/

// Include the header file required to use the generic
// interface for the C++ shared library generated by the
// MATLAB Compiler SDK.
#include "MatlabCppSharedLib.hpp"
#include <iostream>

namespace mc = matlab::cpplib;
namespace md = matlab::data;

std::shared_ptr<mc::MATLABApplication> setup()
{
	auto mode = mc::MATLABApplicationMode::IN_PROCESS;
	// Specify MATLAB startup options
	std::vector<std::u16string> options = {};
	std::shared_ptr<mc::MATLABApplication> matlabApplication = mc::initMATLABApplication(mode, options);
	return matlabApplication;
}

int mainFunc(std::shared_ptr<mc::MATLABApplication> app, const int argc, const char * argv[])
{	
	md::ArrayFactory factory;
	md::TypedArray<double> xIn = factory.createArray<double>({126, 1}, {1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
	try {
		// The path to the CTF (library archive file) passed to 
		// initMATLABLibrary or initMATLABLibraryAsync may be either absolute
		// or relative. If it is relative, the following will be prepended
		// to it, in turn, in order to find the CTF:
		// - the directory named by the environment variable 
		// CPPSHARED_BASE_CTF_PATH, if defined
		// - the working directory
		// - the directory where the executable is located
		// - on Mac, the directory three levels above the directory
		// where the executable is located
		
		// If the CTF is not in one of these locations, do one of the following:
		// - copy the CTF
		// - move the CTF
		// - change the working directory ("cd") to the location of the CTF
		// - set the environment variable to the location of the CTF
		// - edit the code to change the path
		auto lib = mc::initMATLABLibrary(app, u"CEC2011_MATLAB.ctf");
		std::vector<md::Array> inputs{xIn};
		auto result = lib->feval(u"cost_fn", 3, inputs);
	} catch (const std::exception & exc) {
		std::cerr << exc.what() << std::endl;
		return -1;
	}
	return 0;
}

// The main routine. On the Mac, the main thread runs the system code, and
// user code must be processed by a secondary thread. On other platforms, 
// the main thread runs both the system code and the user code.
int main(const int argc, const char * argv[]) 
{
	int ret = 0;
	try {
		auto matlabApplication = setup();
		ret = mc::runMain(mainFunc, std::move(matlabApplication), argc, argv);
		// Calling reset() on matlabApplication allows the user to control
		// when it is destroyed, which automatically cleans up its resources.
		// Here, the object would go out of scope and be destroyed at the end 
		// of the block anyway, even if reset() were not called.
		// Whether the matlabApplication object is explicitly or implicitly
		// destroyed, initMATLABApplication() cannot be called again within
		// the same process.
		matlabApplication.reset();
	} catch(const std::exception & exc) {
		std::cerr << exc.what() << std::endl;
		return -1;
	}
	return ret;
}

