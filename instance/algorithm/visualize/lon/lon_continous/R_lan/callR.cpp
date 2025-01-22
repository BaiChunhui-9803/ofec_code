#include "callR.h"
namespace ofec::lon {
	const std::string RCaller::m_Rscript_path = "D:/R-4.3.0/bin/Rscript";
	void RCaller::runR(const std::string& cmd_line) {
		//	system("D:/R-4.2.1/bin/Rscript E:/Diao_Yiya/paper/fitness_landscape/compare_alg/LONs-Numerical-main/Graph/Example.R");
		std::string cmd = m_Rscript_path + " " + cmd_line;
		system(cmd.c_str());
	}
}