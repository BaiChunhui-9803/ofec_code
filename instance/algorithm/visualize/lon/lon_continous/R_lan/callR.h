#ifndef LON_CALL_R_H
#define LON_CALL_R_H

#include <string>

namespace ofec::lon {
	struct RCaller {

		const static std::string m_Rscript_path;
		static void runR(const std::string& cmd);

	};
}

//extern void CallR(const std::string& R_filePath,
//     
//	);


#endif 