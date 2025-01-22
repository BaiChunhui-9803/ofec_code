#ifndef NBN_IMPROVE_ACCURACY_MULTITHREAD_H
#define NBN_IMPROVE_ACCURACY_MULTITHREAD_H


#include "../nbn_fla/nbn_info.h"


namespace ofec {
	namespace nbn {


		void udpateNBNinfoIterationMultiThread(
			NBNinfo& m_nbnInfo,
			std::vector<bool>& updateIds,
			int from, int to,
			const std::vector<int>& sortedIds,
			ofec::Environment* env, ofec::Random* rnd);


	}
}


#endif