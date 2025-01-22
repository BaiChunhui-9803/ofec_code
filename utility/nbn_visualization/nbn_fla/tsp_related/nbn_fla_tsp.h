#ifndef NBN_FLA_TSP_H
#define NBN_FLA_TSP_H


#include <vector>
#include <set>
#include "../../../../core/random/newran.h"
#include "../../../../core/problem/solution.h"


#include "../nbn_info.h"
#include "../nbn_fla_utility.h"


namespace ofec {
	namespace nbn {
		namespace tsp {



			struct GapInfo {
				double m_maxlength = 0;
				int m_numNodes = 0;
				double m_avgBasinSize = 0;
				double maxNBD = 0;
				double m_minDirectionSons = 0;
				
			};

			// for tsp gap
			//  get narrowGap without parents nodes;
			void getNarrowGap(const NBNinfo& nbnInfo, std::vector<int>& bridgeIds, double bestFitThreadhold, double basinSisethread);

			

			// max gap length
			//void getTSP_narrowGap_maxLength(const NBNinfo& nbnInfo, const std::vector<int>& bridgeIds, std::vec);
			void calculateNarrowGapInfo(const NBNinfo& nbnInfo, const std::vector<int>& bridgeIds, GapInfo& info);


			void TspSolutionsToNBNinfo(NBNinfo& info, const std::string& saveDir, const std::string& filename, ofec::Environment* env);
		

		}
	}
}


#endif 
