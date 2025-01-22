/********* Begin Register Information **********
{
	"name": "MIDACO",
	"identifier": "MIDACO",
	"problem tags": [ "SOP", "MOP", "ConOP", "COP", "MMOP", "GOP" ]
}
*********** End Register Information **********/

#include "../../../../core/algorithm/algorithm.h"

namespace ofec {
	class MIDACO : public Algorithm {
		OFEC_CONCRETE_INSTANCE(MIDACO)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void run_(Environment *env) override;
	private:
		std::unique_ptr<SolutionBase> m_sol;
	};
}
