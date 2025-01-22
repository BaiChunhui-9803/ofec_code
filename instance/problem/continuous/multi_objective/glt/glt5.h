/********* Begin Register Information **********
{
	"name": "MOP_GLT5",
	"identifier": "GLT5",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_GLT5_H
#define OFEC_GLT5_H

#include "glt.h"

namespace ofec {
	class GLT5 : public GLT {
	public:
		void updateOptima() override;
		void sampleParetoSols(size_t sample_num);

		/*void sampleParetoFront(size_t sample_num);
		void loadParetoFront(size_t sample_num);*/
	protected:
		void evaluateObjective(Real* x, std::vector<Real>& obj) override;
	};
}

#endif
