/********* Begin Register Information **********
{
	"name": "MOP_GLT3",
	"identifier": "GLT3",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_GLT3_H
#define OFEC_GLT3_H

#include "glt.h"

namespace ofec {
	class GLT3 : public GLT {
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