/********* Begin Register Information **********
{
	"name": "MOP_GLT4",
	"identifier": "GLT4",
	"problem tags": [ "MOP", "ConOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_GLT4_H
#define OFEC_GLT4_H

#include "glt.h"

namespace ofec {
	class GLT4 : public GLT {
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
