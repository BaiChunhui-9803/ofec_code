/********* Begin Register Information **********
{
	"name": "Modified_uneven_de_maxima_1",
	"identifier": "ModifiedUnevenDeMaxima1",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_MODIFIED_UNEVEN_DE_MAXIMA_1_H
#define OFEC_MODIFIED_UNEVEN_DE_MAXIMA_1_H

#include "../../../init_pop_bounded.h"

namespace ofec {
	class ModifiedUnevenDeMaxima1 : public InitPopBounded {
		OFEC_CONCRETE_INSTANCE(ModifiedUnevenDeMaxima1)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}
#endif // !OFEC_MODIFIED_UNEVEN_DE_MAXIMA_1_H