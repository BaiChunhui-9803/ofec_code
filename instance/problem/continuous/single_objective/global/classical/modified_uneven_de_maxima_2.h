/********* Begin Register Information **********
{
	"name": "Modified_uneven_de_maxima_2",
	"identifier": "ModifiedUnevenDeMaxima2",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_MODIFIED_UNEVEN_DE_MAXIMA_2_H
#define OFEC_MODIFIED_UNEVEN_DE_MAXIMA_2_H

#include "../../../init_pop_bounded.h"

namespace ofec {
	class ModifiedUnevenDeMaxima2 : public InitPopBounded {
		OFEC_CONCRETE_INSTANCE(ModifiedUnevenDeMaxima2)
	protected:
		void addInputParameters() {}
		void initialize_(Environment *env) override;
		void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
	};
}
#endif // !OFEC_MODIFIED_UNEVEN_DE_MAXIMA_2_H