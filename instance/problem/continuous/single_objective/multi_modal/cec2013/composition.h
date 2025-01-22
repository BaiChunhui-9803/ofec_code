#ifndef OFEC_CEC13_COMPOSITION_H
#define OFEC_CEC13_COMPOSITION_H

#include "../../global/cec2005/composition.h"

namespace ofec {
	namespace cec2013 {
		class Composition : public cec2005::Composition {
		public:
			void initialize_(Environment *env) override;
			void updateOptima(Environment *env) override;
		protected:
			void evaluateObjective(Real *x, std::vector<Real> &obj) const override;
		};
	}
}

#endif // !OFEC_CEC13_COMPOSITION_H