/********* Begin Register Information **********
{
	"name": "GTOP-Messenger",
	"identifier": "GTOP_Messenger",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_REALWORLD_TRAJECTORY_MESSENGER_H
#define OFEC_REALWORLD_TRAJECTORY_MESSENGER_H

#include "../../../../core/problem/continuous/continuous.h"

namespace ofec {
	namespace trajectory {
		class Messenger : public Continuous {
			OFEC_CONCRETE_INSTANCE(Messenger)
		protected:

			virtual void initialize_(Environment* env) override;
			void evaluateObjective(Real *x, std::vector<Real> &obj) override;
			void addInputParameters() {}
		};
	}
	using GTOP_Messenger = trajectory::Messenger;
}

#endif // !OFEC_REALWORLD_TRAJECTORY_MESSENGER_H
