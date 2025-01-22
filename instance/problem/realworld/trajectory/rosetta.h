/********* Begin Register Information **********
{
	"name": "GTOP-Rosetta",
	"identifier": "GTOP_Rosetta",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_REALWORLD_TRAJECTORY_ROSETTA_H
#define OFEC_REALWORLD_TRAJECTORY_ROSETTA_H

#include "../../../../core/problem/continuous/continuous.h"

namespace ofec {
	namespace trajectory {
		class Rosetta : public Continuous {
		protected:
			void initialize_() override;
			void evaluateObjective(Real *x, std::vector<Real> &obj) override;
		};
	}
	using GTOP_Rosetta = trajectory::Rosetta;
}

#endif // !OFEC_REALWORLD_TRAJECTORY_ROSETTA_H
