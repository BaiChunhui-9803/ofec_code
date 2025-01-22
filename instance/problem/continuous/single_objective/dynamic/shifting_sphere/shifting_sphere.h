#include "../../global/classical/sphere.h"
#include "../../../../../environment/template/uncertainty/dynamic_problem.h"

namespace ofec {
	class ShiftingSphere : public Sphere, public DynamicProblem {
		OFEC_CONCRETE_INSTANCE(ShiftingSphere)
	public:
		void addInputParameters() {}
		void change(Random* rnd) override;
	};
}
