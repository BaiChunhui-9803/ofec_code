/********* Begin Register Information **********
{
	"name": "AMP-DE",
	"identifier": "AMP_DE",
	"tags": [ "dynamic", "continuous", "single-objective" ]
}
*********** End Register Information **********/

#include "../../../../template/framework/amp/continuous/amp_cont_ind.h"
#include "../../../../template/framework/amp/continuous/amp_cont.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/differential_evolution/population.h"

namespace ofec {
	class AMP_DE : public Algorithm {
		OFEC_CONCRETE_INSTANCE(AMP_DE)
	protected:
		using PopulationType = PopulationDE<IndContAMP<IndividualDE>>;
		std::unique_ptr<ContAMP<PopulationType>> m_multi_pop;
		size_t m_pop_size;
		
		void addInputParamters();
		void run_(Environment *env) override;
	};
}