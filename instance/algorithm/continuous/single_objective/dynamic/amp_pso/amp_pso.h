/********* Begin Register Information **********
{
	"name": "AMP-PSO",
	"identifier": "AMP_PSO",
	"problem tags": [ "ConOP", "GOP", "DOP", "MMOP", "SOP" ]
}
*********** End Register Information **********/

#include "../metrics_dynamic.h"
#include "../../../../template/framework/amp/continuous/amp_cont_ind.h"
#include "../../../../template/framework/amp/continuous/amp_cont.h"
#include "../../../../../../core/algorithm/algorithm.h"
#include "../../../../template/classic/pso/particle.h"
#include "../../../../template/classic/pso/swarm.h"
#ifdef OFEC_DEMO
#include "../dynamic_pso.h"
#endif

namespace ofec {
	class AMP_PSO : public MetricsDynamicConOEA
#ifdef OFEC_DEMO
		, public DynamicPSO
#endif
	{
	protected:
		using PopulationType = Swarm<IndContAMP<Particle>>;
		std::unique_ptr<ContAMP<PopulationType>> m_multi_pop;
		Real m_weight, m_accelerator1, m_accelerator2;
		size_t m_pop_size;

		void initialize_() override;
		void run_() override;
		void initMultiPop();
		void informChange(int rf);
		void measureMultiPop(bool flag);
#ifdef OFEC_DEMO
		void updateBuffer();
		std::vector<bool> getPopHiberState() const override;
#endif

	public:
		void record() override;
	};
}
