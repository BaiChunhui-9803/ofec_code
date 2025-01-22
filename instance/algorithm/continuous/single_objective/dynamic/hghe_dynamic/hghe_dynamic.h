/********* Begin Register Information **********
{
	"name": "HGHE-DYNAMIC",
	"identifier": "HGHE_DYNAMIC",
	"problem tags": [ "MMOP", "ConOP", "SOP", "GOP", "DOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_HGHE_DYNAMIC_H
#define OFEC_HGHE_DYNAMIC_H

#include "../../../../template/framework/hghe/hghe.h"
#include "../../../../template/framework/hghe/dynamic/hghec_dynamic_adaptor.h"
#include "../../../../template/classic/de/population.h"
#include "../metrics_dynamic.h"

namespace ofec {
	class HGHE_DYNAMIC : public HGHE, public MetricsDynamicConOEA{
	protected:
		int m_dynamic_strategy;
		Real m_scaling_factor, m_crossover_rate;
		PopDE<> m_pop;
		std::vector<VariableVector<Real>> m_pre_best_vars;

		void initialize_() override;
		// virtual int exploitHills() override;
		virtual void exploreHills() override;
		virtual void reproduceInHill(size_t id_hill) override;
		virtual int classifyPopulation() override;
		void updatePreBestVars(const std::vector<const SolutionBase* >& sols);
		virtual void selectSurvivors(
			std::vector<const SolutionBase*>& candidates,
			std::vector<std::unique_ptr<SolutionBase>>& population) override;

	public:
		AdaptorDynamicHGHEC * adaptorDynamicHGHEC() const;
		virtual void run_() override;
		virtual int evaluation() override;
		void handleDynamicChange(int rf, int strategy);
		virtual void record() override;
	private:
		int m_last_ner, m_last_nei;
	};


}

#endif // !OFEC_HGHE_DYNAMIC_H