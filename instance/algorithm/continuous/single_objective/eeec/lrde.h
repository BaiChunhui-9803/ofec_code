/********* Begin Register Information **********
{
	"name": "LRDE",
	"identifier": "LRDE",
	"tags": [ "continuous", "single-objective" ]
}
*********** End Register Information **********/

#ifndef OFEC_LOCAL_REPULSION_DE_H
#define OFEC_LOCAL_REPULSION_DE_H

#include "../../../../../core/algorithm/algorithm.h"
#include "../../../template/classic/differential_evolution/population.h"

namespace ofec {
	class LRDE : virtual public Algorithm {
		OFEC_CONCRETE_INSTANCE(LRDE)
	protected:
		Real m_repulsion_radius;
		Real m_repulsion_probabilty;
		size_t m_max_stagant_iters;
		size_t m_pop_size;
		PopulationDE<> m_pop;
		std::list<Solution<>> m_archive;
		std::unique_ptr<Solution<>> m_best;

		void addInputParameters();
		void run_(Environment *env) override;
		virtual bool reject(const Solution<> &sol, Environment *env) const;
		void updateArchive(Environment *env);
		void initializePopulation(Environment *env);
	};
}

#endif // !OFEC_LOCAL_REPULSION_DE_H
