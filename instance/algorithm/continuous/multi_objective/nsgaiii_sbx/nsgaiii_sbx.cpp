#include "nsgaiii_sbx.h"
#include "../../../../../core/problem/continuous/continuous.h"
#include "../../../../record/rcr_vec_real.h"
#include "../../../../../utility/metricsMOP/IGD.h"

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
	void NSGAIII_SBX::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;;
		m_pop_size = v.get<int>("population size");
		if (m_pop_size % 2)
			throw MyExcept("Population size of NSGAII should be even.");
		m_cr = v.has("crossover rate") ? v.get<Real>("crossover rate") : 0.9;
		m_mr = v.has("mutation rate") ? v.get<Real>("mutation rate") : 1.0 / CAST_CONOP(m_problem.get())->numberVariables();
		m_ceta = v.has("crossover eta") ? v.get<Real>("crossover eta") : 20.;
		m_meta = v.has("mutation eta") ? v.get<Real>("mutation eta") : 20.;
	}

	void NSGAIII_SBX::initPop() {
		auto size_var = CAST_CONOP(m_problem.get())->numberVariables();
		auto size_obj = CAST_CONOP(m_problem.get())->numberObjectives();
		m_pop.reset(new PopNSGAIII_SBX(m_pop_size, m_problem.get(), size_var, size_obj, CAST_CONOP(m_problem.get())->optimizeMode()));
		m_pop->setRate(m_cr, m_mr);
		m_pop->setEta(m_ceta, m_meta);
		m_pop->initialize_(m_problem.get(), m_random.get());
		m_pop->evaluate(m_problem.get(), this);
	}

	void NSGAIII_SBX::run_() {
		initPop();
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
			updateBuffer();
#endif
		}
	}


#ifdef OFEC_DEMO
	void NSGAIII_SBX::updateBuffer() {
		m_solution.clear();
		m_solution.resize(1);
		for (size_t i = 0; i < m_pop->size(); ++i)
			m_solution[0].push_back(&m_pop->at(i));
		ofec_demo::g_buffer->appendAlgBuffer(this);
	}
#endif

	void NSGAIII_SBX::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		std::vector<std::vector<Real>> ref_objs;
		for (size_t i = 0; i < m_problem->optimaBase()->numberObjectives(); ++i) {
			ref_objs.push_back(m_problem->optimaBase()->objective(i));
		}
		std::vector<std::vector<Real>> pop_objs;
		for (size_t i = 0; i < m_pop->size(); ++i) {
			pop_objs.push_back(m_pop->at(i).objective());
		}
		Real IGDv = IGD(ref_objs, pop_objs);
		entry.push_back(IGDv);
		dynamic_cast<RecordVectorReal*>(m_record.get())->addEntry(this, entry);
	}

	PopNSGAIII_SBX::PopNSGAIII_SBX(size_t size_pop, Problem *pro, size_t size_var, size_t size_obj, const std::vector<OptimizeMode>& opt_mode) :
		PopSBX<>(size_pop, pro), NSGAIII<Solution<>>(size_pop, size_obj, opt_mode, pro) {}

	void PopNSGAIII_SBX::initialize_(Problem *pro, Random *rnd) {
		initialize(pro, rnd);
		for (auto& i : m_individuals) {
			m_offspring.emplace_back(*i);
		}
		for (auto& i : m_individuals) {
			m_offspring.emplace_back(*i);
		}
	}

	int PopNSGAIII_SBX::evolve(Problem *pro, Algorithm *alg, Random *rnd) {
		if (m_individuals.size() % 2 != 0)
			throw "population size should be even @NSGAIII_SBXRealMu::evolveMO()";
		for (size_t i = 0; i < m_individuals.size(); i += 2) {
			std::vector<size_t> p(2);
			p[0] = tournamentSelection(pro, rnd);
			do { p[1] = tournamentSelection(pro, rnd); } while (p[1] == p[0]);
			crossover(p[0], p[1], m_offspring[i], m_offspring[i + 1], pro, rnd);
			mutate(m_offspring[i], pro, rnd);
			mutate(m_offspring[i + 1], pro, rnd);
		}
		for (auto& i : m_offspring) {
			i.evaluate(pro, alg);
		}
		for (size_t i = 0; i < m_individuals.size(); ++i) {
			m_offspring[i + m_individuals.size()] = *m_individuals[i];
		}
		survivorSelection(m_individuals, m_offspring, pro, rnd);


		m_iteration++;
		return kNormalEval;
	}
}