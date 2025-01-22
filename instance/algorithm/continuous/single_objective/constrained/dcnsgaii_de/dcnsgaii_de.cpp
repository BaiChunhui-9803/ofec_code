#include "dcnsgaii_de.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../record/rcr_vec_real.h"
#include <numeric>

using namespace std::placeholders;

namespace ofec {
	DCNSGAII_DE_pop::DCNSGAII_DE_pop(size_t size_pop, Problem *pro) : Population<DCMOEA_ind<IndDE>>(size_pop, pro, CAST_CONOP(pro)->numberVariables()), DCMOEA<DCMOEA_ind<IndDE>>(), m_rand_seq(size_pop) {
		for (auto& i : m_individuals) {
			i->resize_vio_obj(1);
			m_offspring.emplace_back(new IndividualType(*i));
			m_temp_pop.emplace_back(new IndividualType(*i));
		}
		std::iota(m_rand_seq.begin(), m_rand_seq.end(), 0);
		m_scaling_factor = 0.5;
		m_crossover_rate = 0.6;
		m_recombine_strategy = RecombineDE::Binomial;
		/*	CONTINUOUS_CAST->set_eval_monitor_flag(true);*/
	}

	void DCNSGAII_DE_pop::initialize(Problem *pro, Algorithm *alg, Random *rnd) {
		Population::initialize(pro, rnd);
		DCMOEA::initialize(CAST_CONOP(pro)->numberConstraints());
		evaluate(pro, alg);
		calculate_initial_max_violation(m_individuals);
		m_e = m_max_G;
		calculate_violation_objective(m_individuals);
		mark_Solution_efeasible(m_individuals);
		Real dimension = CAST_CONOP(pro)->numberVariables();
		m_R = 0.5 * pow(2 * dimension / (2 * m_individuals.size() * OFEC_PI), 1.0 / dimension);



	}

	void DCNSGAII_DE_pop::sort(Problem *pro) {
		std::vector<int> ranks;
		std::function<Dominance(IndividualType* const&, IndividualType* const&)> comp = std::bind(&DCNSGAII_DE_pop::Pareto_compare, this, _1, _2, pro);
		std::vector<IndividualType*> pop;
		for (auto& i : m_individuals)
			pop.emplace_back(i.get());
		nd_sort::fastSort<IndividualType*>(pop, ranks, comp);
		/*for (size_t i = 0; i < pop.size(); i++) {
			if (ranks[i] == 0) {
				std::cout << "pop[" << i << "]: ";
				for (size_t j = 0; j < CONTINUOUS_CAST->objective_size(); ++j)
					std::cout << pop[i]->objective(j) << "\t\t";
				std::cout << pop[i]->get_vio_obj()[0] << "\t\t" << pop[i]->get_vio_obj()[1] <<std::endl;
			}
		}*/
	}

	int DCNSGAII_DE_pop::evolve(Problem *pro, Algorithm *alg, Random *rnd)
	{
		int tag(0);
		if (judge_population_efeasible(m_individuals)) {
			m_k++;
			if (m_k > m_max_K + 1) {
				//return Terminate;
				//m_flag = true;
			}
			else {
				reduce_boundary();
				reduce_radius(pro);
				mark_Solution_efeasible(m_individuals);
			}
		}

		//generate offspring pop
		for (size_t i = 0; i < m_individuals.size(); i++) {
			std::vector<int> ridx;
			rnd->uniform.shuffle(m_rand_seq.begin(), m_rand_seq.end());
			for (int idx : m_rand_seq) {
				if (idx != i)					ridx.emplace_back(idx);
				if (ridx.size() == 3)			break;
			}
			m_individuals[i]->mutate(m_scaling_factor, m_individuals[ridx[0]].get(), m_individuals[ridx[1]].get(), m_individuals[ridx[2]].get(), pro);
			m_individuals[i]->recombine(m_crossover_rate, m_recombine_strategy, rnd, pro);
			m_offspring[i] = m_individuals[i]->trial();
		}

		//evaluate offspring pop
		for (auto& i : m_offspring)
			i->evaluate(pro, alg);

		calculate_violation_objective(m_offspring);
		mark_Solution_efeasible(m_offspring);

		std::vector<IndividualType*> pop;
		for (auto& i : m_individuals)
			pop.emplace_back(i.get());
		for (auto& i : m_offspring)
			pop.emplace_back(i.get());
		caculate_nichecount(pop, pro);
		//select next pop
		select_next_parent_population(pop, pro);

		if (tag == 0 )
			m_iteration++;

		return tag;
	}

	void DCNSGAII_DE_pop::select_next_parent_population(std::vector<IndividualType*>& pop, Problem *pro)
	{
		/* Nondominated sorting based on e-Pareto domination  */
		std::vector<int> ranks;
		std::function<Dominance(IndividualType* const&, IndividualType* const&)> comp = std::bind(&DCNSGAII_DE_pop::e_Pareto_compare, this, _1, _2,pro);
		nd_sort::fastSort<IndividualType*>(pop, ranks, comp);
		for (size_t i = 0; i < pop.size(); i++) {
			pop[i]->setFitness(ranks[i]);
		}

		size_t cur_rank = 0;
		size_t id_ind = 0;
		while (true) {
			int count = 0;
			for (size_t i = 0; i < pop.size(); i++)
				if (pop[i]->fitness() == cur_rank)
					count++;
			int size2 = id_ind + count;
			if (size2 > this->m_individuals.size()) {
				break;
			}
			for (size_t i = 0; i < pop.size(); i++)
				if (pop[i]->fitness() == cur_rank) {
					*m_temp_pop[id_ind] = *pop[i];
					++id_ind;
				}
			cur_rank++;
			if (id_ind >= this->m_individuals.size()) break;
		}
		if (id_ind < pop.size()) {
			std::vector<int> list;	// save the Solutions in the overflowed front
			for (size_t i = 0; i < pop.size(); i++)
				if (pop[i]->fitness() == cur_rank)
					list.push_back(i);
			int s2 = list.size();
			std::vector<Real> density(s2);
			std::vector<Real> obj(s2);
			std::vector<int> idx(s2);
			std::vector<int> idd(s2);
			for (size_t i = 0; i < s2; i++) {
				idx[i] = i;
				density[i] = 0;
			}
			for (size_t j = 0; j < 2; j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = pop[list[i]]->get_vio_obj()[0];
				}
				mergeSort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			for (size_t j = 0; j < CAST_CONOP(pro)->numberObjectives(); j++) {
				for (size_t i = 0; i < s2; i++) {
					idd[i] = i;
					obj[i] = pop[list[i]]->objective()[j];
				}
				mergeSort(obj, s2, idd, true, 0, s2 - 1, s2);
				density[idd[0]] += -1.0e+30;
				density[idd[s2 - 1]] += -1.0e+30;
				for (int k = 1; k < s2 - 1; k++)
					density[idd[k]] += -(obj[idd[k]] - obj[idd[k - 1]] + obj[idd[k + 1]] - obj[idd[k]]);
			}
			idd.clear();
			obj.clear();
			int s3 = this->m_individuals.size() - id_ind;
			mergeSort(density, s2, idx, true, 0, s2 - 1, s3);
			for (size_t i = 0; i < s3; i++) {
				*m_temp_pop[id_ind] = *pop[list[idx[i]]];
				++id_ind;
			}
			density.clear();
			idx.clear();
			list.clear();
		}
		for (size_t i = 0; i < m_temp_pop.size(); i++) {
			*m_individuals[i] = *m_temp_pop[i];
		}
	}

	void DCNSGAII_DE_pop::reduce_radius(Problem *pro) {
		Real production = 1.0;
		for (int i = 0; i < CAST_CONOP(pro)->numberVariables(); ++i) {
			production *= (CAST_CONOP(pro)->range(i).second - CAST_CONOP(pro)->range(i).first);
		}
		Real z1 = m_z / (pow(production, (1.0 / CAST_CONOP(pro)->numberVariables())));
		Real C = m_R + z1;
		Real c = sqrt(abs(log2(C / z1) / log2(exp(1))));
		Real D = m_max_K / c;
		Real q = m_k / D;
		Real f = C * exp(-pow(q, 2)) - z1;
		if (abs(f) < m_Nearzero)
			f = 0.0;
		m_r = f;
	}

	void DCNSGAII_DE_pop::caculate_nichecount(std::vector<IndividualType*>& pop, Problem *pro)
	{
		size_t size_pop = pop.size();
		std::vector<Real> nicheCount(size_pop, 0.0);
		for (int i = 0; i < size_pop; ++i) {
			for (int j = i + 1; j < size_pop; ++j) {
				Real sum1 = 0.0;
				for (int k = 0; k < CAST_CONOP(pro)->numberVariables(); ++k) {
					auto difference = pop[i]->variable()[k] - pop[j]->variable()[k];
					sum1 += (difference * difference);
				}
				Real aDist = sqrt(sum1);
				if (aDist < m_R) {
					nicheCount[i] += (1 - (aDist / m_r));
					nicheCount[j] += (1 - (aDist / m_r));
				}
			}
			//pop[i]->get_vio_obj()[1] = nicheCount[i];
		}
	}

	Dominance DCNSGAII_DE_pop::e_Pareto_compare(IndividualType* const& s1, IndividualType* const& s2, Problem *pro)
	{	/* One efeasible one in-efeasible */
		if (s1->get_efeasible() != s2->get_efeasible()) {
			if (s1->get_efeasible())
				return Dominance::kDominant;
			else
				return Dominance::kDominated;
		}

		/* Both efeasible */
		else if (s1->get_efeasible() && s2->get_efeasible()) {
			auto nor_obj_result = objectiveCompare<Real>(s1->objective(), s2->objective(), CAST_CONOP(pro)->optimizeMode());
			auto vio_obj_result = objectiveCompare<Real>(s1->get_vio_obj(), s2->get_vio_obj(), OptimizeMode::kMinimize);

			if (nor_obj_result == Dominance::kDominant && vio_obj_result == Dominance::kEqual)
				return Dominance::kDominant;
			if (nor_obj_result == Dominance::kEqual && vio_obj_result == Dominance::kDominant)
				return Dominance::kDominant;
			if (nor_obj_result == Dominance::kDominant && vio_obj_result == Dominance::kDominant)
				return Dominance::kDominant;

			if (nor_obj_result == Dominance::kDominated && vio_obj_result == Dominance::kEqual)
				return Dominance::kDominated;
			if (nor_obj_result == Dominance::kEqual && vio_obj_result == Dominance::kDominated)
				return Dominance::kDominated;
			if (nor_obj_result == Dominance::kDominated && vio_obj_result == Dominance::kDominated)
				return Dominance::kDominated;

			if (nor_obj_result == Dominance::kDominated && vio_obj_result == Dominance::kDominant)
				return Dominance::kNonDominated;
			if (nor_obj_result == Dominance::kDominant && vio_obj_result == Dominance::kDominated)
				return Dominance::kNonDominated;
			if (nor_obj_result == Dominance::kNonDominated || vio_obj_result == Dominance::kNonDominated)
				return Dominance::kNonDominated;

			if (nor_obj_result == Dominance::kEqual && vio_obj_result == Dominance::kEqual)
				return Dominance::kEqual;
		}
		//return objective_compare<Real>(s1->objective(), s2->objective(), CONTINUOUS_CAST->opt_mode());

	/* Both in-efeasible */
		else {
			if (s1->get_vio_obj()[0] < s2->get_vio_obj()[0])
				return Dominance::kDominant;
			else if (s1->get_vio_obj()[0] > s2->get_vio_obj()[0])
				return Dominance::kDominated;
			else
				return Dominance::kEqual;
		}
	}
	Dominance DCNSGAII_DE_pop::Pareto_compare(IndividualType* const& s1, IndividualType* const& s2, Problem *pro)
	{
		/* Both efeasible */
		if (s1->get_efeasible() && s2->get_efeasible()) {
			auto nor_obj_result = objectiveCompare<Real>(s1->objective(), s2->objective(), CAST_CONOP(pro)->optimizeMode());
			auto vio_obj_result = objectiveCompare<Real>(s1->get_vio_obj(), s2->get_vio_obj(), OptimizeMode::kMinimize);

			if (nor_obj_result == Dominance::kDominant && vio_obj_result == Dominance::kEqual)
				return Dominance::kDominant;
			if (nor_obj_result == Dominance::kEqual && vio_obj_result == Dominance::kDominant)
				return Dominance::kDominant;
			if (nor_obj_result == Dominance::kDominant && vio_obj_result == Dominance::kDominant)
				return Dominance::kDominant;

			if (nor_obj_result == Dominance::kDominated && vio_obj_result == Dominance::kEqual)
				return Dominance::kDominated;
			if (nor_obj_result == Dominance::kEqual && vio_obj_result == Dominance::kDominated)
				return Dominance::kDominated;
			if (nor_obj_result == Dominance::kDominated && vio_obj_result == Dominance::kDominated)
				return Dominance::kDominated;

			if (nor_obj_result == Dominance::kDominated && vio_obj_result == Dominance::kDominant)
				return Dominance::kNonDominated;
			if (nor_obj_result == Dominance::kDominant && vio_obj_result == Dominance::kDominated)
				return Dominance::kNonDominated;
			if (nor_obj_result == Dominance::kNonDominated || vio_obj_result == Dominance::kNonDominated)
				return Dominance::kNonDominated;

			if (nor_obj_result == Dominance::kEqual && vio_obj_result == Dominance::kEqual)
				return Dominance::kEqual;
		}
	}


	void DCNSGAII_DE::initialize_() {
		Algorithm::initialize_();
		auto& v = *m_param;
		m_pop_size = v.get<int>("population size");
		m_pop.reset();
	}

	void DCNSGAII_DE::run_() {
		initPop();
		while (!terminating()) {
			m_pop->evolve(m_problem.get(), this, m_random.get());
		}
		m_pop->sort(m_problem.get());
	}

	void DCNSGAII_DE::initPop() {
		m_pop.reset(new DCNSGAII_DE_pop(m_pop_size, m_problem.get()));
		m_pop->initialize(m_problem.get(), this, m_random.get());
		//m_pop->evaluate(m_problem.get(), this);
	}

	void DCNSGAII_DE::record() {
		std::vector<Real> entry;
		entry.push_back(m_evaluations);
		Real err = std::fabs(m_candidates.front()->objective(0) - CAST_CONOP(m_problem.get())->optima()->objective(0).at(0));
		entry.push_back(err);
		for (size_t i = 0; i < m_problem->numberConstraints(); i++) {
			Real cons = m_candidates.front()->constraint()[i];
			entry.push_back(cons);
		}
		Real vio = (*m_pop)[0].get_vio_obj()[0];
		entry.push_back(vio);
		dynamic_cast<RecordVectorReal*>(m_record.get())->addEntry(this, entry);
	}
}