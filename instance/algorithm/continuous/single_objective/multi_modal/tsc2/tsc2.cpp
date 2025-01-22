#include "tsc2.h"
#include "../../../../../../core/problem/continuous/continuous.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
    void TSC2::addInputParameters() {
        m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 20));
        m_input_parameters.add("number of gradations", new RangedSizeT(m_num_gradations, 1, 15, 2));
        m_input_parameters.add("percent of seeds", new RangedReal(m_percent_seeds, 0, 100, 20));
        m_input_parameters.add("crossover rate", new RangedReal(m_pc, 0, 1, 0.74));
        m_input_parameters.add("mutation rate", new RangedReal(m_pm, 0, 1, 0.65));
    }

    void TSC2::initialize_(Environment *env) {
        Algorithm::initialize_(env);
        m_ms = 0;
        for (size_t j = 0; j < env->problem()->numberVariables(); ++j) {
            auto &range = CAST_CONOP(env->problem())->range(j);
            m_ms += pow(range.second - range.first, 2);
        }
        m_tournament_size = 2;
        m_ms = sqrt(m_ms) * 2;
    }

	void TSC2::run_(Environment *env) {
        m_iteration = 0;
        m_gradations.resize(m_num_gradations);
        for (int i = 0; i < m_num_gradations; i++)
            m_gradations[i] = (double)(1 + i) / (double)(m_num_gradations + 1);
        m_positive = env->problem()->optimizeMode(0) == OptimizeMode::kMaximize ? 1 : -1;
		m_pop.resize(m_pop_size, env, env->problem()->numberVariables());
		m_pop.initialize(env, m_random.get());
		m_pop.evaluate(env);
        for (size_t i = 0; i < m_pop_size; ++i) {
            m_pop[i].setFitness(m_positive * m_pop[i].objective(0));
        }
		while (!terminating()) {
			determineSeeds(env);
#ifdef OFEC_DATUM_MULTI_POP_H
            g_multi_pop.pops.clear();
            g_multi_pop.pops.resize(m_seeds.size() + 1);
            for (size_t i = 0; i < m_pop_size; ++i) {
                g_multi_pop.pops[m_pop[i].seed_id + 1].push_back(&m_pop[i]);
            }
            datumUpdated(env, g_multi_pop);
#endif
			tournamentSelection(env);
			recombination(env);
			normalMutation(env);
			insertSeedsWithoutDuplicates(env);
			integrateSolutionsToExistingSeeds(env);
			if (countFreeSolutions() > 0) {
				integrateFreeSolutions(env);
			}
			m_iteration++;
		}
	}

	void TSC2::determineSeeds(Environment *env) {
		for (size_t i = 0; i < m_pop_size; i++) {
			m_pop[i].marked = false;
			m_pop[i].marked_as_seeds = false;
			if (m_iteration == 0)
				m_pop[i].seed_id = -1;
			m_pop[i].ids.clear();
		}
        int num_seeds = 0;
        decltype(m_pop) temp(m_pop_size, env, env->problem()->numberVariables());
        bool finished_processing = false;
        while (!finished_processing) {
            finished_processing = true;
            int best_id = 0;
            while ((m_pop[best_id].marked) && (best_id < m_pop_size - 1))
                best_id++;
            IndTSC2 best(m_pop[best_id]);
            for (int i = best_id + 1; i < m_pop_size; i++)
                if (m_pop[i].marked == false && best.fitness() < m_pop[i].fitness()) {
                    best = m_pop[i];
                    best_id = i;
                }
            m_pop[best_id].marked = true;
            if (num_seeds == 0) {
                temp[num_seeds] = best;
                temp[num_seeds].seed_id = num_seeds;
                temp[num_seeds].old_seed_id = best.seed_id;
                m_pop[best_id].seed_id = num_seeds;
                num_seeds++;
            }
            else {  //there are several seeds
                bool completely_new = true;  //true if we see that best follows an optimum different from ALL the other "seeds"
                for (int k = 0; k < num_seeds; k++) {
                    if (m_iteration > 0) {
                        if (best.seed_id == temp[k].old_seed_id) {
                            completely_new = false;
                            m_pop[best_id].seed_id = temp[k].seed_id;
                        }
                    }
                    else//it's the first generation or reordering of subpopulations is performed
                        if (!hillValley(best, temp[k], env)) {
                            completely_new = false;
                            m_pop[best_id].ids.push_back(temp[k].seed_id);
                        }
                }
                if (completely_new) {
                    if (num_seeds < m_pop_size * m_percent_seeds / 100) {
                        temp[num_seeds] = best;
                        temp[num_seeds].seed_id = num_seeds;
                        temp[num_seeds].old_seed_id = best.seed_id;
                        m_pop[best_id].seed_id = num_seeds;
                        m_pop[best_id].marked_as_seeds = true;
                        num_seeds++;
                    }
                    else {  //integrate best to the closest seed
                        Real minimum_distance = m_pop[best_id].variableDistance(temp[0], env);
                        int min_id = 0;
                        for (int j = 1; j < num_seeds; j++) {
                            Real new_distance = m_pop[best_id].variableDistance(temp[j], env);
                            if (new_distance < minimum_distance) {
                                minimum_distance = new_distance;
                                min_id = j;
                            }
                        }
                        m_pop[best_id].seed_id = temp[min_id].seed_id;
                    }
                }
            }
            for (int i = 0; i < m_pop_size; i++) {
                if (!m_pop[i].marked) {
                    finished_processing = false;
                }
            }
        }
        m_seeds.resize(num_seeds, env, env->problem()->numberVariables());
        for (int i = 0; i < num_seeds; i++) {
            m_seeds[i] = temp[i];
        }
        if (m_iteration == 0) {
            for (int i = 0; i < m_pop_size; i++) {
                if (!m_pop[i].marked_as_seeds) {
                    for (int j = 0; j < m_seeds.size(); j++) {
                        if (!containsSeed(m_pop[i], m_seeds[j].seed_id)) {
                            if (!hillValley(m_pop[i], m_seeds[j], env)) {
                                m_pop[i].ids.push_back(m_seeds[j].seed_id);
                            }
                        }
                    }
                }
            }
            assignIDs(env);
        }
	}

	void TSC2::tournamentSelection(Environment *env) {
        decltype(m_pop) new_population(m_pop_size, env, env->problem()->numberVariables());
        Real worst_fit = m_pop[0].fitness();
        for (size_t i = 1; i < m_pop.size(); ++i) {
            if (worst_fit > m_pop[i].fitness()) {
                worst_fit = m_pop[i].fitness();
            }
        }
        for (int i = 0; i < m_pop_size; i++) {
            //count_selected counts how many Solutions follow same optima as Solution i
            std::vector<int> id_selected(m_tournament_size), count_selected(m_tournament_size);
            std::vector<Real> evaluations(m_tournament_size);
            for (int k = 0; k < m_tournament_size; k++) {
                id_selected[k] = m_random->uniform.nextNonStd<int>(0, m_pop_size);
                count_selected[k] = 0;
            }
            for (int k = 0; k < m_tournament_size; k++) {
                for (int j = 0; j < m_pop_size; j++) {
                    if (m_pop[j].seed_id == m_pop[id_selected[k]].seed_id) {
                        count_selected[k]++;
                    }
                }
            }
            for (int k = 0; k < m_tournament_size; k++) {
                evaluations[k] = (1 + (m_pop[id_selected[k]].fitness() - worst_fit)) / (Real)count_selected[k];
            }
            int winner_id = id_selected[0];
            Real max_evaluation = evaluations[0];
            for (int k = 1; k < m_tournament_size; k++) {
                if (max_evaluation < evaluations[k]) {
                    max_evaluation = evaluations[k];
                    winner_id = id_selected[k];
                }
            }
            new_population[i] = m_pop[winner_id];
        }
        for (int i = 0; i < m_pop_size; i++) {
            m_pop[i] = new_population[i];
        }
	}

	void TSC2::recombination(Environment *env) {
        decltype(m_pop) parents(m_pop_size, env, env->problem()->numberVariables());
        int num_parents = 0;
        for (int i = 0; i < m_pop_size; i++)
            m_pop[i].selected_for_crossover = false;
        for (int i = 0; i < m_pop_size; i++) {
            if (m_random->uniform.next() < m_pc) {
                parents[num_parents] = m_pop[i];
                m_pop[i].selected_for_crossover = true;
                parents[num_parents++].former_index = i;
            }
        }
        if (num_parents % 2 != 0) {
            if ((m_random->uniform.next() < 0.5) && (num_parents < m_pop_size)) {   //add a parent
                int new_parent_id;
                do {
                    new_parent_id = m_random->uniform.nextNonStd<int>(0, m_pop_size);
                } while (m_pop[new_parent_id].selected_for_crossover == true);
                parents[num_parents] = m_pop[new_parent_id];
                parents[num_parents].former_id = new_parent_id;
                num_parents++;
                m_pop[new_parent_id].selected_for_crossover = true;
            }
            else {  //remove a parent
                int remove_parent = m_random->uniform.nextNonStd<int>(0, num_parents);
                m_pop[parents[remove_parent].former_index].selected_for_crossover = false;
                for (int j = remove_parent; j < num_parents - 1; j++) {
                    parents[j] = parents[j + 1];
                }
                num_parents--;
            }
        }   //made sure num_parents was even
        decltype(m_pop) new_population(m_pop_size, env, env->problem()->numberVariables());
        //add to the new population the Solutions that were not selected for crossover
        int count = 0;
        for (int i = 0; i < m_pop_size; i++) {
            if (m_pop[i].selected_for_crossover == false) {
                new_population[count++] = m_pop[i];
            }
        }

        for (int i = 0; i < num_parents; i += 2) {
            //the offspring
            IndTSC2 off(parents[i]);
            Real modify = m_random->uniform.next();
            for (int k = 0; k < env->problem()->numberVariables(); k++)
                off.variable()[k] = parents[i].variable()[k] + modify * (Real)(parents[i + 1].variable()[k] - parents[i].variable()[k]);
            off.evaluate(env);
            off.setFitness(m_positive * off.objective(0));
            off.seed_id = -1;
            //if both parents have the same seed_id, the offspring will inherit it
            if (parents[i].seed_id == parents[i + 1].seed_id) {
                off.seed_id = parents[i].seed_id;
            }
            if (parents[i].fitness() > parents[i + 1].fitness()) {
                new_population[count] = parents[i];
                count++;
                if (off.fitness() > parents[i + 1].fitness()) {
                    new_population[count] = off;
                    count++;
                }
                else {
                    new_population[count] = parents[i + 1];
                    count++;
                }
            }
            else {  //if(eval(parents[i]) <= eval(parents[i + 1]))
                new_population[count] = parents[i + 1];
                count++;
                if (off.fitness() > parents[i].fitness()) {
                    new_population[count] = off;
                    count++;
                }
                else {
                    new_population[count] = parents[i];
                    count++;
                }
            }
        }
        //replace the population
        for (int i = 0; i < m_pop_size; i++) {
            m_pop[i] = new_population[i];
        }
	}

	void TSC2::normalMutation(Environment *env) {
        for (int i = 0; i < m_pop_size; i++) {
            IndTSC2 temp(m_pop[i]);
            for (int j = 0; j < env->problem()->numberVariables(); j++) {
                auto &range = CAST_CONOP(env->problem())->range(j);
                do {
                    Real mutate_or_not = m_random->uniform.next();
                    if (mutate_or_not < m_pm) {
                        temp.variable()[j] = m_pop[i].variable()[j] + m_ms * m_random->normal.next();
                    }
                    else {
                        temp.variable()[j] = m_pop[i].variable()[j];
                    }
                } while ((temp.variable()[j] > range.second) || (temp.variable()[j] < range.first));
            }
            bool changed = false;   //verify if there is any change to m_pop[i]
            for (int k = 0; k < env->problem()->numberVariables(); k++) {
                if (temp.variable()[k] != m_pop[i].variable()[k]) {
                    changed = true;
                }
            }
            if (changed) {
                temp.evaluate(env);
                temp.setFitness(m_positive * temp.objective(0));
                if (temp.fitness() > m_pop[i].fitness()) {
                    m_pop[i] = temp;
                    m_pop[i].seed_id = -1;
                }
            }   //end if(changed)
        }   //end for(int i = 0; i < m_pop_size; i++)
	}

	void TSC2::insertSeedsWithoutDuplicates(Environment *env) {
        for (int i = 0; i < m_pop_size; i++) {
            m_pop[i].marked = false;
        }
        //insert the seeds
        for (int i = 0; i < m_seeds.size(); i++) {
            //verify whether seed i already exists in the population
            int j1 = 0;
            bool already_exist = false;
            while ((j1 < m_pop_size) && (already_exist == false)) {
                if (m_pop[j1].seed_id == m_seeds[i].seed_id && already_exist == false) {
                    bool identical = true;
                    int k = 0;
                    if (m_seeds[i].fitness() != m_pop[j1].fitness()) {
                        identical = false;
                    }
                    else {
                        while (k < env->problem()->numberVariables() && identical) {
                            if (m_seeds[i].variable()[k] != m_pop[j1].variable()[k])
                                identical = false;
                            k++;
                        }
                    }
                    if (identical) {
                        already_exist = true;
                    }
                }
                j1++;
            }
            if (already_exist) {
                m_pop[j1 - 1].marked = true;
            }
            else {  //conserve the seed as it does not exist in the current population
                int num_neighbours = 0;
                decltype(m_pop) temp(m_pop_size, env, env->problem()->numberVariables());
                for (int j = 0; j < m_pop_size; j++) {
                    if (m_pop[j].seed_id == m_seeds[i].seed_id) {
                        temp[num_neighbours] = m_pop[j];
                        temp[num_neighbours++].former_id = j;
                    }
                }
                decltype(m_pop) temp1(num_neighbours, env, env->problem()->numberVariables());
                for (int s = 0; s < num_neighbours; s++) {
                    temp1[s] = temp[s];
                }
                if (num_neighbours > 0) {
                    IndTSC2 worst(temp1[0]);
                    for (int j = 1; j < num_neighbours; j++)
                        if (worst.fitness() > temp1[j].fitness()) {
                            worst = temp1[j];
                        }
                    if (worst.fitness() < m_seeds[i].fitness()) {
                        m_pop[worst.former_id] = m_seeds[i];
                    }
                    m_pop[worst.former_id].seed_id = m_seeds[i].seed_id;
                    m_pop[worst.former_id].marked = true;
                }
                else {  //if there are no neighbours for m_seeds[i]
                    int unmarked_id = 0;
                    while (m_pop[unmarked_id].marked && unmarked_id < m_pop_size)
                        unmarked_id++;
                    IndTSC2 worst(m_pop[unmarked_id]);
                    worst.former_id = unmarked_id;
                    for (int j = unmarked_id + 1; j < m_pop_size; j++) {
                        if (!m_pop[j].marked && worst.fitness() > m_pop[j].fitness()) {
                            worst  = m_pop[j];
                            worst.former_id = j;
                        }
                    }
                    m_pop[worst.former_id] = m_seeds[i];
                    m_pop[worst.former_id].marked = true;
                }   //end else - there are no neighbours for m_seeds[i]
            }   //end if(!already_exist)
        }   //end for(int i = 0; i < m_seeds.size(); i++)
	}

	void TSC2::integrateSolutionsToExistingSeeds(Environment *env) {
        for (int i = 0; i < m_pop_size; i++) {
            if (m_pop[i].seed_id == -1) {   //if ind i is offspring
                for (int iq = 0; iq < m_seeds.size(); iq++) {
                    m_seeds[iq].found = false;
                }
                int start_seed = 0;
                while (start_seed < m_seeds.size() && m_pop[i].seed_id == -1) {
                    start_seed = 0;
                    while (start_seed < m_seeds.size() && m_seeds[start_seed].found) {
                        start_seed++;
                    }
                    if (start_seed < m_seeds.size()) {
                        Real minimum_distance = m_pop[i].variableDistance(m_seeds[start_seed], env);
                        int id = start_seed;    //the id of the seed that is closer to ind i
                        for (int j = id + 1; j < m_seeds.size(); j++) {
                            Real distance = m_pop[i].variableDistance(m_seeds[j], env);
                            if (!m_seeds[j].found && minimum_distance > distance) {
                                minimum_distance = distance;
                                id = j;
                            }
                        }
                        m_seeds[id].found = true;
                        if (!hillValley(m_pop[i], m_seeds[id], env)) {
                            m_pop[i].seed_id = m_seeds[id].seed_id;
                        }
                    }
                }
            }
        }
	}

	size_t TSC2::countFreeSolutions() {
        size_t num_free_Solutions = 0;
        for (int i = 0; i < m_pop_size; i++) {
            if (m_pop[i].seed_id == -1) {
                num_free_Solutions++;
            }
        }
        return num_free_Solutions;
	}

	void TSC2::integrateFreeSolutions(Environment *env) {
        if (m_seeds.size() < m_pop_size * m_percent_seeds / 100) {  //if you can add more seeds
            int num_free_Solutions = 0;
            for (int i = 0; i < m_pop_size; i++) {
                if (m_pop[i].seed_id == -1) {
                    m_pop[i].marked = false;
                    num_free_Solutions++;
                }
                else {
                    m_pop[i].marked = true;
                }
            }
            int num_seeds = 0;
            decltype(m_pop) temp(num_free_Solutions, env, env->problem()->numberVariables());
            bool finished_processing = false;
            while (!finished_processing) {
                finished_processing = true;
                int best_id = 0;
                while ((m_pop[best_id].marked) && (m_pop[best_id].seed_id != -1) && (best_id < m_pop_size - 1))
                    best_id++;
                IndTSC2 best(m_pop[best_id]);
                for (int i = best_id + 1; i < m_pop_size; i++)
                    if (m_pop[i].marked == false && m_pop[i].seed_id == -1 && best.fitness() < m_pop[i].fitness()) {
                        best = m_pop[i];
                        best_id = i;
                    }
                m_pop[best_id].marked = true;
                bool completely_new = true; // true if we see that best follows an optimum different from ALL the other "seeds"
                if (completely_new && (num_seeds > 0)) {    // verify if current ind does not belong to a newly added seed
                    for (int k = 0; k < num_seeds; k++) {
                        if (!hillValley(best, temp[k], env)) {
                            completely_new = false;
                            m_pop[best_id].seed_id = temp[k].seed_id;
                        }
                    }
                }
                if (completely_new) {   // then it is a new seed
                    temp[num_seeds] = best;
                    temp[num_seeds].seed_id = m_seeds.size() + num_seeds;   // the ID has to be higher than existing ones
                    m_pop[best_id].seed_id = m_seeds.size() + num_seeds;
                    num_seeds++;
                }
                for (int i = 0; i < m_pop_size; i++) {
                    if (m_pop[i].seed_id == -1 && !m_pop[i].marked) {
                        finished_processing = false;
                    }
                }
            }
            int new_num_seeds = num_seeds + m_seeds.size();
            decltype(m_seeds) new_seeds(new_num_seeds, env, env->problem()->numberVariables());
            for (int j = 0; j < m_seeds.size(); j++) {
                new_seeds[j] = m_seeds[j];
            }
            int ii = 0;
            for (int j = m_seeds.size(); j < new_num_seeds; j++) {
                new_seeds[j] = temp[ii++];
            }
            m_seeds.resize(new_num_seeds, env, env->problem()->numberVariables());
            for (int i = 0; i < new_num_seeds; i++) {
                m_seeds[i] = new_seeds[i];
            }
        }
        else {  //if you cannot add more seeds, just attach the -1 Solutions to their closest seeds
            for (int i = 0; i < m_pop_size; i++) {
                if (m_pop[i].seed_id == -1) {
                    int min_id = 0;
                    Real minimum_distance = m_pop[i].variableDistance(m_seeds[0], env);
                    for (int j = 1; j < m_seeds.size(); j++) {
                        Real distance = m_pop[i].variableDistance(m_seeds[j], env);
                        if (minimum_distance > distance) {
                            minimum_distance = distance;
                            min_id = j;
                        }
                    }
                    m_pop[i].seed_id = m_seeds[min_id].seed_id;
                }
            }
        }
	}

    bool TSC2::containsSeed(const IndTSC2 &ind, int seed_id) {
        bool ret = false;
        for (int i = 0; i < ind.ids.size(); i++)
            if (ind.ids[i] == seed_id)
                ret = true;
        return ret;
    }

    bool TSC2::hillValley(const IndTSC2 &ind1, const IndTSC2 &ind2, Environment *env) {
        bool ret = false;
        bool identical_points = true;
        int i1 = 0;
        while (i1 < env->problem()->numberVariables() && identical_points) {
            if (ind1.variable()[i1] != ind2.variable()[i1])
                identical_points = false;
            i1++;
        }
        if (identical_points) {
            ret = false;
        }
        else {
            Real min_eval = ind1.fitness(); //ind1.fitness();
            if (min_eval > ind2.fitness())
                min_eval = ind2.fitness();
            decltype(m_pop) interior_Solution(m_num_gradations, env, env->problem()->numberVariables());
            int i = 0;
            while (i < m_num_gradations && ret == false) {
                for (int j = 0; j < env->problem()->numberVariables(); j++)
                    interior_Solution[i].variable()[j] = ind1.variable()[j] + (ind2.variable()[j] - ind1.variable()[j]) * m_gradations[i];
                interior_Solution[i].evaluate(env);
                interior_Solution[i].setFitness(m_positive * interior_Solution[i].objective(0));
                if (min_eval > interior_Solution[i].fitness()) {
                    ret = true;
                }
                i++;
            }
        }   //end the case when the two points are not identical
        return ret;
    }

	void TSC2::assignIDs(Environment *env) {
		for (int i = 0; i < m_pop_size; i++) {
			if (!m_pop[i].marked_as_seeds) {
                if (m_pop[i].ids.size() == 1) {
                    m_pop[i].seed_id = m_pop[i].ids[0];
                }
                else {
                    if (m_pop[i].ids.size() > 1) {//then find the closest seed and assign the m_pop[i] to it
                        Real t = m_seeds[m_pop[i].ids[0]].variableDistance(m_pop[i], env);
                        int searched_id = 0;
                        for (int j = 1; j < m_pop[i].ids.size(); j++) {
                            Real distance = m_seeds[m_pop[i].ids[j]].variableDistance(m_pop[i], env);
                            if (t > distance) {
                                t = distance;
                                searched_id = j;
                            }
                        }
                        m_pop[i].seed_id = m_pop[i].ids[searched_id];
                    }
                }
			}
		}
	}
}
