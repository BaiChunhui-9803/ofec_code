#include "s_lshade_dp.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../core/environment/environment.h"
#include <algorithm>



#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {
	
    /*
  Return random value from Cauchy distribution with mean "mu" and variance "gamma"
  http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Cauchy
*/
    double cauchy_g(double mu, double gamma, Random* rnd) {
        return mu + gamma * tan(OFEC_PI * (rnd->uniform.next() - 0.5));
    }

    /*
  Return random value from normal distribution with mean "mu" and variance "gamma"
  http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Gauss
*/
    inline double gauss(double mu, double sigma, Random* rnd) {
        return mu + sigma * sqrt(-2.0 * log(rnd->uniform.next())) * sin(2.0 * OFEC_PI * rnd->uniform.next());
    }

    /*
  For each dimension j, if the mutant vector element v_j is outside the boundaries [x_min , x_max], we applied this bound handling method
  If you'd like to know that precisely, please read:
  J. Zhang and A. C. Sanderson, "JADE: Adaptive differential evolution with optional external archive,"
  IEEE Tran. Evol. Comput., vol. 13, no. 5, pp. 945¨C958, 2009.
 */
    void S_LSHADE_DP::modifySolutionWithParentMedium(Individual& child, const Individual& parent, Environment* env)
    {

        ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
        int DIM = con_pro->numberVariables();
        const auto& curDomain = con_pro->domain();


        for (int j = 0; j < DIM; j++)
        {
            if (child.variable()[j] < curDomain[j].limit.first)
            {
                child.variable()[j] = (curDomain[j].limit.first + parent.variable()[j]) / 2.0;
            }
            else if (child.variable()[j] > curDomain[j].limit.second)
            {
                child.variable()[j] = (curDomain[j].limit.second+ parent.variable()[j]) / 2.0;
            }
        }

    }

    void S_LSHADE_DP::reducePopulationWithSort(std::vector<Individual>& pop, std::vector<int>& stagnation) {
        int worst_ind;

        for (int i = 0; i < reduction_ind_num; i++)
        {
            worst_ind = 0;
            for (int j = 1; j < pop.size(); j++)
            {
                if (pop[j].fitness() < pop[worst_ind].fitness()) {
                    worst_ind = j;
                }
            }
            pop.erase(pop.begin() + worst_ind);
            stagnation.erase(stagnation.begin() + worst_ind);
            --m_pop_size;
        }
    }

    void S_LSHADE_DP::operateCurrentToPBest1BinWithArchive(const std::vector<Individual>& pop, 
        Individual& child, int& target, int& p_best_individual, double& scaling_factor, double& cross_rate, const std::vector<Individual>& archive, int& arc_ind_count, Environment* env) {
        ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
        int DIM = con_pro->numberVariables();
        int r1(0), r2(0);

        do
        {
            r1 = m_random->uniform.nextNonStd<int>(0, m_pop_size);
            //r1 = rand() % m_pop_size;
        } while (r1 == target);
        do
        {
            r2 = m_random->uniform.nextNonStd<int>(0, m_pop_size+ arc_ind_count);
            //r2 = rand() % (m_pop_size + arc_ind_count);
        } while ((r2 == target) || (r2 == r1));

        //int random_variable = rand() % problem_size;
        int random_variable = m_random->uniform.nextNonStd<int>(0, DIM);
        if (r2 >= m_pop_size)
        {
            r2 -= m_pop_size;
            for (int i = 0; i < DIM; i++)
            {
                if ((m_random->uniform.next()< cross_rate) || (i == random_variable))
                {
                    child.variable()[i] = pop[target].variable()[i] + scaling_factor * (pop[p_best_individual].variable()[i] - pop[target].variable()[i]) + scaling_factor * (pop[r1].variable()[i] - archive[r2].variable()[i]);
                }
                else
                {
                    child.variable()[i] = pop[target].variable()[i];
                }
            }
        }
        else
        {
            for (int i = 0; i < DIM; i++)
            {
                if ((m_random->uniform.next() < cross_rate) || (i == random_variable))
                {
                    child.variable()[i] = pop[target].variable()[i] + scaling_factor * (pop[p_best_individual].variable()[i] - pop[target].variable()[i]) + scaling_factor * (pop[r1].variable()[i] - pop[r2].variable()[i]);
                }
                else
                {
                    child.variable()[i]= pop[target].variable()[i];
                }
            }
        }

        //If the mutant vector violates bounds, the bound handling method is applied
        modifySolutionWithParentMedium(child, pop[target],env);
    }
    void S_LSHADE_DP::operateTarget1BinWithArchive(const std::vector<Individual>& pop, 
        Individual& child, int& target, double& scaling_factor, double& cross_rate, 
        const std::vector<Individual>& archive, int& arc_ind_count, Environment* env) {
        ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
        int DIM = con_pro->numberVariables();
        int r1, r2;
        do
        {
            r1 = m_random->uniform.nextNonStd<int>(0, m_pop_size);
        } while (r1 == target);
        do
        {
            r2 = m_random->uniform.nextNonStd<int>(0, m_pop_size + arc_ind_count);
        } while ((r2 == target) || (r2 == r1));

        int random_variable = m_random->uniform.nextNonStd<int>(0, DIM);

        if (r2 >= pop.size())
        {
            r2 -= pop.size();
            for (int i = 0; i < DIM; i++)
            {
                if ((m_random->uniform.next()< cross_rate) || (i == random_variable))
                {
                    child.variable()[i] = pop[target].variable()[i] + scaling_factor * (pop[r1].variable()[i] - archive[r2].variable()[i]);
                }
                else
                {
                    child.variable()[i] = pop[target].variable()[i];
                }
            }
        }
        else
        {
            for (int i = 0; i < DIM; i++)
            {
                if ((m_random->uniform.next() < cross_rate) || (i == random_variable))
                {
                    child.variable()[i] = pop[target].variable()[i] + scaling_factor * (pop[r1].variable()[i] - pop[r2].variable()[i]);
                }
                else
                {
                    child.variable()[i] = pop[target].variable()[i];
                }
            }
        }

        //If the mutant vector violates bounds, the bound handling method is applied
        modifySolutionWithParentMedium(child, pop[target],env);
    }

    void S_LSHADE_DP::addInputParameters() {
        m_input_parameters.add("population size", new RangedSizeT(m_pop_size, 5, 1000, 100));
        m_input_parameters.add("memory size", new RangedSizeT(m_memory_size, 0,100, 6));
        m_input_parameters.add("arc rate", new RangedReal(m_arc_rate, 0, 100, 2.6));
        m_input_parameters.add("p best rate", new RangedReal(m_p_best_rate, 0.0001, 1.0, 0.11));
    }

    

    void S_LSHADE_DP::initialize_(Environment* env)
    {
        Algorithm::initialize_(env);
        m_arc_size= (int)round(m_pop_size * m_arc_rate);

        m_eval_fun = [](ofec::SolutionBase& sol, ofec::Environment* env) {
            using namespace ofec;
            sol.evaluate(env, true);
            ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
            sol.setFitness(pos * sol.objective(0));
            };

    }

    void S_LSHADE_DP::run_(Environment* env)
    {
        using namespace std;
        ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
        int DIM = con_pro->numberVariables();
        auto& curdomain = con_pro->domain();

        std::vector<Individual> pop(m_pop_size,
            { con_pro->numberObjectives(), con_pro->numberConstraints(),con_pro->numberVariables()});
        std::vector<Individual> children(m_pop_size,
            { con_pro->numberObjectives(), con_pro->numberConstraints(),con_pro->numberVariables() });

        vector<int> stagnation(m_pop_size, 0);

        for (auto& it : pop) {
            it.initialize(env, m_random.get());
            m_eval_fun(it, env);
        }
        updateDatum(pop, env);
        Individual bsf_solution(con_pro->numberObjectives(), con_pro->numberConstraints(), con_pro->numberVariables());
        bsf_solution.setFitness(std::numeric_limits<double>::lowest());
        updateBsf(pop, bsf_solution);

        ////////////////////////////////////////////////////////////////////////////

//for external archive
        int arc_ind_count = 0;
        int random_selected_arc_ind;
        vector<Individual> archive(m_arc_size,
            { con_pro->numberObjectives(), con_pro->numberConstraints(),con_pro->numberVariables() });

        vector<int> num_success_params(M);
        vector<vector<double>> success_sf(M);
        vector<vector<double>> dif_fitness(M);

        vector<vector<double>> memory_sf;
        for (int i = 0; i < M; i++)
        {
            vector<double> mem_f_i(m_memory_size, 0.5);
            memory_sf.push_back(mem_f_i);
        }

        double temp_sum_sf;
        double sum;
        double weight;

        vector<double> improve_fitness(M, 0.0);
        vector<int> consumed_fes(M, 0);
        int best_op = 0; //rand() % M;

        //memory index counter
        vector<int> memory_pos(M, 0);

        //for new parameters sampling
        double mu_sf;
        int random_selected_period;
        std::vector<int> mut_op(m_pop_size,0);
        //int* mut_op = (int*)malloc(sizeof(int) * m_pop_size);
        std::vector<double> pop_sf (m_pop_size,0);
        std::vector<double> pop_cr(m_pop_size, 0);

        //for current-to-pbest/1
        int p_best_ind;
        int p_num = round(m_pop_size * m_p_best_rate);

        //Fitness* temp_fit = (Fitness*)malloc(sizeof(Fitness) * m_pop_size);

        // for linear population size reduction
        int max_pop_size = m_pop_size;
        int min_pop_size = 4;
        int plan_pop_size;

        //main loop
        int generation = 0;
        while (!terminating())
        {
            generation++;
            std::vector<int> sorted_array(m_pop_size);
            for (int i = 0; i < m_pop_size; i++)
                sorted_array[i] = i;
            
            std::sort(sorted_array.begin(), sorted_array.end(), [&](int a,int b) {
                return pop[a].fitness() > pop[b].fitness();
                });

            for (int target = 0; target < m_pop_size; target++)
            {
                // select a mutation operator
                int op_code = -1;
                double rand_dbl = m_random->uniform.next();
                for (int op = 0; op < M; op++)
                {
                    if (op * gamma <= rand_dbl && rand_dbl < (op + 1) * gamma)
                    {
                        op_code = op;
                        break;
                    }
                }

                if (op_code == -1)
                    op_code = best_op;

                mut_op[target] = op_code;
                consumed_fes[op_code]++;

                //In each generation, CR_i and F_i used by each individual x_i are generated by first selecting an index r_i randomly from [1, H]
                random_selected_period = m_random->uniform.nextNonStd<int>(0, m_memory_size);
                    //rand() % m_memory_size;
                mu_sf = memory_sf[op_code][random_selected_period];

                //generate F_i and repair its value
                do
                {
                    pop_sf[target] = cauchy_g(mu_sf, 0.1,m_random.get());
                } while (pop_sf[target] <= 0);

                if (pop_sf[target] > 1)
                    pop_sf[target] = 1;

                // generate CR_i
                if (evaluations() <= 0.5 * maximumEvaluations() || op_code == 0)
                {
                    pop_cr[target] = 0;
                }
                else
                {
                    pop_cr[target] = m_random->uniform.next();
                }

                if (op_code == 0)
                {
                    operateTarget1BinWithArchive(pop, children[target], target, pop_sf[target], pop_cr[target], archive, arc_ind_count,env);
                }
                else if (op_code == 1)
                {
                    //p-best individual is randomly selected from the top m_pop_size *  p_i members
                    p_best_ind = sorted_array[m_random->uniform.nextNonStd<int>(0,p_num)];
                    operateCurrentToPBest1BinWithArchive(pop, children[target], target, p_best_ind, pop_sf[target], pop_cr[target], archive, arc_ind_count,env);
                }
            }

            // evaluate the children's fitness values
            for (auto& it : children) {
                m_eval_fun(it, env);
            }
            //evaluatePopulation(children, children_fitness);

            /////////////////////////////////////////////////////////////////////////
            //update the bsf-solution and check the current number of fitness evaluations
            // if the current number of fitness evaluations over the max number of fitness evaluations, the search is terminated
            // So, this program is unconcerned about L-SHADE algorithm directly
            updateDatum(children, env);
            updateBsf(children, bsf_solution);

       
            ////////////////////////////////////////////////////////////////////////////

            //generation alternation
            for (int i = 0; i < m_pop_size; i++)
            {
                if (children[i].fitness() == pop[i].fitness())
                {
                    pop[i] = children[i];
              /*      fitness[i] = children_fitness[i];
                    for (int j = 0; j < problem_size; j++)
                    {
                        pop[i][j] = children[i][j];
                    }*/

                    if (evaluations() >= 0.5 * maximumEvaluations())
                    {
                        int count = 0;
                        for (int j = 0; j < DIM; j++)
                        {
                            if (children[i].variable()[j] != pop[i].variable()[j])
                            {
                                count++;
                            }
                        }

                        if (count > 0)
                        {
                            stagnation[i] = 0;
                        }
                        else
                        {
                            stagnation[i]++;
                        }
                    }
                }
                else if (children[i].fitness()>  pop[i].fitness())
                {
                    //parent vectors x_i which were worse than the trial vectors u_i are preserved
                    if (m_arc_size > 1)
                    {
                        if (arc_ind_count < m_arc_size)
                        {
                            archive[arc_ind_count] = pop[i];
                          /*  for (int j = 0; j < DIM; j++)
                                archive[arc_ind_count].variable()[j] = pop[i].variable()[j];*/
                            arc_ind_count++;
                        }
                        //Whenever the size of the archive exceeds, randomly selected elements are deleted to make space for the newly inserted elements
                        else
                        {
                            random_selected_arc_ind = m_random->uniform.nextNonStd<int>(0, m_arc_size);
                            archive[random_selected_arc_ind] = pop[i];
                 /*           for (int j = 0; j < problem_size; j++)
                                archive[random_selected_arc_ind][j] = pop[i][j];*/
                        }
                    }

                    dif_fitness[mut_op[i]].push_back(fabs(pop[i].fitness() - children[i].fitness()));
                    improve_fitness[mut_op[i]] += fabs(pop[i].fitness() - children[i].fitness());

                    //fitness[i] = children_fitness[i];
                    pop[i] = children[i];
               /*     for (int j = 0; j < problem_size; j++)
                        pop[i][j] = children[i][j];*/

                    //successful parameters are preserved in S_F and S_CR
                    success_sf[mut_op[i]].push_back(pop_sf[i]);

                    stagnation[i] = 0;
                }
                else
                {
                    if (evaluations() >= 0.5 * maximumEvaluations())
                    {
                        stagnation[i]++;
                    }
                }
            }

            for (int i = 0; i < M; i++)
                num_success_params[i] = success_sf[i].size();

            // if numeber of successful parameters > 0, historical memories are updated
            for (int m = 0; m < M; m++)
            {
                if (num_success_params[m] > 0)
                {
                    memory_sf[m][memory_pos[m]] = 0;
                    temp_sum_sf = 0;
                    sum = 0;

                    for (int i = 0; i < num_success_params[m]; i++)
                        sum += dif_fitness[m][i];

                    //weighted lehmer mean
                    for (int i = 0; i < num_success_params[m]; i++)
                    {
                        weight = dif_fitness[m][i] / sum;

                        memory_sf[m][memory_pos[m]] += weight * success_sf[m][i] * success_sf[m][i];
                        temp_sum_sf += weight * success_sf[m][i];
                    }

                    memory_sf[m][memory_pos[m]] /= temp_sum_sf;

                    //increment the counter
                    memory_pos[m]++;
                    if (memory_pos[m] >= m_memory_size)
                        memory_pos[m] = 0;

                    //clear out the S_F, S_CR and delta fitness
                    success_sf[m].clear();
                    dif_fitness[m].clear();
                }
            }

            // update best mutation operator
            if (generation % NG == 0)
            {
                // cout << generation << " " << best_op << " " << consumed_fes[0] << "  " << consumed_fes[1] << endl;
                int new_best_op = -1;
                double best_improve_rate = 0;

                for (int i = 0; i < M; i++)
                {
                    if (consumed_fes[i] > 0)
                    {
                        double improve_rate = improve_fitness[i] / consumed_fes[i];

                        if (improve_rate > best_improve_rate)
                        {
                            best_improve_rate = improve_rate;
                            new_best_op = i;
                        }
                    }

                    consumed_fes[i] = 0;
                    improve_fitness[i] = 0;
                }

                if (new_best_op == -1)
                {
                    best_op = 0; //rand() % M;
                }
                else
                {
                    best_op = new_best_op;
                }
            }

            // calculate the population size in the next generation
            plan_pop_size = round((((min_pop_size - max_pop_size) / (double)maximumEvaluations() )* evaluations()) + max_pop_size);

            if (m_pop_size > plan_pop_size)
            {
                reduction_ind_num = m_pop_size - plan_pop_size;
                if (m_pop_size - reduction_ind_num < min_pop_size)
                    reduction_ind_num = m_pop_size - min_pop_size;

                reducePopulationWithSort(pop, stagnation);

                // resize the archive size
                m_arc_size = m_pop_size * m_arc_rate;
                if (arc_ind_count > m_arc_size)
                    arc_ind_count = m_arc_size;

                // resize the number of p-best individuals
                // m_p_best_rate = 0.4 + (0.2 - 0.4) * count_fes * 1.0 / max_num_evaluations;
                p_num = round(m_pop_size * m_p_best_rate);
                if (p_num <= 1)
                    p_num = 2;
            }

            // dynamic perturbation
            for (int i = 0; i < m_pop_size; i++)
            {
                if (stagnation[i] > 100)
                {
                    double alpha = evaluations() * 1.0 / maximumEvaluations();

                    for (int j = 0; j < DIM; j++)
                    {
                        if (m_random->uniform.next() <= 0.5)
                        {
                            double rand_x = curdomain[j].limit.first + m_random->uniform.next() * (curdomain[j].length);
                            pop[i].variable()[j] = alpha * pop[i].variable()[j] + (1.0 - alpha) * rand_x;
                        }
                    }
                    m_eval_fun(pop[i], env);
                    //evaluateIndividual(pop[i], &fitness[i]);
                    stagnation[i] = 0;
                }
            }

            updateDatum(pop, env);
        }


    }


    void S_LSHADE_DP::updateBsf(const std::vector<Individual>& pop, Individual& bsf) {
        for (auto& it : pop) {
            if (it.fitness() > bsf.fitness()) {
                bsf = it;
            }
        }
    }

    void S_LSHADE_DP::updateDatum(const std::vector<Individual>& pop, Environment* env) {

#ifdef OFEC_DATUM_MULTI_POP_H
        g_multi_pop.pops.clear();
        g_multi_pop.pops.resize(1);
        for (size_t i = 0; i < pop.size(); ++i) {
            g_multi_pop.pops[0].push_back(&pop[i]);
        }
        datumUpdated(env, g_multi_pop);
#endif
    }
}