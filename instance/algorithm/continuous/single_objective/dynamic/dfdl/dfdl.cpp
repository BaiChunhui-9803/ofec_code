//
// Created by Mai Peng on 2022/11/16.
//
#ifndef DFDL_H_
#define DFDL_H_
#include "dfdl.h"
#include <map>
#include <memory>
#include <unordered_map>
#include <random>
#include "../../../../../../utility/clustering/nbc.h"
#include <exception>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {
    void DFDL::initialize_() {
        HGHEC::initialize_();
        MetricsDynamicConOEA::initialize_();
        auto& v = *m_param;;
        m_size_pop = v.get<int>("population size");
        m_alpha = v.has("alpha") ?
                  v.get<Real>("alpha") : 0.45;
        m_beta = v.has("beta") ?
                 v.get<Real>("beta") : 0.45;
        m_gama = v.has("gama") ?
                 v.get<Real>("gama") : 0.1;
        m_init_num_ssp = v.has("initial number of subspaces") ?
                       v.get<int>("initial number of subspaces") : 1;
        m_uncertainty_type = v.has("uncertainty type") ?
                             v.get<int>("obsolete factor") : UncertaintyType::linear;
        m_pop_start_type = v.has("Population Strategy") ? static_cast<DFDL::PopStart>(v.get<int>("Population Strategy")) : PopStart::Niching_Lips;
        m_clustering_type = v.has("Clustering Strategy") ? static_cast<DFDL::ClusteringType>(v.get<int>("Clustering Strategy")) : ClusteringType::NBC;
        m_min_active_hill_volume = 0.4 / m_problem->numberVariables();
//        m_obsoleteFactor = v.has("obsolete factor") ?
//                           v.get<int>("obsolete factor") : 1000;
    }


    void DFDL::initSolutions(int Solution_num, int type) {
        const size_t num_obj = CAST_CONOP(m_problem.get())->numberObjectives();
        const size_t num_con = CAST_CONOP(m_problem.get())->numberConstraints();
        const size_t num_var = CAST_CONOP(m_problem.get())->numberVariables();
        if(type == InitSolType::Random) {
            int originSize = m_new_indis.size();
            m_new_indis.resize(Solution_num + originSize);

            for (int i(originSize); i < m_new_indis.size(); ++i) {
                m_new_indis[i] = std::make_unique<Solution_type>(num_obj, num_con, num_var);
                m_new_indis[i]->initialize(i, m_problem.get(), m_random.get());
            }
        }
        else if (type == InitSolType::AtHighUncertainty) {
            m_new_indis.resize(Solution_num);
            int selected_count = 0;
            while(selected_count < Solution_num) {
                Real probability_total = 0.;
                int selected_ssp = 0;
                Real rnd = m_random->uniform.next();
                for (int i = 0; i < m_info_ssp.size(); ++i) {
                    probability_total += m_probability_SSP[i];
                    if(probability_total > rnd) {
                        selected_ssp = i;
                        break;
                    }
                }
                const auto & ssp_boundary = m_ssp_tree->getBox(selected_ssp);
                m_new_indis[selected_count] = std::make_unique<Solution_type>(num_obj, num_con, num_var);
                for (int i = 0; i < ssp_boundary.size(); ++i) {
                    m_new_indis[selected_count] -> variable()[i] = m_random->uniform.nextNonStd(ssp_boundary[i].first, ssp_boundary[i].second);
                }
                m_new_indis[selected_count]->initVelocity(m_problem.get(), m_random.get());
                m_new_indis[selected_count]->initVelocityMax(m_problem.get(), m_random.get());
                ++selected_count;
            }
        }
    }

    int DFDL::getNumberAboveMeanUncertainty() {
        int count = 0;
        Real meanUncertainty = 0.;
        for (const auto &ssp: m_info_ssp) {
            meanUncertainty += ssp.uncertainty;
        }
        meanUncertainty = meanUncertainty/m_info_ssp.size();
        for (const auto &ssp: m_info_ssp) {
            if(ssp.uncertainty > meanUncertainty) {
                ++count;
            }
        }
//        std::cout << count << std::endl;
        return count;
    }

    void DFDL::calculateProbabilitySSP() {
        m_probability_SSP.resize(m_info_ssp.size());
        Real temp = 0.;
        for (const auto &ssp: m_info_ssp) {
            temp += ssp.uncertainty;
        }
        for (int i = 0; i < m_info_ssp.size(); ++i) {
            m_probability_SSP[i] = m_info_ssp[i].uncertainty / temp;
        }
    }
    bool DFDL::checkDiversity(Problem *pro) {
        Real meanRadius = 0.;
        Real meanVelocity = 0.;
        int count = 0;
        for (const auto &pop: m_sub_pop) {
            if(pop->isActive() && !pop->getStagnantedFlag()) {
                pop->updateCurRadius(m_problem.get());
                pop->updateCurVelocity(m_problem.get());
                meanRadius += pop->getCurRadius();
                meanVelocity += pop->getCurVelocity();
                count++;
            }
        }
        if(count > 0) {
            meanRadius = meanRadius / (int) count;
            meanVelocity = meanVelocity / (int) count;
        }
        else return true;
        if(meanRadius < CAST_CONOP(m_problem.get())->domainArea() * 0.0025 * CAST_CONOP(m_problem.get())->numberVariables())
            return true;
        return false;
    }

    // Create Pop for Multi-Population
    int DFDL::createNewSwarms(int Solution_num, int type) {
        int rf = kNormalEval;
        initSolutions(Solution_num, type);
        HSLH<Solution_type> cluster(m_new_indis, m_problem.get());
        cluster.clustering(min_subpopSize, m_problem.get());
        for (int i(0); i < cluster.size(); ++i) {
            if (cluster[i].size() > m_min_num_inds) {
                std::unique_ptr<TypePop> new_pop = std::make_unique<TypePop>(cluster[i].size(), m_problem.get(), DFDLSwarm::PopulationType::Explore);
                auto iter = cluster[i].begin();
                for (int j(0); j < cluster[i].size(); ++j, ++iter) {
                    (*new_pop)[j] = *(iter->second);
                    rf = (*new_pop)[j].evaluate(m_problem.get(), this, true);
                    if (rf != kNormalEval) return rf;
                    archiveSolution((*new_pop)[j], TaskEval::kExplore);
                }
                new_pop -> initializeParameters(m_problem.get(), m_random.get());
                m_sub_pop.append(new_pop);
            } else {
                auto iter = cluster[i].begin();
                for (int j(0); j < cluster[i].size(); ++j, ++iter) {
                    m_rest_indis.emplace_back(new Solution_type(*(iter->second)));
                }
            }
        }
        m_new_indis.clear();
        return m_sub_pop.size();
    }

    // Create Pop for Niching
    int DFDL::createNichingSwarm(int Solution_num, int type) {
        int rf = kNormalEval;
        const size_t num_obj = CAST_CONOP(m_problem.get())->numberObjectives();
        const size_t num_con = CAST_CONOP(m_problem.get())->numberConstraints();
        const size_t num_var = CAST_CONOP(m_problem.get())->numberVariables();
        if(type == InitSolType::Random) {
            m_swarm_Lips = std::make_unique<SwarmLIP>(Solution_num, m_problem.get(), 500000);
            m_swarm_Lips->initialize(m_problem.get(), m_random.get());
            rf = m_swarm_Lips->evaluate(m_problem.get(), this);
            if (rf != kNormalEval) return rf;
            m_swarm_Lips->initPbest(m_problem.get());
        }
        else if(type == InitSolType::AtHighUncertainty) {
            m_new_indis.resize(Solution_num);
            int selected_count = 0;
            int tries = 0;
            std::map<int, bool> selected_map;
            while(selected_count < Solution_num) {
                tries++;
                bool success_flag = true;
                Real probability_total = 0.;
                int selected_ssp = 0;
                Real rnd = m_random->uniform.next();
                // 如果只有一个吸引域，则在吸引域内部产生个体 -> 与之前的做法非常类似
                if(m_num_hills == 1) {
                    for (int i = 0; i < m_info_ssp.size(); ++i) {
                        probability_total += m_probability_SSP[i];
                        if (probability_total > rnd) {
                            selected_ssp = i;
                            if (selected_map.find(i) == selected_map.end()) {
                                selected_map.insert(std::pair<size_t, bool>(i, true));
                                success_flag = true;
                            } else success_flag = false;
                            break;
                        }
                    }
                }
                else if(tries < m_init_num_ssp) {
                    for (int i = 0; i < m_info_ssp.size(); ++i) {
                        probability_total += m_probability_SSP[i];
                        if (probability_total > rnd) {
                            selected_ssp = i;
                            if (m_info_ssp[i].hills_located.size() == 1) {
                                success_flag = false;
                                break;
                            }
                            if (selected_map.find(i) == selected_map.end()) {
                                selected_map.insert(std::pair<size_t, bool>(i, true));
                                success_flag = true;
                            } else success_flag = false;
                            break;
                        }
                    }
                }
                else if(tries >= m_init_num_ssp) {
                    for (int i = 0; i < m_info_ssp.size(); ++i) {
                        probability_total += m_probability_SSP[i];
                        if (probability_total > rnd) {
                            selected_ssp = i;
                            if (selected_map.find(i) == selected_map.end()) {
                                selected_map.insert(std::pair<size_t, bool>(i, true));
                                success_flag = true;
                            } else success_flag = false;
                            break;
                        }
                    }
                }
                if(success_flag) {
                    const auto &ssp_boundary = m_ssp_tree->getBox(selected_ssp);
                    m_new_indis[selected_count] = std::make_unique<Solution_type>(num_obj, num_con, num_var);
                    for (int i = 0; i < ssp_boundary.size(); ++i) {
                        m_new_indis[selected_count]->variable()[i] = m_random->uniform.nextNonStd(
                                ssp_boundary[i].first, ssp_boundary[i].second);
                    }
                    m_new_indis[selected_count]->initVelocity(m_problem.get(), m_random.get());
                    m_new_indis[selected_count]->initVelocityMax(m_problem.get(), m_random.get());
                    ++selected_count;
                }
            }
            m_swarm_Lips = std::make_unique<SwarmLIP>(Solution_num, m_problem.get(), 500000);
            for (int i = 0; i < m_swarm_Lips->size(); ++i) {
                (*m_swarm_Lips)[i].variable() = m_new_indis[i]->variable();
                (*m_swarm_Lips)[i].velocity() = m_new_indis[i]->velocity();
            }
            rf = m_swarm_Lips->evaluate(m_problem.get(), this);
            if (rf != kNormalEval) return rf;
            m_swarm_Lips->initPbest(m_problem.get());
            m_swarm_Lips->initVelocityMax(m_problem.get(), m_random.get());
        }
        std::map<size_t, bool> located_ssp;
        for (int i = 0; i < m_swarm_Lips->size(); ++i) {
            auto id_ssp = m_ssp_tree->getRegionIdx(((*m_swarm_Lips)[i]).variable().vect());
            if(located_ssp.find(id_ssp) == located_ssp.end()) located_ssp.insert(std::pair<size_t, bool>(id_ssp, true));
        }
        m_origin_cover_rate = located_ssp.size();

        while (!terminating() && (!detectPopDivision(m_swarm_Lips, type))) {
            rf = m_swarm_Lips->evolve(m_problem.get(), this, m_random.get());
            if (rf != kNormalEval) return rf;
        }
        return rf;
    }

    bool DFDL::detectPopDivision(std::unique_ptr<SwarmLIP> & niche, int type) {
        m_new_pop_evolve_count++;
        std::map<size_t, bool> located_ssp;
        size_t now_cover_rate = 0;
        for (int i = 0; i < niche->size(); ++i) {
            archiveSolution((*niche)[i], TaskEval::kExplore);
            auto id_ssp = m_ssp_tree->getRegionIdx(((*niche)[i]).variable().vect());
            if(located_ssp.find(id_ssp) == located_ssp.end()) located_ssp.insert(std::pair<size_t, bool>(id_ssp, true));
        }
        now_cover_rate = located_ssp.size();
//        std::cout << "now_cover_rate: " << now_cover_rate << std::endl;
#ifdef OFEC_DEMO
        updateBuffer();
#endif
        // 等种群收敛到m_divide_factor * 原始占比时，开始聚类；或者演化五代之后直接开始聚类
        if((int)now_cover_rate < m_divide_factor * (int)m_origin_cover_rate || m_new_pop_evolve_count > 5) {
            m_new_pop_evolve_count = 0;
            if(m_clustering_type == ClusteringType::HSLH) {
                std::vector<std::unique_ptr<LI_particle>> indis;
                for (int i = 0; i < niche->size(); ++i) {
                    indis.emplace_back(new LI_particle((*niche)[i]));
                }
                HSLH<LI_particle> cluster(indis, m_problem.get());
                cluster.clustering(min_subpopSize, m_problem.get());
                for (int i(0); i < cluster.size(); ++i) {
                    if (cluster[i].size() > m_min_num_inds) {
                        std::unique_ptr<TypePop> new_pop = std::make_unique<TypePop>(cluster[i].size(), m_problem.get(), DFDLSwarm::PopulationType::Explore);
                        auto iter = cluster[i].begin();
                        for (int j(0); j < cluster[i].size(); ++j, ++iter) {
                            (*new_pop)[j] = *(iter->second);
                        }
                        new_pop->initializeNichingPops(m_problem.get(), m_random.get());
                        m_sub_pop.append(new_pop);
                    }

/*            Real radius = 0.;
            std::unique_ptr<Solution<>> center = std::make_unique<DFDLParticle>(CAST_CONOP(m_problem.get())->numberObjectives(), CAST_CONOP(m_problem.get())->numberConstraints(), CAST_CONOP(m_problem.get())->numberVariables());
            // Computer center
            for (size_t k = 0; k < CAST_CONOP(m_problem.get())->numberVariables(); k++) {
                double x = 0.;
                for (int j(0); j < cluster[i].size(); ++j, ++iter) {
                    auto chr = (iter->second).variable();
                    x = chr[k] + x;
                }
                center->variable()[k] = x / cluster[i].size();
            }
            // Computer radius
            iter = cluster[i].begin();
            for (int j = 0; j < cluster[i].size(); j++, ++iter)
                radius += (iter->second).variableDistance(*center, m_problem.get());
            radius = radius / cluster[i].size();
*/
                }
            }
            else if(m_clustering_type == ClusteringType::NBC) {
                NBC nbc(2.0, NBC::UpdateNBD::kByDistMat, NBC::CalThrshld::kByMean);
                nbc.setData(*niche, m_problem.get());
                nbc.clustering(5);
                size_t pop_size = 0;
                for (const auto &cluster: nbc.clusters()) {
                    if(cluster.size() > m_min_num_inds)
                        pop_size++;
                }
                if(type == InitSolType::Random && pop_size < m_min_init_pop_size) {
                    // 保证有足够的种群进行演化
                    int count = 0;
                    int ref = 0;
                    do {
                        count = 0;
                        nbc.setData(*niche, m_problem.get());
                        nbc.clustering(m_min_init_pop_size + ref);
                        for (const auto &cluster: nbc.clusters()) {
                            if(cluster.size() > m_min_num_inds)
                                count++;
                        }
                        ref++;
                    } while(count < m_min_init_pop_size);
                }
                for (const auto &cluster: nbc.clusters()) {
                    if (cluster.size() > m_min_num_inds) {
                        std::unique_ptr<TypePop> new_pop = std::make_unique<TypePop>(cluster.size(), m_problem.get(),
                                                                                     DFDLSwarm::PopulationType::Explore);
                        for (int i = 0; i < cluster.size(); ++i) {
                            (*new_pop)[i] = niche->at(cluster[i]);
                        }
                        new_pop->initializeNichingPops(m_problem.get(), m_random.get());
                        m_sub_pop.append(new_pop);
                    }
                }
            }
            niche.reset();
            return true;
        }
        else return false;
    }

    int DFDL::createExploreSwarm() {
        updateNumSolutions();
        int rf = createSubSwarm(m_num_increase);
        return rf;
    }

    int DFDL::createSupplementarySwarm(int pop_size) {
        int rf = kNormalEval;
        for (int i = 0; i < m_num_hills; ++i) {
            if(!isAnySolutionInHill(i)) {
//                Real volume = 0.;
//                for (const auto &index: m_unique_ssp_hills[i]) {
//                    volume += m_info_ssp[index].volume;
//                }
                // 如果吸引域的占比体积太小，考虑在该吸引域里采一次样
                if (m_unique_ssp_hills[i].size() == 1) {
                    const auto &ssp_boundary = m_ssp_tree->getBox(m_unique_ssp_hills[i][0]);
                    auto *sample_sol = new Solution<>(CAST_CONOP(m_problem.get())->numberObjectives(),
                                                      CAST_CONOP(m_problem.get())->numberConstraints(),
                                                      CAST_CONOP(m_problem.get())->numberVariables());
                    for (int k = 0; k < CAST_CONOP(m_problem.get())->numberVariables(); ++k) {
                        sample_sol->variable()[k] = m_random->uniform.nextNonStd(ssp_boundary[k].first,
                                                                                         ssp_boundary[k].second);
                    }
                    rf = sample_sol->evaluate(m_problem.get(), this, true);
                    if(rf != kNormalEval) break;
                    archiveSolution(sample_sol, TaskEval::kExplore);
                }
            }
        }
        // 采样后再进行吸引域的聚类
        groupSubspace();
        for (int i = 0; i < m_num_hills; ++i) {
            if(isAllow2SupplySwarm(i)) {
                /* 在该满足条件的吸引域中补充探索种群 */
                std::unique_ptr<TypePop> new_pop = std::make_unique<TypePop>(pop_size, m_problem.get(), DFDLSwarm::PopulationType::Explore);
                for (int j = 0; j < pop_size; ++j) {
                    auto index = m_random->uniform.nextNonStd(0, (int) m_unique_ssp_hills[i].size());
                    const auto &ssp_boundary = m_ssp_tree->getBox(m_unique_ssp_hills[i][index]);
                    for (int k = 0; k < ssp_boundary.size(); ++k) {
                        (*new_pop)[j].variable()[k] = m_random->uniform.nextNonStd(ssp_boundary[k].first, ssp_boundary[k].second);
                    }
                    (*new_pop)[j].initVelocity(m_problem.get(), m_random.get());
                    (*new_pop)[j].initVelocityMax(m_problem.get(), m_random.get());
                    rf = (*new_pop)[j].evaluate(m_problem.get(), this, true);
                    if (rf != kNormalEval) return rf;
                    archiveSolution((*new_pop)[j], TaskEval::kExplore);
                }
                new_pop->initializeParameters(m_problem.get(), m_random.get());
                m_sub_pop.append(new_pop);
            }
        }
        return rf;
    }

    int DFDL::createExploitSwarm(int pop_size) {
        int rf = kNormalEval;
        if(m_num_hills == 1) return rf;
        for (int i = 0; i < m_num_hills; ++i) {
            Real hill_volume = 0.;
            for (int j = 0; j < m_unique_ssp_hills[i].size(); ++j) {
                hill_volume += m_info_ssp[m_unique_ssp_hills[i][j]].volume;
            }
//            std::cout << "hill_volume: " << hill_volume << std::endl;
            // 在体积小于1%的吸引域中不产生开发群体
            if(hill_volume <= 1) continue;
            // 依据不确定度计算吸引域中子空间产生解的概率
            std::vector<int> index;
            std::vector<Real> prob;
            std::map<int, Real> pro_ssp;
            Real temp = 0.;
            for (int j = 0; j < m_unique_ssp_hills[i].size(); ++j) {
                temp += m_info_ssp[m_unique_ssp_hills[i][j]].uncertainty;
            }
            for (int j = 0; j < m_unique_ssp_hills[i].size(); ++j) {
                index.push_back(m_unique_ssp_hills[i][j]);
                prob.push_back(m_info_ssp[m_unique_ssp_hills[i][j]].uncertainty/temp);
            }

            std::vector<int> selected;
            std::vector<int> selected_ssp;
            selected.resize(pop_size);
            std::vector<Real> cumulative_probabilities(prob.size());
            // 累计概率分布
            std::partial_sum(prob.begin(), prob.end(),cumulative_probabilities.begin());
            // 随机抽样
            std::random_device rd;
            std::mt19937 gen(rd());
            std::discrete_distribution<int> dist(cumulative_probabilities.begin(), cumulative_probabilities.end());
            std::generate(selected.begin(), selected.end(), [&] { return dist(gen); });
            for (const auto &item: selected) {
                selected_ssp.push_back(index[item]);
            }
            // 产生pop_size大小的种群
            int count = 0;
            std::unique_ptr<TypePop> new_pop = std::make_unique<TypePop>(pop_size, m_problem.get(), DFDLSwarm::PopulationType::Exploit);
            for (const auto &ssp: selected_ssp) {
                const auto &ssp_boundary = m_ssp_tree->getBox(ssp);
                for (int k = 0; k < ssp_boundary.size(); ++k) {
                    (*new_pop)[count].variable()[k] = m_random->uniform.nextNonStd(ssp_boundary[k].first, ssp_boundary[k].second);
                }
                (*new_pop)[count].initVelocity(m_problem.get(), m_random.get());
                (*new_pop)[count].initVelocityMax(m_problem.get(), m_random.get());
                rf = (*new_pop)[count].evaluate(m_problem.get(), this, true);
                if (rf != kNormalEval) return rf;
                archiveSolution((*new_pop)[count], TaskEval::kExplore);
                count++;
            }
            new_pop->initializeParameters(m_problem.get(), m_random.get());
            m_sub_pop.append(new_pop);
        }
        return rf;
    }

    void DFDL::awakeSwarms() {
        for (size_t i = 0; i < m_sub_pop.size(); i++) {
            if(!m_sub_pop[i].isActive()) {
                m_sub_pop[i].setActive(true);
                m_sub_pop[i][m_sub_pop[i].bestIndex(m_problem.get())].brownianMove(m_problem.get(), m_random.get(), m_sub_pop[i].getCurRadius() * 10);
                m_sub_pop[i].setConvergeFlag(false);
                m_sub_pop[i].setStagnantedFlag(false);
                m_sub_pop[i].m_stagnanted_count = 0;
            }
        }
    }

    void DFDL::run_() {
        int rf = kNormalEval;
        initSubSwarm();
        int count = 0;
        while (!terminating()) {
            rf = popEvolve();
            if (rf != kNormalEval) {
                handleChange(rf);
            }
            checkOverCrowdAndConverge();
            rf = checkStagnantion();
            if (rf != kNormalEval) {
                handleChange(rf);
            }
            groupSubspace();
            // 吸引域聚类之后对新出现的吸引域 立刻补充种群
            createSupplementarySwarm(4);
            avoidRepeatSearch();
            updatePotentialExplore();
            if(checkDiversity(m_problem.get())) {
                //++count;
                //std::cout << "activate diversity method " << count << std::endl;
                //std::cout << "at evaluate time:  " << this->m_evaluations << std::endl;
                awakeSwarms();
                updateHills();
                updateUncertaintySSP();
                calculateProbabilitySSP();
                removeOutdatedHisSolutions(m_outdated_evaluation);
                rf = createExploreSwarm();
                if (rf != kNormalEval) {
                    handleChange(rf);
                }
                rf = createExploitSwarm(4);
                if (rf != kNormalEval) {
                    handleChange(rf);
                }
            }

#ifdef OFEC_DEMO
            updateBuffer();
#endif
        }
    }

    void DFDL::handleChange(int rf) {
        if (rf & kChangeCurEval) {
            measureMultiPops(true);
            m_archive_best.clear();
            // 更新子空间的最佳适应度值
            // for (int i(0); i < m_info_ssp.size(); ++i) {
            //     m_info_ssp[i].best_sol
            // }
        }
    }

    void DFDL::measureMultiPops(bool flag) {
        for (int i(0); i < m_sub_pop.size(); ++i) {
            for (int j(0); j < m_sub_pop[i].size(); j++) {
                m_sub_pop[i][j].pbest().evaluate(m_problem.get(), this, flag);
            }
        }
    }

    void DFDL::archiveEvolveSolution() {
        for (size_t i = 0; i < m_sub_pop.size(); i++) {
            if(!m_sub_pop[i].isActive()) continue;
            for (size_t j = 0; j < m_sub_pop[i].size(); j++) {
                archiveSolution(m_sub_pop[i][j], TaskEval::kExplore);
            }
        }
    }

    void DFDL::updateUncertaintySSP() {
        switch (m_uncertainty_type) {
            case UncertaintyType::linear :
                calculateUncertaintySSPbyLinear();
                break;
            case UncertaintyType::exponential :
                calculateUncertaintySSPbyExponential();
                break;
        }
    }

    void DFDL::updateInfoSSP(const SolutionType *p_sol, TaskEval task) {
        auto id_ssp = m_ssp_tree->getRegionIdx(p_sol->variable().vect());
        if (!m_info_ssp[id_ssp].best_sol || p_sol->dominate(*m_info_ssp[id_ssp].best_sol, m_problem.get())) {
            m_info_ssp[id_ssp].best_sol = p_sol;
            m_info_ssp[id_ssp].best_in_exploit = task == TaskEval::kExploit;
        }
        if(m_info_ssp[id_ssp].his_sparse_sols.empty()) m_info_ssp[id_ssp].his_sparse_sols.push_back(p_sol);
        else {
            bool insert_flag = true;
            for (int i = 0; i < m_info_ssp[id_ssp].his_sparse_sols.size(); ++i) {
                if (p_sol->variableDistance(*m_info_ssp[id_ssp].his_sparse_sols[i], m_problem.get()) < converge_radius) {
                    if(p_sol->dominate(*m_info_ssp[id_ssp].his_sparse_sols[i], m_problem.get())) {
                        m_info_ssp[id_ssp].his_sparse_sols[i] = p_sol;
                        insert_flag = false;
                        break;
                    }
                }
            }
            if(insert_flag) m_info_ssp[id_ssp].his_sparse_sols.push_back(p_sol);
        }
        m_info_ssp[id_ssp].his_sols.push_back(p_sol);
        ++m_info_ssp[id_ssp].num_sample;
        if (task == TaskEval::kExplore)
            m_info_ssp[id_ssp].his_explore.push_back(p_sol);
        else
            m_info_ssp[id_ssp].his_exploit.push_back(p_sol);
    }

    void DFDL::removeOutdatedHisSolutions(size_t outdated_evaluation) {
        for (int i = m_his_sols.size() - 1; i >= 0; --i) {
            auto iter = m_his_sols.begin();
            std::advance(iter, i);
            // 如果这个解的评价时间大于约定的过期期限，则将该解直接清除
            if((GET_ALG(this)->evaluations() - (*iter)->timeEvaluate()) > outdated_evaluation)
                m_his_sols.erase(iter);
        }
//        std::cout << "----------------------------number of his_solutions: " << m_his_sols.size() << std::endl;
        for (auto &ssp: m_info_ssp) {
            // his_sols
            for (int i = ssp.his_sols.size() - 1; i >= 0; --i) {
                // 如果这个解的评价时间大于约定的过期期限，则将该解直接清除
                if((GET_ALG(this)->evaluations() - ssp.his_sols[i]->timeEvaluate()) > outdated_evaluation)
                    ssp.his_sols.erase(ssp.his_sols.begin() + (int)i);
            }
            // his_sparse_sols
            for (int i = ssp.his_sparse_sols.size() - 1; i >= 0; --i) {
                // 如果这个解的评价时间大于约定的过期期限，则将该解直接清除
                if((GET_ALG(this)->evaluations() - ssp.his_sparse_sols[i]->timeEvaluate()) > outdated_evaluation)
                    ssp.his_sparse_sols.erase(ssp.his_sparse_sols.begin() + (int)i);
            }
            // his_exploit
            for (int i = ssp.his_exploit.size() - 1; i >= 0; --i) {
                auto iter = ssp.his_exploit.begin();
                std::advance(iter, i);
                // 如果这个解的评价时间大于约定的过期期限，则将该解直接清除
                if((GET_ALG(this)->evaluations() - (*iter)->timeEvaluate()) > outdated_evaluation)
                    ssp.his_exploit.erase(iter);
            }
            // his_explore
            for (int i = ssp.his_explore.size() - 1; i >= 0; --i) {
                // 如果这个解的评价时间大于约定的过期期限，则将该解直接清除
                if((GET_ALG(this)->evaluations() - ssp.his_explore[i]->timeEvaluate()) > outdated_evaluation)
                    ssp.his_explore.erase(ssp.his_explore.begin() + (int)i);
            }
        }

        // update best_sol in each ssp
        for (auto &ssp: m_info_ssp) {
            if(ssp.his_sols.empty())
                ssp.best_sol = nullptr;
            if (!ssp.his_sols.empty()) {
                const SolutionType *best_sol = ssp.his_sols.front();
                for (auto &sol: ssp.his_sols) {
                    if (sol != best_sol && sol->dominate(*best_sol, m_problem.get()))
                        best_sol = sol;
                }
                ssp.best_sol = best_sol;
            }
        }
    }

    void DFDL::distinctHisSols() {
        std::vector<std::vector<const SolutionType*>> his_sols_without_repeat;
        his_sols_without_repeat.resize(m_info_ssp.size());
        for (int i = 0; i < m_info_ssp.size(); ++i) {
            for (const auto &sol: m_info_ssp[i].his_sols) {
                for (const auto &item: his_sols_without_repeat[i]) {
                    if(sol->variableDistance(*item, m_problem.get()) < converge_radius && sol->objectiveDistance(*item) < converge_radius) continue;
                    else his_sols_without_repeat[i].push_back(sol);
                }
            }
        }
    }
    void DFDL::calculateUncertaintySSPbyLinear() {
        std::vector<Real> ssp_frequency;
        std::vector<Real> ssp_rough;
        std::vector<Real> ssp_coverage;
        std::vector<Real> ssp_fitness;
        ssp_frequency.resize(m_info_ssp.size());
        ssp_coverage.resize(m_info_ssp.size());
        ssp_rough.resize(m_info_ssp.size());
        ssp_fitness.resize(m_info_ssp.size());
//        std::vector<std::vector<const SolutionType*>> his_sols_without_repeat;
//        his_sols_without_repeat.resize(m_info_ssp.size());
//        for (int i = 0; i < m_info_ssp.size(); ++i) {
//            去除重复的记录
//            for (const auto &sol: m_info_ssp[i].his_sols) {
//                bool sol_flag = true;
//                if(his_sols_without_repeat[i].empty()) his_sols_without_repeat[i].push_back(sol);
//                for (const auto &item: his_sols_without_repeat[i]) {
//                    if(sol->variableDistance(*item, m_problem.get()) < converge_radius) {
//                        sol_flag = false;
//                        break;
//                    }
//                }
//                if(sol_flag) his_sols_without_repeat[i].push_back(sol);
//            }
//            std::cout << "before:" << m_info_ssp[i].his_sols.size() << '\t' << "after:" <<  m_info_ssp[i].his_sparse_sols.size() << std::endl;
//        }
        // Sample
        for (int i = 0; i < m_info_ssp.size(); ++i) {
            ssp_frequency[i] = m_info_ssp[i].his_sparse_sols.size() / m_info_ssp[i].volume;
        }
        Real max_sample_rate = *std::max_element(ssp_frequency.begin(), ssp_frequency.end());
        Real min_sample_rate = *std::min_element(ssp_frequency.begin(), ssp_frequency.end());
        // Coverage
//        std::cout << "---------------------------------" << std::endl;
        for (int i = 0; i < m_info_ssp.size(); ++i) {
            auto &box = m_ssp_tree->getBox(i);
            ssp_tree_for_coverage.reset(new nanoflann::PartitioningKDTree<Real>(m_problem->numberVariables(), 2));
            ssp_tree_for_coverage->setInitBox(box);
            size_t initial_num = 100;
            // 将子空间划分成100份
            ssp_tree_for_coverage->inputRatioData(std::vector<Real>(100, 1. / 100));
            ssp_tree_for_coverage->buildIndex();
            std::map<size_t, bool> located_ssp;
            for (auto & his_sol : m_info_ssp[i].his_sols) {
                auto id_ssp = ssp_tree_for_coverage->getRegionIdx(his_sol->variable().vect());
                if(located_ssp.find(id_ssp) == located_ssp.end()) located_ssp.insert(std::pair<size_t, bool>(id_ssp, true));
            }
            ssp_coverage[i] = located_ssp.size();
//            std::cout << "coverage " << i << ":" << ssp_coverage[i] << std::endl;
        }
        // Fitness
        for (int i = 0; i < m_info_ssp.size(); ++i) {
            ssp_fitness[i] = m_info_ssp[i].obj;
        }
        getBestAndWorstSSPIndex();
        Real bestSSPObj = m_info_ssp[m_best_ssp_index].obj;
        Real worstSSPObj = m_info_ssp[m_worst_ssp_index].obj;

        // Time： 离得越远的评价越没有用
        for (int i = 0; i < m_info_ssp.size(); ++i) {
            m_info_ssp[i].uncertainty =
                    m_alpha * ((0.8 / (min_sample_rate - max_sample_rate)) * ssp_frequency[i] + (1 - 0.8*min_sample_rate/(min_sample_rate - max_sample_rate)))
                    + m_beta * (100 - ssp_coverage[i]) * 0.01
                    + m_gama * (ssp_fitness[i] - worstSSPObj)/(bestSSPObj - worstSSPObj);
        }
        // Rugged
/*        for (int i = 0; i < m_info_ssp.size(); ++i) {
            Real rough = 0.0;
            if(m_info_ssp[i].his_sols.empty() || m_info_ssp[i].his_sols.size() == 1) ssp_rough[i] = 0.0;
            else {
                for (int j = 0; j < m_info_ssp[i].his_sols.size() - 1; j++) {
                    for (int k = j + 1; k < m_info_ssp[i].his_sols.size(); k++) {
                        Real rough_temp = 0.0;
                        Real d_dis = m_info_ssp[i].his_sols[j]->variableDistance(*m_info_ssp[i].his_sols[k], m_problem.get());
                        Real d_obj = m_info_ssp[i].his_sols[j]->objective(0) - m_info_ssp[i].his_sols[k]->objective(0);
                        if (d_dis < 10e-3 && d_obj < 10e-3) rough_temp = 0.0;
                        else rough_temp = d_obj / d_dis;
                        if (rough_temp > rough) rough = rough_temp;
                        if(rough_temp>10) std::cout << "d_obj" << d_obj << "\t" << "d_dis" << d_dis << std::endl;
                    }
                }
                ssp_rough[i] = rough;
            }
        }
        for (int i = 0; i < m_info_ssp.size(); ++i) {
            m_info_ssp[i].uncertainty = ssp_rough[i];
        }*/

    }

    void DFDL::calculateUncertaintySSPbyExponential() {
        getBestAndWorstSSPIndex();
        Real bestSSPObj = m_info_ssp[m_best_ssp_index].obj;
        Real worstSSPObj = m_info_ssp[m_worst_ssp_index].obj;
        size_t curEvals = this->evaluations();
        for (auto &ssp: m_info_ssp) {
            Real reduceFactor;
            if (!ssp.his_sols.empty()) {
                reduceFactor = 1;
            } else {
                reduceFactor = 0.8;
            }
            int evalTime = 0;
            Real obsolete;
            if(!ssp.his_sols.empty()) {
                for (const auto &sol: ssp.his_sols) {
                    evalTime += sol->timeEvaluate();
                }
                obsolete = evalTime * 1.0 / ssp.his_sols.size();
            }
            else {
                obsolete = 1.0;
            }
            ssp.uncertainty = m_alpha * (2 / (1 + exp(-ssp.volume * 1.0 * CAST_CONOP(m_problem.get())->numberVariables()/ ssp.his_sols.size())) - 1) +
                              m_beta * (ssp.obj - worstSSPObj)/(bestSSPObj - worstSSPObj) * reduceFactor;
//                    m_gama * obsolete/ m_obsoleteFactor;
        }
    }

    int DFDL::popEvolve() {
        int rf = m_sub_pop.evolve(m_problem.get(), this, m_random.get());
        if (rf != kNormalEval) return rf;
        archiveEvolveSolution();
        return rf;
    }

    void DFDL::updateMemory() {
        if (m_num_hills >= m_pop_indis.size()) {
            unsigned size = m_pop_indis.size();
            m_pop_indis.resize(m_num_hills + 1);
            for (unsigned i = size; i < m_pop_indis.size(); i++) {
                m_pop_indis[i].push_back(0);
                m_pop_indis[i].push_back(0);
            }
        }
        m_pop_indis[m_num_hills].push_back(m_total_indi_size);
        double mean = m_pop_indis[m_num_hills][0];
        m_pop_indis[m_num_hills][0] = (mean * (m_pop_indis[m_num_hills].size() - 3) + m_total_indi_size) / (m_pop_indis[m_num_hills].size() - 2);
        m_pop_indis[m_num_hills][1] = 0;
        for (unsigned k = 2; k < m_pop_indis[m_num_hills].size(); k++) {
            m_pop_indis[m_num_hills][1] += (m_pop_indis[m_num_hills][k] - m_pop_indis[m_num_hills][0]) * (m_pop_indis[m_num_hills][k] - m_pop_indis[m_num_hills][0]);
        }
        m_pop_indis[m_num_hills][1] = sqrt(m_pop_indis[m_num_hills][1] / (m_pop_indis[m_num_hills].size() - 2));
    }

    void DFDL::updateNumSolutions() {
        m_cur_hill_size = m_num_hills;
        m_total_indi_size = 0;
        for (auto & pop : m_sub_pop) {
            m_total_indi_size += pop->size();
        }

        if (m_first_update) {
            m_pre_hill_size = m_cur_hill_size;
            m_pre_indi_size = m_total_indi_size;
            m_first_update = false;
        }
        updateMemory();
        int indisNum = static_cast<int>(m_random->normal.nextNonStd(m_pop_indis[m_cur_hill_size][0], m_pop_indis[m_cur_hill_size][1]));
        size_t dif = m_cur_hill_size - m_pre_hill_size;
        double ratio = fabs(dif * 1. / mc_off_peak);
        if (ratio >= 1) {
            m_next_indi_size = indisNum + mc_step_indis * dif;
        }
        else if (ratio == 0) {
            m_next_indi_size = indisNum;
        }
        else {
            double p = m_random->uniform.next();
            if (p <= ratio) m_next_indi_size = indisNum + mc_step_indis * dif;
            else m_next_indi_size = indisNum;
        }

        m_pre_hill_size = m_cur_hill_size;

        if (m_next_indi_size <= m_total_indi_size + mc_step_indis * 2) {
            m_next_indi_size = m_total_indi_size + mc_step_indis * 2;
        }

        if (m_next_indi_size > mc_max_idxs) m_next_indi_size = mc_max_idxs;
        if (m_next_indi_size < mc_min_idxs) m_next_indi_size = mc_min_idxs;
        m_pre_indi_size = m_next_indi_size;
//        std::cout << "m_next_indi_size: " << m_next_indi_size << std::endl;
        m_num_increase = m_next_indi_size - m_total_indi_size;
        if (m_num_increase >= mc_step_indis * 6) m_num_increase = mc_step_indis * 6;
    }

    void DFDL::checkOverCrowdAndConverge() {
        for (int i = 0; i < m_sub_pop.size(); ++i) {
            m_sub_pop[i].checkOverCrowd(m_problem.get(), m_random.get());
            m_sub_pop[i].checkConverged(m_problem.get());
            if(m_sub_pop[i].getConvergeFlag()){
                 m_sub_pop[i].setActive(false);
            }
            if(m_sub_pop[i].getSeedFlag()) {
                // archive best_sols
                // 先查询archive里有没有这个解？
                bool save_flag = true;
                for (const auto &archive_best: m_archive_best) {
                    Real dis = m_sub_pop[i][m_sub_pop[i].bestIndex(m_problem.get())].variableDistance(*archive_best, m_problem.get());
                    if(dis < seed_radius) {
                        save_flag = false;
                        break;
                    }
                }
                if(save_flag) {
                    m_archive_best.emplace_back(new Solution<>(m_sub_pop[i][m_sub_pop[i].bestIndex(m_problem.get())]));
//                    std::cout << "m_archive_best size : " << m_archive_best.size() << std::endl;
//                    updateHills();
                }
            }
        }
    }

    int DFDL::checkStagnantion() {
        int rf = kNormalEval;
        Real avgR = 0.;
        if(m_sub_pop.size() > 0) {
            for (int i = 0; i < m_sub_pop.size(); ++i) {
                m_sub_pop[i].updateCurRadius(m_problem.get());
                avgR += m_sub_pop[i].getCurRadius();
            }
            avgR = avgR / m_sub_pop.size();
            for (int i = 0; i < m_sub_pop.size(); ++i) {
                m_sub_pop[i].checkStagnanted(avgR);
                if(m_sub_pop[i].getStagnantedFlag()) {
                    rf = m_sub_pop[i][m_sub_pop[i].bestIndex(m_problem.get())].cauchyMove(m_problem.get(), this, m_random.get());
                    if (rf != kNormalEval) break;
//                    m_sub_pop[i].setStagnantedFlag(false);
//                    m_sub_pop[i].m_stagnanted_count = 0;
                }
            }
        }
        return rf;
    }

    int DFDL::checkOverlapping() {
        for (size_t i = 0; i < m_sub_pop.size(); i++) {
            if (this->m_sub_pop[i].size() == 0) continue;
            for (size_t j = i + 1; j < m_sub_pop.size(); j++) {
                if (this->m_sub_pop[j].size() == 0) continue;
                const size_t best_i_inx = m_sub_pop[i].bestIndex(m_problem.get());
                const size_t best_j_inx = m_sub_pop[j].bestIndex(m_problem.get());
                Real dist = m_problem->variableDistance(m_sub_pop[i][best_i_inx], m_sub_pop[j][best_j_inx]);
                if (dist < m_sub_pop[i].getInitRadius() && dist < m_sub_pop[j].getInitRadius()) {
                    int c1(0);
                    int c2(0);
                    for (size_t k = 0; k < m_sub_pop[j].size(); k++) {
                        dist = m_problem->variableDistance(m_sub_pop[i][best_i_inx], m_sub_pop[j][k]);
                        if (dist < m_sub_pop[i].getInitRadius()) c1++;
                    }

                    for (size_t k = 0; k < m_sub_pop[i].size(); k++) {
                        dist = m_problem->variableDistance(m_sub_pop[j][best_j_inx], m_sub_pop[i][k]);
                        if (dist < m_sub_pop[j].getInitRadius()) c2++;	//i->j Fixed in 24/10/2021
                    }
                    if(c1 > 0 && c2 > 0){
                        //i->j j->i Fixed in 24/10/2021
                        int idx = -1;
                        if (m_sub_pop[i][best_i_inx].pbest().dominate(m_sub_pop[j][best_j_inx].pbest(), m_problem.get())) {
                            m_sub_pop[i].merge(m_sub_pop[j],m_problem.get(), m_random.get());
                            m_sub_pop[i].checkOverCrowd(m_problem.get(), m_random.get());
                            m_sub_pop.remove(m_sub_pop.begin() + j);
                            idx = j;
                        }
                        else{
                            m_sub_pop[j].merge(m_sub_pop[i], m_problem.get(), m_random.get());
                            m_sub_pop[j].checkOverCrowd(m_problem.get(), m_random.get());
                            m_sub_pop.remove(m_sub_pop.begin() + i);
                            idx = i;
                        }
                        return idx;
                    }
                }
            }
        }
        return -1;
    }

    void DFDL::hillShielding() {

        /* 1. 排斥掉进入同一个吸引域的探索中群 */
        // 更新每个种群落入哪个吸引域
        updateSwarmHillLocated();
        // 更新每个独立吸引域中包含那些种群
        updateUniqueHillsSwarm();
        // 如果独立吸引域中包含超过一个探索种群，将会移除较差的那些探索种群
        std::vector<int> need_to_remove;
        for (auto &hill_pops: m_hills_pops) {
            if(hill_pops.size() > 1) {
                int n = hill_pops.size();
                for (int i = 0; i < n-1; i++) {
                    int max_idx = i;
                    for (int j = i+1; j < n; j++) {
                        if (m_sub_pop[hill_pops[j]][m_sub_pop[hill_pops[j]].bestIndex(m_problem.get())].dominate(
                                m_sub_pop[hill_pops[max_idx]][m_sub_pop[hill_pops[max_idx]].bestIndex(m_problem.get())], m_problem.get())) {
                            max_idx = j;
                        }
                    }
                    std::swap(hill_pops[i], hill_pops[max_idx]);
                }
                need_to_remove.insert(need_to_remove.end(), hill_pops.begin() + 1, hill_pops.end());
            }
        }
        std::sort(need_to_remove.begin(), need_to_remove.end(), [](int a, int b) { return a > b; });
        for (const auto &remove_index: need_to_remove) {
            if (m_sub_pop[remove_index].m_pop_type == DFDLSwarm::PopulationType::Explore)
                m_sub_pop.remove(m_sub_pop.begin() + remove_index);
        }

        /* 2. 排斥掉进入吸引域中种群排斥半径里的种群 */
        // 更新每个种群落入哪个吸引域
        updateSwarmHillLocated();
        // 更新每个独立吸引域中包含那些种群
        updateUniqueHillsSwarm();
        // 如果独立吸引域中包含开发种群，则让其靠近该吸引域中最好的种群的排斥半径中排斥掉
        std::vector<size_t> need_to_exclusion;
        int num_hill = 0;
        for (auto &hill_pops: m_hills_pops) {
            // 计算该吸引域的体积所占的大小
            Real volume = 0.;
            for (const auto &index: m_unique_ssp_hills[num_hill]) {
                volume += m_info_ssp[index].volume;
            }
            if(volume < m_min_active_hill_volume) volume = m_min_active_hill_volume;
            // 确定排斥半径 [0.5 * X/(p^(1/d))] * volume%
            Real radius = volume * 0.01 * 1 * (CAST_CONOP(m_problem.get())->range(0).second - CAST_CONOP(m_problem.get())->range(0).first)/pow((double)m_num_hills, 1./CAST_CONOP(m_problem.get())->numberVariables());
//            std::cout << "exclusion radius of hill " << num_hill << " is: " << radius << std::endl;
            // 当种群中的最好个体靠近吸引域中最好粒子的排斥半径内被排除
            if(hill_pops.size() > 0) {
                for (int i = 0; i < hill_pops.size(); ++i) {
                    for (int j = 0; j < m_sub_pop.size(); ++j) {
                        if(hill_pops[i] == j) continue;
                        if(m_sub_pop[hill_pops[i]][m_sub_pop[hill_pops[i]].bestIndex(m_problem.get())].variableDistance(m_sub_pop[j][m_sub_pop[j].bestIndex(m_problem.get())],m_problem.get()) <= radius) {
                            if (m_sub_pop[hill_pops[i]][m_sub_pop[hill_pops[i]].bestIndex(m_problem.get())].dominate(
                                    m_sub_pop[j][m_sub_pop[j].bestIndex(m_problem.get())], m_problem.get())) {
                                need_to_exclusion.push_back(j);
                            }
                        }
                    }
                }
            }
            num_hill++;
        }

        /* 3. 排斥掉互相进入排斥半径的较差种群 */
        for (int i = 0; i < m_sub_pop.size(); ++i) {
            for (int j = 0; j < m_sub_pop.size(); ++j) {
                if(i == j) continue;
                // 确定排斥半径 [0.5 * X/(p^(1/d))] * 1 %
                Real radius = 0.01 * m_min_active_hill_volume * 1 * (CAST_CONOP(m_problem.get())->range(0).second - CAST_CONOP(m_problem.get())->range(0).first)/pow((double)m_num_hills, 1./CAST_CONOP(m_problem.get())->numberVariables());
                if(m_sub_pop[i][m_sub_pop[i].bestIndex(m_problem.get())].variableDistance(m_sub_pop[j][m_sub_pop[j].bestIndex(m_problem.get())],m_problem.get()) <= radius) {
                    if (m_sub_pop[i][m_sub_pop[i].bestIndex(m_problem.get())].dominate(
                            m_sub_pop[j][m_sub_pop[j].bestIndex(m_problem.get())], m_problem.get())) {
                        need_to_exclusion.push_back(j);
                    }
                }
            }
        }


        std::sort(need_to_exclusion.begin(), need_to_exclusion.end());
        auto last = std::unique(need_to_exclusion.begin(), need_to_exclusion.end());
        need_to_exclusion.erase(last, need_to_exclusion.end());
        std::sort(need_to_exclusion.rbegin(), need_to_exclusion.rend());
//        std::cout << "Need To Exclusion" << need_to_exclusion.size() << std::endl;
        for (const auto &exclude_index: need_to_exclusion) {
            if (m_sub_pop[exclude_index].m_pop_type == DFDLSwarm::PopulationType::Exploit)
                m_sub_pop.remove(m_sub_pop.begin() + exclude_index);
        }
    }

    void DFDL::updateSwarmHillLocated() {
        for (size_t i = 0; i < m_sub_pop.size(); i++) {
            std::map<int, bool> located;
            for (size_t j = 0; j < m_sub_pop[i].size(); j++) {
                int idx = SolutionLocatedInUniqueHills(m_sub_pop[i][j]);
                if(located.find(idx) == located.end()) {
                    located.insert(std::pair<int, bool>(idx, true));
                }
            }
            if(located.size() == 1 && located.find(-1) == located.end()) {
                m_sub_pop[i].setHillLocated(located.begin()->first);
            }
            else m_sub_pop[i].setHillLocated(-1);
        }
    }

    void DFDL::updateUniqueHillsSwarm() {
        m_hills_pops.clear();
        m_hills_pops.resize(m_num_hills);
        for (size_t i = 0; i < m_sub_pop.size(); i++) {
            if(m_sub_pop[i].getHillLocated() != -1) {
                m_hills_pops[m_sub_pop[i].getHillLocated()].emplace_back(i);
            }
        }
    }

    bool DFDL::isAllow2SupplySwarm(int hill) {
        Real volume = 0.;
        for (const auto &index: m_unique_ssp_hills[hill]) {
            volume += m_info_ssp[index].volume;
        }
        if(volume < m_min_active_hill_volume)
            return false;
        return !isAnySolutionInHill(hill);
    }

    bool DFDL::isAnySolutionInHill(int hill) {
        for (size_t i = 0; i < m_sub_pop.size(); i++) {
            for (size_t j = 0; j < m_sub_pop[i].size(); j++) {
                auto id_ssp = m_ssp_tree->getRegionIdx(m_sub_pop[i][j].variable().vect());
                if(std::find(m_unique_ssp_hills[hill].begin(), m_unique_ssp_hills[hill].end(), id_ssp) != m_unique_ssp_hills[hill].end()) {
                    return true;
                }
            }
        }
        return false;
    }

    void DFDL::avoidRepeatSearch() {
        if(m_num_hills == 1)
            while (checkOverlapping() != -1);
        else
            hillShielding();
    }

    void DFDL::getBestAndWorstSSPIndex() {
        m_best_ssp_index = 0;
        m_worst_ssp_index = 0;
        if (m_problem->optimizeMode()[0] == OptimizeMode::kMaximize) {
            for (int i = 0; i < m_info_ssp.size(); ++i) {
                if (m_info_ssp.at(i).obj > m_info_ssp[m_best_ssp_index].obj)
                    m_best_ssp_index = i;
                if(m_info_ssp.at(i).obj < m_info_ssp[m_worst_ssp_index].obj)
                    m_worst_ssp_index = i;
            }
        }
        else if(m_problem->optimizeMode()[0] == OptimizeMode::kMinimize) {
            for (int i = 0; i < m_info_ssp.size(); ++i) {
                if (m_info_ssp.at(i).obj < m_info_ssp[m_best_ssp_index].obj)
                    m_best_ssp_index = i;
                if (m_info_ssp.at(i).obj > m_info_ssp[m_worst_ssp_index].obj)
                    m_worst_ssp_index = i;
            }
        }
    }

    void DFDL::updateHills() {
        if (!m_his_sols.empty()) {
            std::vector<float> input_size;
            for (size_t id_hill = 0; id_hill < m_num_hills; ++id_hill) {
                std::vector<const SolutionType*> his_sols;
                getHisSolsInHill(his_sols, id_hill);
                std::vector<const SolutionType*> seed_sols;
                for (auto &pop: m_sub_pop) {
                    pop->checkConverged(m_problem.get());
                    if(pop->getSeedFlag()) {
                        seed_sols.push_back(&pop->at(pop->bestIndex(m_problem.get())));
                    }
                }
//                for (const auto &archive_best: m_archive_best) {
//                    seed_sols.push_back(archive_best.get());
//                }
                subdivideSpace(seed_sols);
            }
        }
        groupSubspace();
        updatePotentialExplore();
    }

    void DFDL::initSubSwarm() {
        if(m_pop_start_type == PopStart::MultiPopulation)
            createNewSwarms(m_size_pop, InitSolType::Random);
        else if(m_pop_start_type == PopStart::Niching_Lips)
            createNichingSwarm(m_size_pop, InitSolType::Random);
    }

    int DFDL::createSubSwarm(int num) {
        int rf = kNormalEval;
        if(m_pop_start_type == PopStart::MultiPopulation)
            rf = createNewSwarms(num, InitSolType::AtHighUncertainty);
        else if(m_pop_start_type == PopStart::Niching_Lips)
            rf = createNichingSwarm(num, InitSolType::AtHighUncertainty);
        return rf;
    }

    void DFDL::groupSubspace() {
        HGHEC::assignSubspaceFitness();
        HGHEC::groupSubspace();
        m_unique_ssp_hills = getUniqueSSPinHills(m_ssp_hills);
    }

    int DFDL::SolutionLocatedInUniqueHills(const Solution_type & Solution) {
        auto id_ssp = m_ssp_tree->getRegionIdx(Solution.variable().vect());
        for (int i = 0; i < m_unique_ssp_hills.size(); i++) {
            for (int j = 0; j < m_unique_ssp_hills[i].size(); j++) {
                if (m_unique_ssp_hills[i][j] == id_ssp) {
                    return i;
                }
            }
        }
        return -1;
    }

    std::vector<std::vector<size_t>> DFDL::getUniqueSSPinHills(const std::vector<std::vector<size_t>> & arrays) {
        using std::vector;
        std::unordered_map<int, std::unordered_map<int, bool>> freq;
        for (int i = 0; i < arrays.size(); i++) {
            vector<size_t> array = arrays[i];
            for (int element : array) {
                if (freq.count(element)) {
                    freq[element][i] = true;
                } else {
                    freq[element] = std::unordered_map<int, bool>();
                    freq[element][i] = true;
                }
            }
        }
        vector<vector<size_t>> unique_arrays(arrays.size());
        for (auto pair : freq) {
            int element = pair.first;
            std::unordered_map<int, bool> freq_arr = pair.second;
            int count = 0;
            for (int i = 0; i < arrays.size(); i++) {
                if (freq_arr.count(i) && freq_arr[i]) {
                    count++;
                }
            }
            if (count == 1) {
                for (int i = 0; i < arrays.size(); i++) {
                    if (freq_arr.count(i) && freq_arr[i]) {
                        unique_arrays[i].push_back(element);
                    }
                }
            }
        }
        return unique_arrays;
    }
#ifdef OFEC_DEMO
    void DFDL::updateBuffer() {
        if (drawSolutionMode == "HisSolution") {
            m_solution.clear();
            m_solution.resize(1);
            for (const auto &sol : m_his_sols)
                m_solution[0].push_back(&(*sol));
        }
        else if (drawSolutionMode == "CurSolution") {
            m_solution.clear();
            m_solution.resize(m_sub_pop.size());
            if (m_swarm_Lips != nullptr) {
                m_solution.resize(m_sub_pop.size() + 1);
                for (size_t i = 0; i < m_swarm_Lips->size(); ++i) {
                    m_solution[m_solution.size() - 1].push_back(&m_swarm_Lips->at(i));
                }
            }
            for (size_t k = 0; k < m_sub_pop.size(); ++k) {
                for (size_t i = 0; i < m_sub_pop[k].size(); ++i) {
                    m_solution[k].push_back(&m_sub_pop[k][i]);
                }
            }
        }
        ofec_demo::g_buffer->appendAlgBuffer(this);
    }
#endif
}

#endif