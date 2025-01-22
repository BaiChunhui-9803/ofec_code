/******************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*******************************************************************************
* Author: Changhe Li
* Email: changhe.lw@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*
*  see https://github.com/Changhe160/OFEC for more information
*
*-------------------------------------------------------------------------------
* heuristic single-linkage hierarchical clustering
*********************************************************************************/
#ifndef OFEC_HSLH_CLUSTERING_H
#define OFEC_HSLH_CLUSTERING_H

#include "group.h"
#include <utility>
#include <vector>
#include <algorithm>
#include "../../core/algorithm/population.h"		//

namespace ofec {
	template<typename TInd>
	class HSLH {
		std::vector<std::vector<size_t>> m_clusters;
		std::vector<Group<TInd>> m_group;
		std::vector<std::vector<Real> > m_dis;	//distance between all objects

		Real m_inter_dis;			// inter distance of all group
		Real m_intra_dis;			// intra distance of all group
		int m_space;				//objective space or solution space

		typedef typename std::vector<Group<TInd>>::iterator mem_iter;

	public:
		HSLH(int space = 0) : m_space(space) {}

		HSLH(Population<TInd> &p, Environment *env, int space = 0) :
			m_space(space), 
			m_group(p.size()), 
			m_dis(p.size(), std::vector<Real>(p.size())) 
		{	//construct function	
			setData(p, env);
		}

		HSLH(std::vector<std::unique_ptr<TInd>>& p, Environment *env, int space = 0) :
			m_space(space), 
			m_group(p.size()), 
			m_dis(p.size(), std::vector<Real>(p.size())) 
		{	//construct function	
			setData(p, env);
		}

		void setData(Population<TInd>& pop, Environment *env) {
			if (m_group.size() != pop.size())
				m_group.resize(pop.size());
			if (m_dis.size() != pop.size()) {
				m_dis.assign(pop.size(), std::vector<Real>(pop.size()));
			}
			for (int k = 0; k < pop.size(); ++k) {
				m_group[k].initialize(std::make_pair(&pop[k], k));
			}
			for (int i = 0; i < pop.size(); i++) {					//initialize m_dis
				m_dis[i][i] = 0;
				for (int j = 0; j < i; j++) {
					if (m_space == 0)
						m_dis[i][j] = m_dis[j][i] = m_group[i].begin()->second->variableDistance(*m_group[j].begin()->second, env);
					else
						m_dis[i][j] = m_dis[j][i] = m_group[i].begin()->second->objectiveDistance(*m_group[j].begin()->second);
				}
			}
			updateDistance();
		}

		void setData(std::vector<std::unique_ptr<TInd>>& inds, Environment *env) {
			if (m_group.size() != inds.size())
				m_group.resize(inds.size());
			if (m_dis.size() != inds.size()) {
				m_dis.assign(inds.size(), std::vector<Real>(inds.size()));
			}
			for (int k = 0; k < inds.size(); ++k) {
				m_group[k].initialize(std::make_pair(inds[k].get(), k));
			}
			for (int i = 0; i < inds.size(); i++) {					//initialize m_dis
				m_dis[i][i] = 0;
				for (int j = 0; j < i; j++) {
					if (m_space == 0)
						m_dis[i][j] = m_dis[j][i] = m_group[i].begin()->second->variableDistance(*m_group[j].begin()->second, env);
					else
						m_dis[i][j] = m_dis[j][i] = m_group[i].begin()->second->objectiveDistance(*m_group[j].begin()->second);
				}
			}
			updateDistance();
		}

		void setData(std::vector<TInd*>& inds, Environment *env) {
			if (m_group.size() != inds.size())
				m_group.resize(inds.size());
			if (m_dis.size() != inds.size()) {
				m_dis.assign(inds.size(), std::vector<Real>(inds.size()));
			}
			for (int k = 0; k < inds.size(); ++k) {
				m_group[k].initialize(std::make_pair(inds[k], k));
			}
			for (int i = 0; i < inds.size(); i++) {					//initialize m_dis
				m_dis[i][i] = 0;
				for (int j = 0; j < i; j++) {
					if (m_space == 0)
						m_dis[i][j] = m_dis[j][i] = m_group[i].begin()->second->variableDistance(*m_group[j].begin()->second, env);
					else
						m_dis[i][j] = m_dis[j][i] = m_group[i].begin()->second->objectiveDistance(*m_group[j].begin()->second);
				}
			}
			updateDistance();
		}

		void setData(std::vector<const TInd*> &inds, Environment *env) {
			if (m_group.size() != inds.size())
				m_group.resize(inds.size());
			if (m_dis.size() != inds.size()) {
				m_dis.assign(inds.size(), std::vector<Real>(inds.size()));
			}
			for (int k = 0; k < inds.size(); ++k) {
				m_group[k].initialize(std::make_pair(inds[k], k));
			}
			for (int i = 0; i < inds.size(); i++) {					//initialize m_dis
				m_dis[i][i] = 0;
				for (int j = 0; j < i; j++) {
					if (m_space == 0)
						m_dis[i][j] = m_dis[j][i] = m_group[i].begin()->second->variableDistance(*m_group[j].begin()->second, env);
					else
						m_dis[i][j] = m_dis[j][i] = m_group[i].begin()->second->objectiveDistance(*m_group[j].begin()->second);
				}
			}
			updateDistance();
		}

		void updateDistance() {
			m_inter_dis = m_intra_dis = 0;
			//calculate inter-distance
			for (unsigned i = 0; i < m_group.size(); i++) {
				for (unsigned j = 0; j < i; j++) {
					m_inter_dis += m_dis[m_group[i].best()][m_group[j].best()];
				}
			}
		
			m_inter_dis /= m_group.size()*(m_group.size() - 1) / 2;		//average inter distance

			// calulate intra-distance
			for (unsigned i = 0; i < m_group.size(); i++) {
				m_intra_dis += m_group[i].intraDis();
			}

			m_intra_dis /= m_group.size();				//average intra distance
		}
		
		void clustering2(size_t num_clusters, Environment *env) {
			while (m_group.size() > num_clusters) {
				auto p = nearest_group(-1);		//find nearest two group

				if ((p.first == m_group.end()))
					break;

				p.first->merge(*p.second, m_dis, env);
				m_group.erase(p.second);

				updateDistance();				//update distance

				if (m_group.size() <= num_clusters)
					break;
			}
			m_clusters.clear();
			m_clusters.resize(m_group.size());
			for (size_t i = 0; i < m_group.size(); i++) {
				for (auto iter = m_group[i].begin(); iter != m_group[i].end(); ++iter)
					m_clusters[i].push_back(iter->first);
			}
		}

		void clustering(int minsize, Environment *env) {
			while (1) {
				mem_iter i = m_group.begin();
				while (i != m_group.end() && i->size() >= minsize) 
					i++;
				if (i == m_group.end())
					break;

				auto p = nearest_group(-1);		//find nearest two group

				if ((p.first == m_group.end())) 
					break;

				p.first->merge(*p.second, m_dis, env);
				m_group.erase(p.second);

				updateDistance();				//update distance

				if (m_inter_dis <= m_intra_dis)
					break;

			}
			m_clusters.clear();
			m_clusters.resize(m_group.size());
			for (size_t i = 0; i < m_group.size(); i++)	{
				for (auto iter = m_group[i].begin(); iter != m_group[i].end(); ++iter)
					m_clusters[i].push_back(iter->first);
			}
		}

		void clustering(int maxsize, int minsize, Environment *env) {
			while (1) {
				mem_iter i = m_group.begin();
				while (i != m_group.end() && i->size() >= minsize) 
					i++;
				if (i == m_group.end()) 
					break;

				auto &&p = nearest_group(maxsize);		//find nearest two group

				if ((p.first == m_group.end()) ) 
					break;
				p.first->merge(*p.second, m_dis, env);
				m_group.erase(p.second);

				updateDistance();				//update distance

				if (m_inter_dis <= m_intra_dis)
					break;
			}
			m_clusters.clear();
			m_clusters.resize(m_group.size());
			for (size_t i = 0; i < m_group.size(); i++) {
				for (auto iter = m_group[i].begin(); iter != m_group[i].end(); ++iter)
					m_clusters[i].push_back(iter->first);
			}
		}

		const std::vector<std::vector<size_t>>& clusters() const {
			return m_clusters; 
		}

		Group<TInd>& operator[](const int i) {
			return m_group[i];
		}

		int size() {					//size of group
			return m_group.size();
		}

		int size_above(int n) {
			int count = 0;
			for (auto i = m_group.begin(); i != m_group.end(); ++i) {
				if (i->size() > n) ++count;
			}
			return count;
		}

		std::pair<mem_iter, mem_iter> nearest_group(int maxsize) {
			
			Real Min_dis = std::numeric_limits<Real>::max(), dist;
			auto g1 = m_group.end(), g2 = m_group.end();

			for (auto i = m_group.begin(); i != m_group.end(); ++i) {
				// can't merge two mp_groups whose total m_number are greater than m_maxsize
				for (auto j = i+1; j != m_group.end(); ++j) {		
					if (maxsize > 0 && (i->size() + j->size()) > maxsize) continue;

					dist = m_dis[i->best()][j->best()];
					if (Min_dis > dist) {
						Min_dis = dist;
						g1 = i;
						g2 = j;
					}
				}
			}

			return std::make_pair(g1, g2);
		}

		std::pair<mem_iter, mem_iter> nearest_group_rough(int maxsize) {

			Real Min_dis = std::numeric_limits<Real>::max(), dist;
			auto g1 = m_group.end(), g2 = m_group.end();

			for (auto i = m_group.begin(); i != m_group.end(); ++i) {
				// can't merge two mp_groups whose total m_number are greater than m_maxsize
				for (auto j = i + 1; j != m_group.end(); ++j) {
					if (maxsize > 0 && (i->size() + j->size()) > maxsize) continue;

					dist = m_dis[i->begin()->first][j->begin()->first];
					for (auto pi = i->begin(); pi != i->end(); ++pi) {
						for (auto pj = j->begin(); pj != j->end(); ++pj) {
							if(dist < m_dis[pi->first][pj->first]){
								dist = m_dis[pi->first][pj->first];
							}
						}
					}

					if (Min_dis > dist) {
						Min_dis = dist;
						g1 = i;
						g2 = j;
					}
				}
			}

			return std::make_pair(g1, g2);
		}

		void clusteringRough(int maxsize, int minsize, Environment *env) {
			while (1) {
				mem_iter i = m_group.begin();
				while (i != m_group.end() && i->size() > 1)
					i++;
				if (i == m_group.end())
					break;

				auto&& p = nearest_group_rough(maxsize);		//find nearest two group

				if ((p.first == m_group.end()))
					break;
				p.first->merge(*p.second, m_dis, env);
				m_group.erase(p.second);
			}
			m_clusters.clear();
			m_clusters.resize(m_group.size());
			for (size_t i = 0; i < m_group.size(); i++) {
				for (auto iter = m_group[i].begin(); iter != m_group[i].end(); ++iter)
					m_clusters[i].push_back(iter->first);
			}
		}

        void clusteringInFixedNum(int cluster_num, Environment *env) {
            while (1) {
                mem_iter i = m_group.begin();

                auto p = nearest_group(-1);		//find nearest two group

                if ((p.first == m_group.end()))
                    break;

                p.first->merge(*p.second, m_dis, env);
                m_group.erase(p.second);
/*
                // 每个粒子的平均类内距离 a
                std::vector<std::vector<Real>> intra_dis;
                std::vector<std::vector<Real>> intramax_dis;
                for (int n = 0; n < m_group.size(); ++n) {
                    intra_dis.resize(m_group.size());
                    intramax_dis.resize(m_group.size());
                    for (auto it = m_group[n].begin(); it != m_group[n].end(); it++) {
                        std::vector<Real> dis_record;
                        Real dis = 0.;
                        for (auto it2 = m_group[n].begin(); it2 != m_group[n].end(); it2++) {
                            if(it != it2) {
                                dis += m_dis[it->first][it2->first];
                                dis_record.push_back(m_dis[it->first][it2->first]);
                            }
                        }
                        if(m_group[n].size() == 1) {
                            intra_dis[n].push_back(dis / m_group[n].size());
                        }
                        else {
                            intra_dis[n].push_back(dis / (m_group[n].size() -1));
                        }
                        if(dis_record.size() == 0) {
                            intramax_dis[n].push_back(0);
                        }
                        else {
                            intramax_dis[n].push_back(*std::max_element(dis_record.begin(),dis_record.end()));
                        }
                    }
                }
                // 每个粒子的平均类间 b
                std::vector<std::vector<Real>> inter_dis;
                std::vector<std::vector<Real>> intermax_dis;
                inter_dis.resize(m_group.size());
                intermax_dis.resize(m_group.size());
                for (int n = 0; n < m_group.size(); ++n) {
                    for (auto it = m_group[n].begin(); it != m_group[n].end(); it++) {
                        Real dis = 0.;
                        std::vector<Real> dis_record;
                        int count = 0;
                        for (int m = 0; m < m_group.size(); ++m) {
                            if (n != m) {
                                for (auto it2 = m_group[m].begin(); it2 != m_group[m].end(); it2++) {
                                    count++;
                                    dis += m_dis[it->first][it2->first];
                                    dis_record.push_back(m_dis[it->first][it2->first]);
                                }
                            }
                        }
                        if(m_group.size() == 1) {
                            intermax_dis[n].push_back(0);
                        }
                        else {
                            intermax_dis[n].push_back(*std::max_element(dis_record.begin(),dis_record.end()));
                        }
                        if(count == 0) {
                            inter_dis[n].push_back(0);
                        }
                        else {
                            inter_dis[n].push_back(dis / count);
                        }
                    }
                }

                // (b-a)/max{b-a}
                int count = 0;
                Real sc = 0.;
                for (int j = 0; j < inter_dis.size(); ++j) {
                    for (int k = 0; k < inter_dis[j].size(); ++k) {
                        count++;
                        sc += (inter_dis[j][k] - intra_dis[j][k]) / std::max(intermax_dis[j][k], intramax_dis[j][k]);
                    }
                }
                sc = sc/count;
                std::cout << "silhouette coefficient: " << sc << "\t" << m_group.size() << std::endl;
*/
                updateDistance();				//update distance

                if (m_group.size() == cluster_num)
                    break;
            }
            m_clusters.clear();
            m_clusters.resize(m_group.size());
            for (size_t i = 0; i < m_group.size(); i++)	{
                for (auto iter = m_group[i].begin(); iter != m_group[i].end(); ++iter)
                    m_clusters[i].push_back(iter->first);
            }
        }
	};
}

#endif // OFEC_HSLH_CLUSTERING_H
