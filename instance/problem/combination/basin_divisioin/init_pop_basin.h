#ifndef OFEC_INIT_POP_BASIN_H
#define OFEC_INIT_POP_BASIN_H

#include "../../../../core/problem/problem.h"
#include "../../../../utility/hnsw/hnsw_nbn/hnsw_model.h"


namespace ofec {

#define CAST_NBN_Basin(pro) dynamic_cast<InitNBN_BasinBase*>(pro)

    class InitNBN_BasinBase : virtual public Problem {
        OFEC_ABSTRACT_INSTANCE(InitNBN_BasinBase)
    protected:

        nbn::HnswModel m_hnswModel;
        std::vector<ofec::Real> m_dis2parent;
        std::vector<int> m_belong;
        std::vector<std::vector<int>> m_basinIds;
        std::vector<int> m_belongBasinIds;
        //std::vector<std::unique_ptr<const SolutionBase>> m_sols_init_pop; // solutions for intializing population
        int m_currentBasinId = 0;

        std::vector<ofec::Real> m_basinBestFit;
    public:

        size_t numberSolutions()const {
            return m_belong.size();
        }


        void addInputParameters() {
            m_input_parameters.add("Search basin Id", new RangedInt(m_currentBasinId, 0, 100, 0));
        }
        

        void setBasins(const std::vector<std::vector<int>>& basins) {


     
            m_basinIds = basins;
            auto& sols = m_hnswModel.solutions();
            m_belongBasinIds.resize(m_dis2parent.size());
            m_basinBestFit.resize(m_basinIds.size());
            std::fill(m_basinBestFit.begin(), m_basinBestFit.end(), std::numeric_limits<Real>::lowest());
            for (int idBelong(0); idBelong < m_basinIds.size(); ++idBelong) {
                for (auto& it : m_basinIds[idBelong]) {
                    m_belongBasinIds[it] = idBelong;
                    auto& cursol = sols[it];
                    m_basinBestFit[idBelong] = std::min(m_basinBestFit[idBelong], cursol->fitness());
                }
                
            }
            //m_basins.resize(basins.size());
            //for (int idx(0); idx < m_basins.size(); ++idx) {
            //    m_basins[idx].resize(basins[idx].size());
            //    for (int idy(0); idy < basins[idx].size(); ++idy) {
            //        m_basins[idx][idy].reset(createSolution(*basins[idx][idy]));;
            //    }
            //}
         //   m_sols_init_pop.emplace_back(createSolution(sol));
        }

        void setBasinId(int basinId) {
            m_currentBasinId = basinId;
        }

        int currentBasinId() {
            return m_currentBasinId;
        };
        
        const nbn::HnswModel& hnswModel() const{
            return m_hnswModel;
        }
        nbn::HnswModel& hnswModel() {
            return m_hnswModel;
        }
        const std::vector<ofec::Real>& dis2parent()const {
            return m_dis2parent;
        }

        std::vector<ofec::Real>& dis2parent(){
            return m_dis2parent;
        }
        const std::vector<int>& belong()const {
            return m_belong;
        }

        std::vector<int>& belong(){
            return m_belong;
        }
        const std::vector<std::vector<int>>& basinIds()const {
            return m_basinIds;
        }
        const std::vector<int>& belongBasinIds()const {
            return m_belongBasinIds;
        }
        void setOneBasin() {
            m_basinIds.resize(1);
            m_basinIds.front().clear();
            for (int idx(0); idx < m_dis2parent.size(); ++idx) {
                m_basinIds.front().push_back(idx);
            }
        }
        
    };
}


#endif