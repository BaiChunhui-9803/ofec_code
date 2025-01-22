
#pragma once

#include <memory>
#include <queue>
#include <set>
#include <list>
//#include "spdlog/spdlog.h"


#include "heuristic_sol.h"
#include "../hnsw_origin/common.h"
#include "../../function/visit_lazy_list.h"
#include "../../../utility/general_multithread/general_multithread.h"




namespace n2 {

    
    class HnswSolModel {

    protected:

        class GlobalThreadInfo : public ofec::GeneralMultiThreadInfo {
        public:
            std::mutex m_info_mtx;
            std::set<HnswSolutionNode*> m_calculatingNodes;
        };


        class ThreadLocalInfo : public ofec::GeneralMultiThreadInfo {
        public:
            custom_fun::VisitedLazyList m_visisted;
            custom_fun::VisitedLazyList m_distanceVisited;
            std::vector<double> m_distanceMat;
        };


        class LocalThreadInfo : public ofec::GeneralMultiThreadInfo {
        public:
            int m_curNodeId = 0;
        };
    public:
        HnswSolModel();
        virtual ~HnswSolModel() {
            clear();
        };

        HnswSolModel(const HnswSolModel&) = delete;
        void operator=(const HnswSolModel&) = delete;

        void initialize(ofec::Environment *env, ofec::Random* rnd);

        int addDataSingleThread(ofec::SolutionBase* sol);
        void addDataMutliThread(
            std::vector<ofec::SolutionBase*> sols, 
            std::vector<int>& solId);
        
        void SetConfigs(const std::vector<std::pair<std::string, std::string>>& configs);

        void PrintDegreeDist() const;
        void PrintConfigs() const;
        const std::vector<EdgeInfo>& getNeighbor(int idx) const{
            return nodes_[idx]->getFriends(0);
        }

        int numberNodes()const {
            return nodes_.size();
        }



        void configsToParameters(std::list<ofec::ParameterVariant>& params);
        void configsfromParameters(std::list<ofec::ParameterVariant>& params);
        
        void nodesToParameters(std::list<ofec::ParameterVariant>& params);

        
        void getSolutions(std::vector<const ofec::SolutionBase*>& sols)const;
        void setSolutions(const std::vector<ofec::SolutionBase*>& sols);


        const std::vector<EdgeInfo>& getNeighbors(int nodeId, int level) const{
            return nodes_[nodeId]->getFriends(level);
        }
        
        void getTotalNeighbors(int nodeId, 
            std::vector<std::vector<int>>& neighbors) {
            nodes_[nodeId]->getNeighborIds(neighbors);
        }

    protected:

        void nodesToParameters(const HnswSolutionNode* node, std::list<ofec::ParameterVariant>& params);
        void nodesFromParameters(HnswSolutionNode* node, std::list<ofec::ParameterVariant>& params);


        void assignMemory(int size);
        void clear() {
            for (auto& it : nodes_) {
                delete it;
                it = nullptr;
            }
            nodes_.clear();
        }

        void SetConfigs(int m, int max_m0, int ef_construction, 
            float mult,
            NeighborSelectingPolicy neighbor_selecting, GraphPostProcessing graph_merging);

        void setInfo();

        int GetRandomNodeLevel(ofec::Random* rnd);
        void BuildGraph(bool reverse);

        virtual void InsertNodeSingleThread(HnswSolutionNode* qnode,
            custom_fun::VisitedLazyList* visited_list,
            custom_fun::VisitedLazyList* distanceVisited,
            std::vector<double>& distanceMemory) ;
        virtual void SearchAtLayerSingleThread(
            HnswSolutionNode* qnode, 
            const std::vector<HnswSolutionNode*>& enterpoints, 
            int level, custom_fun::VisitedLazyList* visited_list, 
            custom_fun::VisitedLazyList* distanceVisited,
            std::vector<double>& distanceMemory,
            std::priority_queue<FurtherSolFirst>& result) const;
        
        virtual void SearchAtLayerMultiThread(
            HnswSolutionNode* qnode, 
            const std::vector<HnswSolutionNode*>& enterpoints, 
            int level, custom_fun::VisitedLazyList* visited_list,
            custom_fun::VisitedLazyList* distanceVisited,
            std::vector<double>& distanceMemory,
            std::priority_queue<FurtherSolFirst>& result) const;


        virtual void Link(HnswSolutionNode* source, HnswSolutionNode* target, double dis, int level);
        virtual void LinkSingleThread(HnswSolutionNode* source, HnswSolutionNode* target, double dis, int level) ;
        virtual void MergeEdgesOfTwoGraphsSingleThread(HnswSolutionNode* cur,int level=0) ;

        void addDataTask(
            std::unique_ptr<ofec::GeneralMultiThreadInfo>& curInfo,
            std::unique_ptr<ofec::GeneralMultiThreadInfo>& threadInfo,
            std::unique_ptr<ofec::GeneralMultiThreadInfo>& globalInfo) ;
        


        double getDistance(
            HnswSolutionNode* curNode, HnswSolutionNode* otherNode, 
            ofec::Environment *env, 
            custom_fun::VisitedLazyList* visited,
            std::vector<double>& distance) const {
            if (visited->NotVisited(otherNode->GetId())) {
                visited->MarkAsVisited(otherNode->GetId());
                distance[otherNode->GetId()] = curNode->getOriginData().variableDistance(otherNode->getOriginData(), env);
            }
            return distance[otherNode->GetId()];
        }

    protected:
        //   std::shared_ptr<spdlog::logger> logger_;

        const std::string n2_signature = "TOROS_N2@N9R4";
        size_t m_ = 12;
        size_t max_m_ = 12;
        size_t max_m0_ = 24;
        size_t ef_construction_ = 150;
        ofec::Real level_mult_ = 1 / log(1.0 * m_);
       // int n_threads_ = 1;
        NeighborSelectingPolicy neighbor_selecting_ = NeighborSelectingPolicy::HEURISTIC;
        NeighborSelectingPolicy post_neighbor_selecting_ = NeighborSelectingPolicy::HEURISTIC_SAVE_REMAINS;
        GraphPostProcessing post_graph_process_ = GraphPostProcessing::SKIP;

        int max_level_ = 0;
        HnswSolutionNode* enterpoint_ = nullptr;
        mutable std::shared_mutex max_level_guard_;

        std::vector<HnswSolutionNode*> nodes_;


        bool is_naive_ = false;
        std::unique_ptr<BaseNeighborSelectingPoliciesSol> selecting_policy_;
        std::unique_ptr<BaseNeighborSelectingPoliciesSol> post_selecting_policy_;


        std::shared_ptr<ofec::Environment> m_env;
        std::shared_ptr<ofec::Random> m_rnd;
        custom_fun::VisitedLazyList m_visitedList;
        custom_fun::VisitedLazyList m_distanceCalculated;
        std::vector<double> m_distanceMemory;
    };

}