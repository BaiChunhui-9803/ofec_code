
#pragma once

#include <memory>
#include <queue>

//#include "spdlog/spdlog.h"

#include "../common.h"
#include "heuristic_sol.h"
#include "../../function/visit_lazy_list.h"
#include "hnsw_sol_node.h"

namespace n2 {




    class HnswSolModel {

    protected:
        void clear() {
            for (auto& it : nodes_) {
                delete it;
                it = nullptr;
            }
            nodes_.clear();
        }
    public:
        HnswSolModel();
        virtual ~HnswSolModel() {
            clear();
        };

        HnswSolModel(const HnswSolModel&) = delete;
        void operator=(const HnswSolModel&) = delete;

        void initialize(ofec::Problem* pro, ofec::Random* rnd);

        int AddData(ofec::SolutionBase* sol);
        void SetConfigs(const std::vector<std::pair<std::string, std::string>>& configs);

        void PrintDegreeDist() const;
        void PrintConfigs() const;
        const std::vector<EdgeInfo>& getNeighbor(int idx) {
            return nodes_[idx]->GetFriends(0);
        }

        

    protected:
        void SetConfigs(int m, int max_m0, int ef_construction, int n_threads, float mult,
            NeighborSelectingPolicy neighbor_selecting, GraphPostProcessing graph_merging);


        int GetRandomNodeLevel(ofec::Random* rnd);
        void BuildGraph(bool reverse);

        virtual void InsertNode(HnswSolutionNode* qnode, custom_fun::VisitedLazyList* visited_list) ;
        virtual void SearchAtLayer(HnswSolutionNode* qnode, const std::vector<HnswSolutionNode*>& enterpoints, int level,
            custom_fun::VisitedLazyList* visited_list, std::priority_queue<FurtherSolFirst>& result) ;
        virtual void Link(HnswSolutionNode* source, HnswSolutionNode* target,double dis, int level) ;
        virtual void MergeEdgesOfTwoGraphs(HnswSolutionNode* cur) ;




    protected:
        //   std::shared_ptr<spdlog::logger> logger_;

        const std::string n2_signature = "TOROS_N2@N9R4";
        size_t m_ = 12;
        size_t max_m_ = 12;
        size_t max_m0_ = 24;
        size_t ef_construction_ = 150;
        float level_mult_ = 1 / log(1.0 * m_);
        int n_threads_ = 1;
        NeighborSelectingPolicy neighbor_selecting_ = NeighborSelectingPolicy::HEURISTIC;
        NeighborSelectingPolicy post_neighbor_selecting_ = NeighborSelectingPolicy::HEURISTIC_SAVE_REMAINS;
        GraphPostProcessing post_graph_process_ = GraphPostProcessing::SKIP;

        int max_level_ = 0;
        HnswSolutionNode* enterpoint_ = nullptr;
        mutable std::mutex max_level_guard_;

        std::vector<HnswSolutionNode*> nodes_;
        int num_nodes_ = 0;



        bool is_naive_ = false;
        std::unique_ptr<BaseNeighborSelectingPoliciesSol> selecting_policy_;
        std::unique_ptr<BaseNeighborSelectingPoliciesSol> post_selecting_policy_;


        std::shared_ptr<ofec::Problem> m_pro;
        std::shared_ptr<ofec::Random> m_rnd;
        custom_fun::VisitedLazyList m_visitedList;
    };

}