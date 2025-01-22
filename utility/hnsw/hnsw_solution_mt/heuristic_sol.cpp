#include "heuristic_sol.h"
#include <vector>
#include <queue>


namespace n2 {

	void NaiveNeighborSelectingPoliciesSol::Select(size_t m, bool select_nn, 
        ofec::Environment *env,
        std::priority_queue<FurtherSolFirst>& result)
	{

		while (result.size() > m) {
			result.pop();
		}
	}

	void HeuristicNeighborSelectingPoliciesSol::Select(size_t m, bool select_nn,
        ofec::Environment *env,
        std::priority_queue<FurtherSolFirst>& result)
    {
        if (result.size() <= m) return;

        size_t nn_num = 0;  // # of nearest neighbors
        if (select_nn) {
            nn_num = (size_t)(m * 0.25);
            // m - nn_num neighbors will be chosen as usual with the heuristic algorithm
        }
        size_t nn_picked_num = 0;  // # of nearest neighbors also picked by the heuristic algorithm
        // nn_num - nn_picked_num = # of nearest neighbors selected unconditionally but not picked by the heuristic algorithm 

        std::vector<FurtherSolFirst> neighbors, picked;
       // MinHeap<double, HnswSolutionNode*> skipped;
        std::priority_queue<n2::CloserSolFirst> skipped;
         

        while (!result.empty()) {
            neighbors.push_back(result.top());
            result.pop();
        }

        for (auto it = neighbors.rbegin(); it != neighbors.rend(); it++) {
            double cur_dist = it->GetDistance();
            HnswSolutionNode* cur_node = it->GetNode();
            //    _mm_prefetch(cur_node->GetData(), _MM_HINT_T0);
            bool nn_selected = false;
            if (result.size() < nn_num) {
                result.emplace(*it);
                nn_selected = true;
            }

            bool skip = false;
            for (size_t j = 0; j < picked.size(); ++j) {
                if (cur_node->getOriginData().variableDistance(
                    picked[j].GetNode()->getOriginData(), env) < cur_dist) {
                    skip = true;
                    break;
                }
            }

            if (!skip) {
                picked.push_back(*it);
                if (nn_selected) {
                    // nearest neighbors included in result & picked by the heuristic algorithm
                    ++nn_picked_num;
                }
            }
            else if (!nn_selected && save_remains_) {
                skipped.emplace(cur_node, cur_dist);
            }

            if (picked.size() - nn_picked_num == m - nn_num)
                // check if # of neighbors exclusively picked by the heuristic algorithm equals m - nn_num
                break;
        }

        for (size_t i = nn_picked_num; i < picked.size(); ++i) {
            result.emplace(picked[i]);
        }

        if (save_remains_) {
            while (result.size() < m && skipped.size()) {
                result.emplace(skipped.top().GetNode(), skipped.top().GetDistance());
                skipped.pop();
            }
        }
    }
    void selectNeighbors(std::vector<EdgeInfo>& neighbors, 
        int level, int maxM, int maxM0, 
        bool is_naive_, 
        BaseNeighborSelectingPoliciesSol* select_policy, 
        ofec::Environment *env)
    {

        bool shrink = (level > 0 && neighbors.size() > maxM) ||
            (level <= 0 && neighbors.size() > maxM0);
        if (!shrink) return;
        if (is_naive_) {
            //auto max = dist_func_(source, neighbors[0], data_dim_);
            auto max = neighbors[0].GetDistance();
            int maxi = 0;
            for (size_t i = 1; i < neighbors.size(); ++i) {
                //float curd = dist_func_(source, neighbors[i], data_dim_);
                auto curd = neighbors[i].GetDistance();
                if (curd > max) {
                    max = curd;
                    maxi = i;
                }
            }
            std::swap(neighbors[maxi], neighbors.back());
            neighbors.pop_back();
            //neighbors.erase(neighbors.begin() + maxi);
        }
        else {
            std::priority_queue<FurtherSolFirst> tempres;

            for (const auto& neighbor : neighbors) {
                tempres.emplace(neighbor);
            }
            select_policy->Select(tempres.size() - 1, level == 0, env, tempres);
            neighbors.clear();
            while (tempres.size()) {
                neighbors.emplace_back(tempres.top().GetNode(), tempres.top().GetDistance());
                tempres.pop();
            }
        }
    }
}
