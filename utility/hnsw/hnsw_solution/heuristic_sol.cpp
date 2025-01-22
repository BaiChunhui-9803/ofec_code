#include "heuristic_sol.h"
#include <vector>

#include "../min_heap.h"
namespace n2 {

	void NaiveNeighborSelectingPoliciesSol::Select(size_t m, bool select_nn, 
        ofec::Problem* pro,
        std::priority_queue<FurtherSolFirst>& result)
	{

		while (result.size() > m) {
			result.pop();
		}
	}

	void HeuristicNeighborSelectingPoliciesSol::Select(size_t m, bool select_nn,
        ofec::Problem* pro,
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
        MinHeap<double, HnswSolutionNode*> skipped;
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
                //if (j < picked.size() - 1) {
                //    //      _mm_prefetch(picked[j+1].GetNode()->GetData(), _MM_HINT_T0);
                //}
                ////     _mm_prefetch(cur_node->GetData(), _MM_HINT_T1);
                //if (dist_func_(cur_node, picked[j].GetNode(), dim) < cur_dist) {
                //    skip = true;
                //    break;
                //}

                //auto dis = cur_node->getOriginData()->variableDistance(
                //    *picked[j].GetNode()->getOriginData(), pro);
             //   cur_node->getOriginData()->variableDistance()
                if (cur_node->getOriginData()->variableDistance(
                    *picked[j].GetNode()->getOriginData(), pro) < cur_dist) {
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
                skipped.push(cur_dist, cur_node);
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
                result.emplace(skipped.top().data, skipped.top().key);
                skipped.pop();
            }
        }
    }
}
