#include "hnsw_sol_node.h"
#include <queue>
#include "sort_sol.h"
#include "heuristic_sol.h"
#include <mutex>
#include <memory>
#include <shared_mutex>

void n2::HnswSolutionNode::link(HnswSolutionNode* target, double dis, 
	int level, 
	bool is_naive_,
	BaseNeighborSelectingPoliciesSol* select_policy,
	ofec::Environment *env) {
	using namespace std;

	//std::vector<EdgeInfo> neighbors;
	{
		std::unique_lock lock(mutex_);
		auto& neighbors = friends_at_layer_[level];
		neighbors.emplace_back(target, dis);
		selectNeighbors(neighbors, level, max_m_, max_m0_, is_naive_,
			select_policy, env);

	}
	//

	//{
	//	SetFriends(level, neighbors);
	//}
}

void n2::HnswSolutionNode::linkSingleThread(HnswSolutionNode* target, double dis, 
	int level,
	bool is_naive_,
	BaseNeighborSelectingPoliciesSol* select_policy,
	ofec::Environment *env)
{
	using namespace std;
	auto& neighbors = friends_at_layer_[level];
	neighbors.emplace_back(target, dis);

	selectNeighbors(neighbors, level, max_m_, max_m0_, is_naive_,
		select_policy, env);
}


