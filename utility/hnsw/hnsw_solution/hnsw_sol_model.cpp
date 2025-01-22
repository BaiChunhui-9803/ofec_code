#include "hnsw_sol_model.h"
#include "../../myexcept.h"
#include <unordered_set>
#include "../../../core/problem/problem.h"

#include "../../../utility/random/newran.h"
namespace n2 {
	HnswSolModel::HnswSolModel()
	{
	}
	void HnswSolModel::initialize(ofec::Problem* pro, ofec::Random* rnd)
	{
		clear();
		m_pro = pro->getSharedPtr();
		m_rnd = rnd->getSharedPtr();
		m_visitedList.initialize();
		if (neighbor_selecting_ == NeighborSelectingPolicy::HEURISTIC) {
			selecting_policy_ = std::make_unique<HeuristicNeighborSelectingPoliciesSol>(false);
			is_naive_ = false;
		}
		else if (neighbor_selecting_ == NeighborSelectingPolicy::HEURISTIC_SAVE_REMAINS) {
			selecting_policy_ = std::make_unique<HeuristicNeighborSelectingPoliciesSol>(true);
			is_naive_ = false;
		}
		else if (neighbor_selecting_ == NeighborSelectingPolicy::NAIVE) {
			selecting_policy_ = std::make_unique<NaiveNeighborSelectingPoliciesSol>();
			is_naive_ = true;
		}
		if (post_neighbor_selecting_ == NeighborSelectingPolicy::HEURISTIC) {
			post_selecting_policy_ = std::make_unique<HeuristicNeighborSelectingPoliciesSol>(false);
		}
		else if (post_neighbor_selecting_ == NeighborSelectingPolicy::HEURISTIC_SAVE_REMAINS) {
			post_selecting_policy_ = std::make_unique<HeuristicNeighborSelectingPoliciesSol>(true);
		}
		else if (post_neighbor_selecting_ == NeighborSelectingPolicy::NAIVE) {
			post_selecting_policy_ = std::make_unique<NaiveNeighborSelectingPoliciesSol>();
		}
	}
	int HnswSolModel::AddData(ofec::SolutionBase* sol)
	{

		int id = nodes_.size();
		int level = GetRandomNodeLevel(m_rnd.get());
		m_visitedList.insert(1);
		if (id == 0) {
			HnswSolutionNode* first = new HnswSolutionNode(0, sol, level, max_m_, max_m0_);
			//nodes_[0] = first;
			nodes_.push_back(first);
			max_level_ = level;
			enterpoint_ = first;
		}
		else {

			HnswSolutionNode* qnode = new HnswSolutionNode(id, sol, level, max_m_, max_m0_);
			nodes_.push_back(qnode);
			InsertNode(qnode,&m_visitedList );
		}
		return id;
	}
	void HnswSolModel::SetConfigs(const std::vector<std::pair<std::string, std::string>>& configs)
	{
		int m = -1, max_m0 = -1, ef_construction = -1, n_threads = -1;
		float mult = -1;
		NeighborSelectingPolicy neighbor_selecting = NeighborSelectingPolicy::HEURISTIC;
		GraphPostProcessing graph_merging = GraphPostProcessing::SKIP;

		for (const auto& c : configs) {
			if (c.first == "M") {
				m = stoi(c.second);
			}
			else if (c.first == "MaxM0") {
				max_m0 = stoi(c.second);
			}
			else if (c.first == "efConstruction") {
				ef_construction = stoi(c.second);
			}
			else if (c.first == "NumThread") {
				n_threads = stoi(c.second);
			}
			else if (c.first == "Mult") {
				mult = stof(c.second);
			}
			else if (c.first == "NeighborSelecting") {
				if (c.second == "heuristic") {
					neighbor_selecting = NeighborSelectingPolicy::HEURISTIC;
				}
				else if (c.second == "heuristic_save_remains") {
					neighbor_selecting = NeighborSelectingPolicy::HEURISTIC_SAVE_REMAINS;
				}
				else if (c.second == "naive") {
					neighbor_selecting = NeighborSelectingPolicy::NAIVE;
				}
				else {
					throw ofec::MyExcept("[Error] at HnswSolModel:: Invalid configuration value for NeighborSelecting: " + c.second);
				}
			}
			else if (c.first == "GraphMerging") {
				if (c.second == "skip") {
					graph_merging = GraphPostProcessing::SKIP;
				}
				else if (c.second == "merge_level0") {
					graph_merging = GraphPostProcessing::MERGE_LEVEL0;
				}
				else {
					throw ofec::MyExcept("[Error]@ HnswSolModel:: Invalid configuration value for GraphMerging: " + c.second);
				}
			}
			else if (c.first == "EnsureK") {
			}
			else {
				throw ofec::MyExcept("[Error]@HnswSolModel:: Invalid configuration key: " + c.first);
			}
		}

		SetConfigs(m, max_m0, ef_construction, n_threads, mult, neighbor_selecting, graph_merging);

	}
	void HnswSolModel::PrintDegreeDist() const
	{
	}
	void HnswSolModel::PrintConfigs() const
	{
	}
	void HnswSolModel::SetConfigs(int m, int max_m0, int ef_construction, int n_threads, float mult, NeighborSelectingPolicy neighbor_selecting, GraphPostProcessing graph_merging)
	{
		if (m > 0) max_m_ = m_ = m;
		if (max_m0 > 0) max_m0_ = max_m0;
		if (ef_construction > 0) ef_construction_ = ef_construction;
		if (n_threads > 0) n_threads_ = n_threads;
		level_mult_ = mult > 0 ? mult : 1 / log(1.0 * m_);
		neighbor_selecting_ = neighbor_selecting;
		post_graph_process_ = graph_merging;
	}
	int HnswSolModel::GetRandomNodeLevel(ofec::Random* rnd)
	{
		double r = rnd->uniform.next();

		if (r < std::numeric_limits<double>::epsilon())
			r = 1.0;
		return (int)(-log(r) * level_mult_);
	}
	void HnswSolModel::BuildGraph(bool reverse)
	{
	}
	void HnswSolModel::InsertNode(HnswSolutionNode* qnode, custom_fun::VisitedLazyList* visited_list)
	{
		using namespace std;
		int cur_level = qnode->GetLevel();

		int max_level_copy = max_level_;
		vector<HnswSolutionNode*> enterpoints;

		if (cur_level < max_level_copy) {
			HnswSolutionNode* cur_node = enterpoint_;
		//	double d = dist_func_(qnode, cur_node, data_dim_);
			double d = qnode->getOriginData()->variableDistance(*cur_node->getOriginData(), m_pro.get());
			double cur_dist = d;
			for (auto i = max_level_copy; i > cur_level; --i) {
				bool changed = true;
				while (changed) {
					changed = false;
				//	unique_lock<mutex> local_lock(cur_node->GetAccessGuard());
					const auto& neighbors = cur_node->GetFriends(i);
					for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
						//d = dist_func_(qnode, *iter, data_dim_);

						d = qnode->getOriginData()->variableDistance(*iter->GetNode()->getOriginData(), m_pro.get());

						if (d < cur_dist) {
							cur_dist = d;
							cur_node = iter->GetNode();
							changed = true;
						}
					}
				}
			}
			enterpoints.push_back(cur_node);
		}
		else {
			enterpoints.push_back(enterpoint_);
		}

		//  _mm_prefetch(&selecting_policy_, _MM_HINT_T0);
		for (auto i = min(max_level_copy, cur_level); i >= 0; --i) {
			priority_queue<FurtherSolFirst> result;
			SearchAtLayer(qnode, enterpoints, i, visited_list, result);

			enterpoints.clear();
			priority_queue<FurtherSolFirst> next_enterpoints = result;
			while (next_enterpoints.size() > 0) {
				auto* top_node = next_enterpoints.top().GetNode();
				next_enterpoints.pop();
				enterpoints.push_back(top_node);
			}

			selecting_policy_->Select(m_,  i == 0, m_pro.get(), result);
			while (result.size() > 0) {
				auto* top_node = result.top().GetNode();
		
				Link(top_node, qnode,result.top().GetDistance(), i);
				Link(qnode, top_node, result.top().GetDistance(), i);
				result.pop();
			}
		}

		if (cur_level > enterpoint_->GetLevel()) {
			enterpoint_ = qnode;
			max_level_ = cur_level;
		}
	}
	void HnswSolModel::SearchAtLayer(HnswSolutionNode* qnode, const std::vector<HnswSolutionNode*>& enterpoints, int level, custom_fun::VisitedLazyList* visited_list, std::priority_queue<FurtherSolFirst>& result)
	{
		using namespace std;
		priority_queue<CloserSolFirst> candidates;
		visited_list->Reset();

		for (const auto& enterpoint : enterpoints) {
			visited_list->MarkAsVisited(enterpoint->GetId());
			//float d = dist_func_(qnode, enterpoint, data_dim_);


			double d = qnode->getOriginData()->variableDistance(*enterpoint->getOriginData(), m_pro.get());

			result.emplace(enterpoint, d);
			candidates.emplace(enterpoint, d);
		}

		while (!candidates.empty()) {
			const CloserSolFirst& candidate = candidates.top();
			float lower_bound = result.top().GetDistance();
			if (candidate.GetDistance() > lower_bound)
				break;

			HnswSolutionNode* candidate_node = candidate.GetNode();
		//	unique_lock<mutex> lock(candidate_node->GetAccessGuard());
			const auto& neighbors = candidate_node->GetFriends(level);
			candidates.pop();

			for (const auto& neighbor : neighbors) {
				int id = neighbor.GetNode()->GetId();
				if (visited_list->NotVisited(id)) {
					//      _mm_prefetch(neighbor->GetData(), _MM_HINT_T0);
					visited_list->MarkAsVisited(id);
					double d = qnode->getOriginData()->variableDistance(*neighbor.GetNode()->getOriginData(), m_pro.get());
					//float d = dist_func_(qnode, neighbor, data_dim_);
					if (result.size() < ef_construction_ || result.top().GetDistance() > d) {
						result.emplace(neighbor.GetNode(), d);
						candidates.emplace(neighbor.GetNode(), d);
						if (result.size() > ef_construction_)
							result.pop();
					}
				}
			}
		}
	}
	void HnswSolModel::Link(HnswSolutionNode* source, HnswSolutionNode* target, double dis, int level)
	{
		using namespace std;
		auto& neighbors = source->GetFriends(level);
		neighbors.emplace_back(target,dis);
		bool shrink = (level > 0 && neighbors.size() > source->GetMaxM()) ||
			(level <= 0 && neighbors.size() > source->GetMaxM0());
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
			neighbors.erase(neighbors.begin() + maxi);
		}
		else {
			priority_queue<FurtherSolFirst> tempres;
			//for (const auto& neighbor : neighbors) {
			//	//    _mm_prefetch(neighbor->GetData(), _MM_HINT_T0);
			//}

			for (const auto& neighbor : neighbors) {
				tempres.emplace(neighbor);
				//tempres.emplace(neighbor, dist_func_(source, neighbor, data_dim_));
			}
			selecting_policy_->Select(tempres.size() - 1, level == 0,m_pro.get(), tempres);
			neighbors.clear();
			while (tempres.size()) {
				neighbors.emplace_back(tempres.top().GetNode(),tempres.top().GetDistance());
				tempres.pop();
			}
		}
	}
	void HnswSolModel::MergeEdgesOfTwoGraphs(HnswSolutionNode* cur)
	{
		//if(cur->GetId()
		int id = cur->GetId();
		const auto& neighbors1 = nodes_[id]->GetFriends(0);
		const auto& neighbors2 = cur->GetFriends(0);
		std::unordered_set<EdgeInfo> merged_neighbor_id_set = std::unordered_set<EdgeInfo>();
		for (const auto& cur : neighbors1) {
			merged_neighbor_id_set.insert(cur);
		}
		for (const auto& cur : neighbors2) {
			merged_neighbor_id_set.insert(cur);
		}
		std::priority_queue<FurtherSolFirst> temp_res;
		for (auto& cur : merged_neighbor_id_set) {

			temp_res.emplace(cur);
			//temp_res.emplace(nodes_[cur], dist_func_(nodes_[cur],
				//nodes_[i], data_dim_));
			//temp_res.emplace(nodes_[cur], dist_func_(nodes_[cur].GetRawData(), 
			 //                                        data_list_[i].GetRawData(), data_dim_));
		}

		// Post Heuristic
		post_selecting_policy_->Select(max_m0_,  true, m_pro.get(),temp_res);
		std::vector<EdgeInfo> merged_neighbors;
		while (!temp_res.empty()) {
			merged_neighbors.emplace_back(temp_res.top());
			temp_res.pop();
		}
		nodes_[id]->SetFriends(0, merged_neighbors);;
	}
}