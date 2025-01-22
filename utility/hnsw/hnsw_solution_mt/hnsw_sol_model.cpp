#include "hnsw_sol_model.h"
#include <unordered_set>
#include "../../../core/exception.h"
#include "../../../core/problem/problem.h"
#include "../../../core/random/newran.h"
namespace n2 {
	HnswSolModel::HnswSolModel()
	{

	}
	void HnswSolModel::initialize(ofec::Environment *env, ofec::Random* rnd)
	{
		clear();

		m_env = env->getSharedPtr();
		m_rnd = rnd->getSharedPtr();

		setInfo();
		
	}

	void HnswSolModel::setInfo() {
		m_visitedList.initialize();
		m_distanceCalculated.initialize();
		
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

	void HnswSolModel::configsToParameters(std::list<ofec::ParameterVariant>& params) {
		params.emplace_back(m_);
		params.emplace_back(max_m_);
		params.emplace_back(max_m0_);
		params.emplace_back(ef_construction_);
		params.emplace_back(level_mult_);
		//params.emplace_back(n_threads_);

		params.emplace_back(static_cast<int>(neighbor_selecting_));
		params.emplace_back(static_cast<int>(post_neighbor_selecting_));
		params.emplace_back(static_cast<int>(post_graph_process_));

		if (enterpoint_ == nullptr) params.emplace_back(int(-1));
		else params.emplace_back(enterpoint_->GetId());
	
		params.emplace_back(nodes_.size());
	}
	 
	 
	void HnswSolModel::configsfromParameters(std::list<ofec::ParameterVariant>& params) {
		using namespace ofec;
	

		size_t node_size(0);
		int enterId(-1);
		int policyType(0);
		
		getFrom(params.back(), node_size);
		params.pop_back();
		getFrom(params.back(), enterId);
		params.pop_back();


		getFrom(params.back(), policyType);
		params.pop_back();
		post_graph_process_ = static_cast<GraphPostProcessing>(policyType);

		getFrom(params.back(), policyType);
		params.pop_back();
		post_neighbor_selecting_ = static_cast<NeighborSelectingPolicy>(policyType);

		getFrom(params.back(), policyType);
		params.pop_back();
		neighbor_selecting_ = static_cast<NeighborSelectingPolicy>(policyType);

		getFrom(params.back(), level_mult_);
		params.pop_back();

		getFrom(params.back(), ef_construction_);
		params.pop_back();

		getFrom(params.back(), max_m0_);
		params.pop_back();

		getFrom(params.back(), max_m_);
		params.pop_back();

		getFrom(params.back(), m_);
		params.pop_back();

		setInfo();
		
		assignMemory(node_size);
		
		if (enterId == -1) enterpoint_ = nullptr;
		else enterpoint_ = nodes_[enterId];
	}


	void HnswSolModel::nodesToParameters(std::list<ofec::ParameterVariant>& params) {
		

	
	}



	void HnswSolModel::nodesToParameters(const HnswSolutionNode* node, std::list<ofec::ParameterVariant>& params) {
		auto& friends_at_layer_ = node->getTotalFriends();
		params.emplace_back(node->GetLevel());
		//params.emplace_back(friends_at_layer_.size());

		std::vector<ofec::Real> value;
		std::vector<int> nodeId;
		for (size_t idx(0); idx < friends_at_layer_.size(); ++idx) {
			auto& curlayer = friends_at_layer_[idx];

			value.resize(curlayer.size());
			nodeId.resize(curlayer.size());
			for (int idy(0); idy < curlayer.size(); ++idy) {
				value[idy] = curlayer[idy].GetDistance();
				nodeId[idy] = curlayer[idy].GetNode()->GetId();
			}
			params.emplace_back(value);
			params.emplace_back(nodeId);
			
		}
	}

	void HnswSolModel::nodesFromParameters(HnswSolutionNode* node, std::list<ofec::ParameterVariant>& params) {
		using namespace ofec;
		auto& friends_at_layer_ = node->getTotalFriends();


		std::vector<ofec::Real> value;
		std::vector<int> nodeId;
		for (size_t idx(friends_at_layer_.size()-1); idx >=0; --idx) {

			getFrom(params.back(), nodeId);
			params.pop_back();
			getFrom(params.back(), value);
			params.pop_back();

			auto& curlayer = friends_at_layer_[idx];

			for (int idy(0); idy < curlayer.size(); ++idy) {
				curlayer[idy].set(nodes_[nodeId[idy]], value[idy]);

			}
		}

		getFrom(params.back(), node->level());
		params.pop_back();

	}

	void HnswSolModel::assignMemory(int size) {
		clear();
		m_visitedList.resize(size);;
		m_distanceCalculated.resize(size);
		m_distanceMemory.resize(size);

		nodes_.resize(size);
		int level = 0;
		for (int idx(0); idx < size; ++idx) {
			
			int id = nodes_.size();
			HnswSolutionNode* curnode = new HnswSolutionNode(0, nullptr, level, max_m_, max_m0_);
			nodes_[idx] = curnode;
		}
	}

	void HnswSolModel::getSolutions(std::vector<const ofec::SolutionBase*>& sols)const {
		sols.resize(nodes_.size());
		for (int idx(0); idx < sols.size(); ++idx) {
			sols[idx] = &nodes_[idx]->getOriginData();
		}
	}
	void HnswSolModel::setSolutions(const std::vector<ofec::SolutionBase*>& sols) {

		for (int idx(0); idx < sols.size(); ++idx) {
			nodes_[idx]->setSolution(sols[idx]);
		}
	}
	

	int HnswSolModel::addDataSingleThread(ofec::SolutionBase* sol)
	{

		int id = nodes_.size();
		int level = GetRandomNodeLevel(m_rnd.get());
		m_visitedList.insert(1);
		m_distanceCalculated.insert(1);
		m_distanceMemory.push_back(0);
		m_distanceCalculated.Reset();
		
	
		if (id == 0) {
			HnswSolutionNode* first = new HnswSolutionNode(0, sol, level, max_m_, max_m0_);
			nodes_.push_back(first);
			{
			//	std::unique_lock lock(max_level_guard_);
				max_level_ = level;
				enterpoint_ = first;
			}
		}
		else {

			HnswSolutionNode* qnode = new HnswSolutionNode(id, sol, level, max_m_, max_m0_);
			nodes_.push_back(qnode);
			InsertNodeSingleThread(qnode,&m_visitedList, &m_distanceCalculated, m_distanceMemory);
		}
		return id;
	}
	void HnswSolModel::addDataMutliThread(std::vector<ofec::SolutionBase*> sols, std::vector<int>& solId) {
		//std::vector<ofec::Node*> sols;
		int fromId = nodes_.size();
		nodes_.resize(nodes_.size() + sols.size());
		int toId = nodes_.size();
		solId.clear();
		for (int id(fromId);  id< toId; ++id) {
			int level = GetRandomNodeLevel(m_rnd.get());
			HnswSolutionNode* qnode = new HnswSolutionNode(id, sols[id-fromId], level, max_m_, max_m0_);
			nodes_[id]=qnode;
			solId.push_back(id);
		}


		if (fromId == 0) {
			max_level_ = nodes_.front()->GetLevel();
			enterpoint_ = nodes_.front();
			++fromId;
		}

		ofec::GeneralMultiThread curThread;
		curThread.ms_global_info.reset(new GlobalThreadInfo());
		auto& globalInfo(dynamic_cast<GlobalThreadInfo&>(*curThread.ms_global_info));


		ThreadLocalInfo threadinfo;
		threadinfo.m_visisted.resize(nodes_.size());
		threadinfo.m_distanceVisited.resize(nodes_.size());
		threadinfo.m_distanceMat.resize(nodes_.size());
		for (auto& it : curThread.ms_theadInfo) {

			it.reset(new ThreadLocalInfo(threadinfo));
		}
		
		LocalThreadInfo localInfo;

		for (int id(fromId); id < toId; ++id) {
			localInfo.m_error_name = "calculating sol model idx " + std::to_string(id);
			localInfo.m_curNodeId = id;
			curThread.ms_info_buffer.emplace_back(new LocalThreadInfo(localInfo));
		}


		curThread.ms_fun = std::bind(&HnswSolModel::addDataTask,this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
		// out(totalInfo.m_save_dir + "error.txt");
		curThread.ms_cout_flag = false;
		curThread.runNoException();
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
			//else if (c.first == "NumThread") {
			//	n_threads = stoi(c.second);
			//}
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
					throw ofec::Exception("[Error] at HnswSolModel:: Invalid configuration value for NeighborSelecting: " + c.second);
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
					throw ofec::Exception("[Error]@ HnswSolModel:: Invalid configuration value for GraphMerging: " + c.second);
				}
			}
			else if (c.first == "EnsureK") {
			}
			else {
				throw ofec::Exception("[Error]@HnswSolModel:: Invalid configuration key: " + c.first);
			}
		}

		SetConfigs(m, max_m0, ef_construction, mult, neighbor_selecting, graph_merging);

	}
	void HnswSolModel::PrintDegreeDist() const
	{
	}
	void HnswSolModel::PrintConfigs() const
	{
	}
	void HnswSolModel::SetConfigs(int m, int max_m0, int ef_construction, float mult, NeighborSelectingPolicy neighbor_selecting, GraphPostProcessing graph_merging)
	{
		if (m > 0) max_m_ = m_ = m;
		if (max_m0 > 0) max_m0_ = max_m0;
		if (ef_construction > 0) ef_construction_ = ef_construction;
		//if (n_threads > 0) n_threads_ = n_threads;
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
	void HnswSolModel::InsertNodeSingleThread(HnswSolutionNode* qnode,
		custom_fun::VisitedLazyList * visited_list,
		custom_fun::VisitedLazyList* distanceVisited,
		std::vector<double>& distanceMemory)
	{
		using namespace std;
		int cur_level = qnode->GetLevel();
		int max_level_copy = 0;
		{
		//	std::shared_lock lock(max_level_guard_);
			max_level_copy = max_level_;
		}


		vector<HnswSolutionNode*> enterpoints;
		//std::vector<EdgeInfo> neighbors;
		if (cur_level < max_level_copy) {
			HnswSolutionNode* cur_node = enterpoint_;
			
			double d = getDistance(qnode, cur_node, m_env.get(), distanceVisited, distanceMemory);
			//double d = qnode->getOriginData().variableDistance(cur_node->getOriginData(), m_env.get());
			double cur_dist = d;
			for (auto i = max_level_copy; i > cur_level; --i) {
				bool changed = true;
				while (changed) {
					changed = false;
				//	unique_lock<mutex> local_lock(cur_node->GetAccessGuard());
					const auto& neighbors = cur_node->getFriends(i);
					//cur_node->getFriends(neighbors, i);
					//auto& neighbors = cur_node->getFriends(i);
					for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
						//d = dist_func_(qnode, *iter, data_dim_);

						d = getDistance(qnode, iter->GetNode(), m_env.get(), distanceVisited, distanceMemory);
						//d = qnode->getOriginData().variableDistance(iter->GetNode()->getOriginData(), m_env.get());

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
			SearchAtLayerSingleThread(qnode, enterpoints, i, visited_list, distanceVisited,distanceMemory, result);

			enterpoints.clear();
			priority_queue<FurtherSolFirst> next_enterpoints = result;
			while (next_enterpoints.size() > 0) {
				auto* top_node = next_enterpoints.top().GetNode();
				next_enterpoints.pop();
				enterpoints.push_back(top_node);
			}

			selecting_policy_->Select(m_,  i == 0, m_env.get(), result);

			


			while (result.size() > 0) {
				auto* top_node = result.top().GetNode();
		
				LinkSingleThread(top_node, qnode,result.top().GetDistance(), i);
				LinkSingleThread(qnode, top_node, result.top().GetDistance(), i);
				result.pop();
			}
		}

		if (cur_level > enterpoint_->GetLevel()) {

			{
			//	std::unique_lock lock(max_level_guard_);
				enterpoint_ = qnode;
				//max_level_ = cur_level;
			}
		}
	}
	void HnswSolModel::SearchAtLayerSingleThread(HnswSolutionNode* qnode,
		const std::vector<HnswSolutionNode*>& enterpoints, 
		int level, custom_fun::VisitedLazyList* visited_list, 
		custom_fun::VisitedLazyList* distanceVisited,
		std::vector<double>& distanceMemory,
		std::priority_queue<FurtherSolFirst>& result)const
	{
		using namespace std;
		priority_queue<CloserSolFirst> candidates;
		visited_list->Reset();

		for (const auto& enterpoint : enterpoints) {
			visited_list->MarkAsVisited(enterpoint->GetId());
			//float d = dist_func_(qnode, enterpoint, data_dim_);
			//double d = qnode->getOriginData().variableDistance(enterpoint->getOriginData(), m_env.get());
			double d = getDistance(qnode, enterpoint, m_env.get(), distanceVisited, distanceMemory);
			result.emplace(enterpoint, d);
			candidates.emplace(enterpoint, d);
		}
		//std::vector<EdgeInfo> neighbors;
		while (!candidates.empty()) {
			const CloserSolFirst& candidate = candidates.top();
			float lower_bound = result.top().GetDistance();
			if (candidate.GetDistance() > lower_bound)
				break;

			HnswSolutionNode* candidate_node = candidate.GetNode();
		//	unique_lock<mutex> lock(candidate_node->GetAccessGuard());
			//const auto& neighbors = candidate_node->GetFriends(level);
			//candidate_node->getFriends(neighbors, level);
			const auto& neighbors = candidate_node->getFriends(level);
			candidates.pop();
			for (const auto& neighbor : neighbors) {
				int id = neighbor.GetNode()->GetId();
				if (visited_list->NotVisited(id)) {
					visited_list->MarkAsVisited(id);
					

					double d = getDistance(qnode, neighbor.GetNode(), m_env.get(), distanceVisited, distanceMemory);
					//double d = qnode->getOriginData().variableDistance(neighbor.GetNode()->getOriginData(), m_env.get());
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
	
	void HnswSolModel::SearchAtLayerMultiThread(HnswSolutionNode* qnode,
		const std::vector<HnswSolutionNode*>& enterpoints,
		int level, custom_fun::VisitedLazyList* visited_list,
		custom_fun::VisitedLazyList* distanceVisited,
		std::vector<double>& distanceMemory,
		std::priority_queue<FurtherSolFirst>& result)const
	{
		using namespace std;
		priority_queue<CloserSolFirst> candidates;
		visited_list->Reset();

		for (const auto& enterpoint : enterpoints) {
			visited_list->MarkAsVisited(enterpoint->GetId());
			//float d = dist_func_(qnode, enterpoint, data_dim_);
			double d = getDistance(qnode, enterpoint, m_env.get(), distanceVisited, distanceMemory);
		//	double d = qnode->getOriginData().variableDistance(enterpoint->getOriginData(), m_env.get());
			result.emplace(enterpoint, d);
			candidates.emplace(enterpoint, d);
		}
		std::vector<EdgeInfo> neighbors;
		while (!candidates.empty()) {
			const CloserSolFirst& candidate = candidates.top();
			float lower_bound = result.top().GetDistance();
			if (candidate.GetDistance() > lower_bound)
				break;

			HnswSolutionNode* candidate_node = candidate.GetNode();
			//const auto& neighbors = candidate_node->getFriends(level);
			candidate_node->getFriends(neighbors, level);
			candidates.pop();
			for (const auto& neighbor : neighbors) {
				int id = neighbor.GetNode()->GetId();
				if (visited_list->NotVisited(id)) {
					visited_list->MarkAsVisited(id);

					double d = getDistance(qnode, neighbor.GetNode(), m_env.get(), distanceVisited, distanceMemory);
					//double d = qnode->getOriginData().variableDistance(neighbor.GetNode()->getOriginData(), m_env.get());
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

		source->link(target, dis, level, is_naive_, selecting_policy_.get(), m_env.get());
	}

	void HnswSolModel::LinkSingleThread(HnswSolutionNode* source, HnswSolutionNode* target, double dis, int level)
	{
		source->linkSingleThread(target, dis, level,  is_naive_, selecting_policy_.get(), m_env.get());
	}
	void HnswSolModel::MergeEdgesOfTwoGraphsSingleThread(
		HnswSolutionNode* cur, int level)
	{
		//if(cur->GetId()
		int id = cur->GetId();

		const auto& neighbors1 = nodes_[id]->getFriends(level);
		const auto& neighbors2= cur->getFriends(level);

		//std::vector<EdgeInfo> neighbors1, neighbors2;
		//nodes_[id]->getFriends(neighbors1,0);
		//cur->getFriends(neighbors2, 0);
		//const auto& neighbors1 = nodes_[id]->GetFriends(0);
		//const auto& neighbors2 = cur->GetFriends(0);
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
		}

		// Post Heuristic
		post_selecting_policy_->Select(max_m0_,  true, m_env.get(),temp_res);
		std::vector<EdgeInfo> merged_neighbors;
		while (!temp_res.empty()) {
			merged_neighbors.emplace_back(temp_res.top());
			temp_res.pop();
		}
		nodes_[id]->setFriendsSingleThread(level, merged_neighbors);;
	}
	void HnswSolModel::addDataTask(
		std::unique_ptr<ofec::GeneralMultiThreadInfo>& curThreadInfo,
		std::unique_ptr<ofec::GeneralMultiThreadInfo>& threadLocalInfo,
		std::unique_ptr<ofec::GeneralMultiThreadInfo>& globalThreadInfo) {
		using namespace std;


		auto& localInfo = dynamic_cast<LocalThreadInfo&>(*curThreadInfo);
		auto& threadInfo = dynamic_cast<ThreadLocalInfo&>(*threadLocalInfo);
		auto& globalInfo = dynamic_cast<GlobalThreadInfo&>(*globalThreadInfo);


		vector<HnswSolutionNode*> calculatingNodes;
		custom_fun::VisitedLazyList& curVisited = threadInfo.m_visisted;
		custom_fun::VisitedLazyList& distanceVisited = threadInfo.m_distanceVisited;
		std::vector<double>& distanceMemory = threadInfo.m_distanceMat;

		auto qnode = nodes_[localInfo.m_curNodeId];

		

		{
			unique_lock<mutex> lock(globalInfo.m_info_mtx);
			for (auto& it : globalInfo.m_calculatingNodes) {
				calculatingNodes.push_back(it);
			}
			globalInfo.m_calculatingNodes.insert(qnode);
		}

		distanceVisited.Reset();

		{

			int cur_level = qnode->GetLevel();
			int max_level_copy = 0;
			HnswSolutionNode* enterpoint_copy = nullptr;
			{
				std::shared_lock lock(max_level_guard_);
				max_level_copy = max_level_;
				enterpoint_copy = enterpoint_;
			}


			vector<HnswSolutionNode*> enterpoints;
			std::vector<EdgeInfo> neighbors;


			if (cur_level < max_level_copy) {
				HnswSolutionNode* cur_node = enterpoint_copy;
				double d = getDistance(qnode, cur_node, m_env.get(), &distanceVisited, distanceMemory);
				//double d = qnode->getOriginData().variableDistance(cur_node->getOriginData(), m_env.get());
				double cur_dist = d;
				for (auto i = max_level_copy; i > cur_level; --i) {
					bool changed = true;
					while (changed) {
						changed = false;
						cur_node->getFriends(neighbors, i);
						for (auto iter = neighbors.begin(); iter != neighbors.end(); ++iter) {
							double d = getDistance(qnode, iter->GetNode(), m_env.get(), &distanceVisited, distanceMemory);
							//	d = qnode->getOriginData().variableDistance(iter->GetNode()->getOriginData(), m_env.get());

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
				enterpoints.push_back(enterpoint_copy);
			}
			std::vector<EdgeInfo> beforeNeighbors;
			std::vector<std::vector<EdgeInfo>> level_neighbors(cur_level+1);
			for (auto i = cur_level; i >= 0; --i) {
				priority_queue<FurtherSolFirst> result;

				if (i <= max_level_copy) {
					SearchAtLayerMultiThread(qnode, enterpoints, i, &curVisited, &distanceVisited, distanceMemory, result);
					enterpoints.clear();

					priority_queue<FurtherSolFirst> next_enterpoints = result;
					while (next_enterpoints.size() > 0) {
						auto* top_node = next_enterpoints.top().GetNode();
						next_enterpoints.pop();
						enterpoints.push_back(top_node);
					}

				}


				{
					for (auto& otherNode : calculatingNodes) {
						if (curVisited.NotVisited(otherNode->GetId()) 
							&& otherNode->GetLevel() >= i) {
							double d = getDistance(qnode, otherNode, m_env.get(), &distanceVisited, distanceMemory);
							if (result.size() < ef_construction_ || result.top().GetDistance() > d) {
								result.emplace(otherNode, d);
								if (result.size() > ef_construction_)
									result.pop();
							}
							curVisited.MarkAsVisited(otherNode->GetId());

						}
					}
				}
				{
					std::unique_lock lock(qnode->mutex());
					beforeNeighbors = std::move(qnode->getFriends(i));
				}
				for (auto& it : beforeNeighbors) {
					if (curVisited.NotVisited(it.GetNode()->GetId())){
						result.emplace(it.GetNode(), it.GetDistance());
						curVisited.MarkAsVisited(it.GetNode()->GetId());
					}

				}
		
				selecting_policy_->Select(m_, i == 0, m_env.get(), result);

				neighbors.clear();
				{
					std::unique_lock lock(qnode->mutex());
					for (auto& it : qnode->getFriends(i)) {
						if (curVisited.NotVisited(it.GetNode()->GetId())) {
							result.emplace(it.GetNode(), it.GetDistance());
							curVisited.MarkAsVisited(it.GetNode()->GetId());
						}
						//result.emplace(it.GetNode(), it.GetDistance());
					}
					selecting_policy_->Select(m_, i == 0, m_env.get(), result);
					{
						while (result.size() > 0) {
							auto topNode = result.top();
							neighbors.emplace_back(topNode.GetNode(), topNode.GetDistance());
							result.pop();
						}
					}
					qnode->setFriendsSingleThread(i, neighbors);
				}

				level_neighbors[i] = std::move(neighbors);

			}
			// update network
			{
				curVisited.Reset();
				for (auto& it : calculatingNodes) {
					curVisited.MarkAsVisited(it->GetId());
				}
				for (auto i = min(max_level_copy, cur_level); i >= 0; --i) {
					auto& neighbors = level_neighbors[i];
					for (auto& it : neighbors) {
						if (it.GetNode()->GetId() < qnode->GetId()|| curVisited.Visited(it.GetNode()->GetId())) {
							Link(it.GetNode(), qnode, it.GetDistance(), i);
						}
					}
				}
			}

			if (cur_level > max_level_copy) {

				{
					std::unique_lock lock(max_level_guard_);
					if (cur_level > max_level_) {
						enterpoint_ = qnode;
						max_level_ = cur_level;
					}
	
				}
			}
		}

		{
			unique_lock<mutex> lock(globalInfo.m_info_mtx);
			globalInfo.m_calculatingNodes.erase(qnode);
			//curVisited = globalInfo.m_visitedLists.front();
		}
	}
}