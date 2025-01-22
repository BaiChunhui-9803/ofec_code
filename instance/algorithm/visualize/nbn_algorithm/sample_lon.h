#ifndef SAMPLELON_H
#define SAMPLELON_H

#include <vector>
#include <queue>
#include <core/algorithm/Solution.h>
#include <utility/function/custom_function.h>

namespace ofec_demo {
	struct SampleLON {
	//	using TSolution = Solution<>;
		static const int mc_num_sample;
		static const int mc_num_exploration;
		static const int mc_num_exploration_sample;
		static const double mc_num_exploration_radius;
		struct NetworkNearNodeInfo {
			int m_id = 0;
			double m_value = 0;
			int m_updated_times = 0;
			
			NetworkNearNodeInfo() = default;
			NetworkNearNodeInfo(int update_time) :m_updated_times(update_time) {}
			void set(int id, double value, int update_times) {
				m_id = id;
				m_value = value;
				m_updated_times = update_times;
			}

			void set(int id, double value) {
				m_id = id;
				m_value = value;
			}
			friend bool operator==(const NetworkNearNodeInfo& lhs, const NetworkNearNodeInfo& rhs)
			{
				return lhs.m_id== rhs.m_id;
			}
		};

		/// <summary>
		///  to be improve 
		/// </summary>
		struct NetworkNearNodeQue {
			std::vector<NetworkNearNodeInfo> m_nodes;
			int m_max_size;
		};
		struct ComSmaller {
			bool operator()(const NetworkNearNodeInfo& a, const NetworkNearNodeInfo& b) {
				if (a.m_value == b.m_value) {
					return a.m_updated_times > b.m_updated_times;
				}
				else return a.m_value < b.m_value;
			}
		};
		struct ComGreater {  
			bool operator()(const NetworkNearNodeInfo& a, const NetworkNearNodeInfo& b) {
				return a.m_value > b.m_value;
			}
		};

		struct NNnetworkQueue {

			static const int mc_queue_size;
			std::priority_queue<NetworkNearNodeInfo, std::vector<NetworkNearNodeInfo>, ComSmaller> m_neighbors;
			int m_updateTime = -1;
			bool pushNode(const NetworkNearNodeInfo& node, int updateTime) {
				bool ret_flag(true);
				m_neighbors.push(node);
				if (m_neighbors.size() > mc_queue_size) {
					if (m_neighbors.top() == node) {
						ret_flag = false;
					}
					m_neighbors.pop();
				}
				if (ret_flag) {
					m_updateTime = updateTime;
				}
				return ret_flag;
			}

			void getAllInfo(std::vector<NetworkNearNodeInfo>& info) {
				info.clear();
				std::priority_queue<NetworkNearNodeInfo, std::vector<NetworkNearNodeInfo>, ComSmaller> temp_neighbors;

				while (!m_neighbors.empty()) {
					info.push_back(m_neighbors.top());
					temp_neighbors.push(m_neighbors.top());
					m_neighbors.pop();
				}
				swap(temp_neighbors, m_neighbors);
			}
		};



		struct NNnetworkNode {

			double m_search_radiu = 1.0;
			double m_min_distance = 1e9;
			int m_belong_idx = -1;
		

			NNnetworkQueue m_nearest_network;
			NNnetworkQueue m_better_network;

			void updateInfo(const NetworkNearNodeInfo & nodeInfo) {
				if (m_min_distance > nodeInfo.m_value) {
					m_min_distance = nodeInfo.m_value;
					m_belong_idx = nodeInfo.m_id;
				}
			}
		};

		template<typename TSolution = ofec::Solution<>>
		struct DivideGroupInfo {
			int from_idx = 0;
			int updateTime = 0;
			std::vector<std::shared_ptr<TSolution>> samples;
			std::vector<NNnetworkNode> networks;
			std::pair<int, int> update_idxs_from_to;
			std::vector<int> update_centers_id;
			std::vector<std::vector<int>> group_idxs;
			std::vector<int> group_belong;
			int center_belong_id = -1;
			

			void initGroupInfo() {
				++updateTime;
				from_idx = samples.size();
				update_idxs_from_to.first = 0;
				update_idxs_from_to.second = 0;
				update_centers_id.clear();
				group_idxs.clear();
				group_belong.clear();
				center_belong_id = -1;
			}

			void addSols(const std::vector<std::shared_ptr<ofec::SolutionBase>>& sols) {
				int groupSize = sols.size();
				std::pair<int, int> from_to = { update_idxs_from_to.second,update_idxs_from_to.second + groupSize };
				update_idxs_from_to.second += groupSize;
				std::pair<int, int> sample_from_to;
				sample_from_to.first = samples.size();
				samples.resize(samples.size() + groupSize);
				sample_from_to.second = samples.size();
				for (int idx(sample_from_to.first); idx < sample_from_to.second; ++idx) {
					samples[idx].reset(new TSolution());
					*samples[idx] = *sols[idx];
				}
				networks.resize(networks.size() + groupSize);
				auto& cur_gidxs(group_idxs.back());
				
				for (int idx(sample_from_to.first); idx < sample_from_to.second; ++idx) {
					cur_gidxs.push_back(idx);
				}
				group_belong.resize(from_to.second);
				std::fill(group_belong.begin() + from_to.first,
					group_belong.begin() + from_to.second,
					center_belong_id
				);
			}

			void addGroupInfo(int groupSize,int center_idx=-1) {
				std::pair<int, int> from_to = { update_idxs_from_to.second,update_idxs_from_to.second+groupSize};
				update_idxs_from_to.second += groupSize;
				std::pair<int, int> sample_from_to;
				sample_from_to.first = samples.size();
				samples.resize(samples.size() + groupSize);
				sample_from_to.second = samples.size();
				networks.resize(networks.size() + groupSize);
				update_centers_id.push_back(center_idx);
				group_idxs.push_back(std::vector<int>());
				auto& cur_gidxs(group_idxs.back());
				for (int idx(sample_from_to.first); idx < sample_from_to.second; ++idx) {
					cur_gidxs.push_back(idx);
				}
				if (center_idx != -1) {
					cur_gidxs.push_back(center_idx);
				}
				group_belong.resize(from_to.second);
				++center_belong_id;
				std::fill(group_belong.begin() + from_to.first, group_belong.begin() + from_to.second, 
					center_belong_id
					);
			}

			
		};


		static void toIdx(int idx, std::pair<int, int>& idxs);
		static void toId(const std::pair<int, int>& idxs, int& id);


		template<typename TSolution = ofec::Solution<>>
		static void sample_sols_NN_divide(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			std::vector<int>& belong_idx,
			Problem *pro, Random *rnd);

		template<typename TSolution = ofec::Solution<>>
		static void sample_sols_NN(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			std::vector<int>& belong_idx,
			Problem *pro, Random *rnd);
		template<typename TSolution = ofec::Solution<>>
		static void sample_sols_NN_threads(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, 
			std::vector<int>& belong_idx,
			Problem *pro, Random *rnd);
		template<typename TSolution = ofec::Solution<>>
		static void getSampleResult(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			std::vector<int>& belong_idxs,
			std::vector<std::shared_ptr<TSolution>> &samples,
			const std::vector<NNnetworkNode>& newtork
			);
		

		template<typename TSolution = ofec::Solution<>>
		static void sample_sols(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			Problem *pro, Random *rnd);	


		template<typename TSolution = ofec::Solution<>>
		static void sample_sols(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			std::vector<int>&belong_idx,
			Problem *pro, Random *rnd);


		template<typename TSolution= ofec::Solution<>>
		static void setFitness(ofec::SolutionBase& sol, Problem *pro);

		template<typename TSolution = ofec::Solution<>>
		static void generate_sols(std::shared_ptr<TSolution>& sol,
			int sol_id, Problem *pro , Random *rnd,
			const std::shared_ptr<TSolution>& center, 
			double radius);



		template<typename TSolution = ofec::Solution<>>
		static void generate_sols_threadTask(std::shared_ptr<TSolution>& sol,
			int sol_id, Problem *pro, Random *rnd,
			const std::shared_ptr<TSolution>& center,
			double radius);

		template<typename TSolution = ofec::Solution<>>
		static void generate_sols(std::shared_ptr<TSolution>& sol,
			const std::vector<std::shared_ptr<TSolution>>& update_centers,
		    const std::vector<double>& update_radius,
			const std::pair<int,int>& from_to,
			int from_idx,
			Problem *pro, Random *rnd);

		template<typename TSolution = ofec::Solution<>>
		static void generate_sols_threadTask_newly(
			DivideGroupInfo<TSolution>& group_info,
			const std::pair<int, int>& from_to,
			Problem *pro, Random *rnd);
		




		//template<typename TSolution = Solution<>>
		//static void generate_sols(std::shared_ptr<TSolution>& sol,
		//	int sol_id, Problem *pro, Random& random,
		//	const std::shared_ptr<TSolution>& center,
		//	double radius);
		template<typename TSolution = ofec::Solution<>>
		static void create_thread_generateSols(
			DivideGroupInfo<TSolution>& group_info,
			const std::vector<int>& temp_id_rnd,
			Problem *pro);

		template<typename TSolution = ofec::Solution<>>
		static void create_thread_updateNetworks(
			DivideGroupInfo<TSolution>& group_info,
			Problem *pro);


		template<typename TSolution = ofec::Solution<>>
		static void generateNNetwork(
			std::vector<std::shared_ptr<TSolution>>& sols,
			const std::vector<int>& idxs,
			std::vector<NNnetworkNode>& newtork,
			Problem *pro,
			int udpateTime = 0
		);


		template<typename TSolution = ofec::Solution<>>
		static void generateNNetwork_ThreadTask(
			DivideGroupInfo<TSolution>& groupInfo,
			const std::pair<int, int>& from_to,
			Problem *pro
		);


		


		template<typename TSolution = ofec::Solution<>>
		static void generateNNetwork(
			DivideGroupInfo<TSolution>& groupInfo,
			const std::pair<int, int>& from_to,
			Problem *pro
		);


		template<typename TSolution = ofec::Solution<>>
		static void generateNNetworkThreadTask(
			std::vector<std::shared_ptr<TSolution>>& sols,
			const std::vector<int>& idxs,
			const std::pair<int, int>& from_to,
			std::vector<NNnetworkNode>& newtork,
			Problem *pro,
			int udpateTime = 0
		);


		template<typename TSolution = ofec::Solution<>>
		static void generateNNetworkThreadTask(
			std::vector<std::shared_ptr<TSolution>>& sols,
			std::vector<NNnetworkNode>& newtork,
			const std::pair<int,int>& from_to,
			int from_idx,Problem *pro, int udpateTime = 0
		);


		template<typename TSolution = ofec::Solution<>>
		static void updateNNetworkInfo(
			std::vector<std::shared_ptr<TSolution>>& sols,
			std::vector<NNnetworkNode>& newtork,
			int sol_idx,
			int updateTime,
			Problem *pro
		);




	};


	/*
	* 		template<typename TSolution = ofec::Solution<>>
		static void sample_sols_NN(
			std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
			const std::vector<int>& belong_idx,
			Problem *pro, Random *rnd);
	*/


	template<typename TSolution>
	void SampleLON::sample_sols_NN(std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
		std::vector<int>& belong_idx,
		Problem *pro, Random *rnd) {
		std::vector<std::shared_ptr<TSolution>> samples;
		std::vector<NNnetworkNode> networks;
		std::vector<int> sort_sidx;
		std::vector<int> center_idxs;
		int updateTime = 0;
		std::vector<int> idxs;
		std::shared_ptr<TSolution> center = nullptr;
		double radius = 1.0;
		{
			samples.resize(mc_num_sample);
			networks.resize(mc_num_sample);
			idxs.resize(samples.size());
			for (int idx(0); idx < samples.size(); ++idx) {
				generate_sols<TSolution>(samples[idx], idx, pro, rnd, center, radius);
				idxs[idx] = idx;
			}
		}

		generateNNetwork<TSolution>(samples, idxs, networks, pro, updateTime);

		bool updateFlag(true);
		while (updateFlag) {

			center_idxs.clear();
			updateFlag = false;
			sort_sidx.resize(samples.size());
			for (int idx(0); idx < sort_sidx.size(); ++idx) {
				sort_sidx[idx] = idx;
			}
			std::sort(sort_sidx.begin(), sort_sidx.end(), [&](int a, int b) {
				return networks[a].m_min_distance > networks[b].m_min_distance;
			});

			std::priority_queue<NetworkNearNodeInfo, std::vector<NetworkNearNodeInfo>, ComGreater> centers;

			int cur_idx(0);
			NetworkNearNodeInfo node(updateTime);
			while (centers.size() < mc_num_exploration || centers.top().m_value < networks[sort_sidx[cur_idx]].m_min_distance) {
				auto& cur_id(sort_sidx[cur_idx]);

				updateNNetworkInfo<TSolution>(samples, networks, cur_id, updateTime, pro);
				node.set(cur_id, networks[cur_id].m_min_distance);
				centers.push(node);
				if (centers.size() > mc_num_exploration) centers.pop();
				++cur_idx;
			}
			center_idxs.clear();
			while (!centers.empty()) {
				center_idxs.push_back(centers.top().m_id);
				centers.pop();
			}

			++updateTime;
			while (!center_idxs.empty()) {

				if (networks[center_idxs.back()].m_search_radiu > 0.1) {
					updateFlag = true;
					idxs.clear();
					idxs.push_back(center_idxs.back());
					center = samples[center_idxs.back()];


					int from(samples.size());
					samples.resize(samples.size() + mc_num_exploration_sample);
					int to(samples.size());


					networks.resize(samples.size());

					for (int idx(from); idx < to; ++idx) {
						generate_sols<TSolution>(samples[idx],idx, pro, rnd, center, networks[center->id()].m_search_radiu);
						idxs.push_back(idx);
					}
					generateNNetwork<TSolution>(samples, idxs, networks, pro, updateTime);

					networks[center->id()].m_search_radiu /= 2.0;
				}

				center_idxs.pop_back();

			}
		}


		getSampleResult<TSolution>(sols, belong_idx, samples,networks);

	}


	template<typename TSolution>
	void SampleLON::sample_sols_NN_divide(
		std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
		std::vector<int>& belong_idx,
		Problem *pro, Random *rnd) {


	}

	template<typename TSolution>
	void SampleLON::getSampleResult(
		std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
		std::vector<int>& belong_idx,
		std::vector<std::shared_ptr<TSolution>>& samples,
		const std::vector<NNnetworkNode>& networks
	) {
		sols.resize(samples.size());
		belong_idx.resize(samples.size());
		for (int idx(0); idx < sols.size(); ++idx) {
			if (networks[samples[idx]->id()].m_belong_idx == -1) {
				belong_idx[idx] = idx;
			}
			else {
				belong_idx[idx] = networks[samples[idx]->id()].m_belong_idx;
			}
			sols[idx] = std::move(samples[idx]);
		}
	}



	template<typename TSolution >
	void SampleLON::create_thread_generateSols(
		DivideGroupInfo<TSolution>& group_info,
		const std::vector<int>& temp_id_rnd,
		Problem *pro) {

		using namespace ofec;

		{

			std::pair<int, int> from_to = group_info.update_idxs_from_to;

			SampleLON::generate_sols_threadTask_newly<TSolution>(
				group_info,from_to,
				pro, temp_id_rnd[0]);
		}


			//{
			//	std::vector<std::thread> thrds;

			//	int num_samples = group_info.update_idxs_from_to.second - group_info.update_idxs_from_to.first;
			//	std::vector<int> tasks;

			//	int num_task = temp_id_rnd.size();
			//	UTILITY::assignThreads(num_samples, num_task, tasks);
			//	std::pair<int, int> from_to;
			//	std::vector<int> temp_rnd(num_task);
			//	for (size_t i = 0; i < num_task; ++i) {

			//		from_to.first = tasks[i];
			//		from_to.second = tasks[i + 1];


			//		thrds.push_back(std::thread(
			//			&SampleLON::generate_sols_threadTask_newly<TSolution>,
			//			std::ref(group_info),
			//			std::cref(from_to),
			//			pro, temp_id_rnd[i]));
			//	}
			//	for (auto& thrd : thrds)
			//		thrd.join();

			//}
	}

	template<typename TSolution >
	void SampleLON::create_thread_updateNetworks(
		DivideGroupInfo<TSolution>& group_info,
		Problem *pro) {
		using namespace ofec;
			{

				std::vector<std::thread> thrds;

				int num_task = std::thread::hardware_concurrency();
				int num_samples = group_info.update_idxs_from_to.second - group_info.update_idxs_from_to.first;
				std::vector<int> tasks;
				UTILITY::assignThreads(num_samples, num_task, tasks);
				std::pair<int, int> from_to;
				std::vector<int> temp_rnd(num_task);

				std::vector<std::pair<int, int>> from_to_record;
				for (size_t i = 0; i < num_task; ++i) {

					from_to.first = tasks[i];
					from_to.second = tasks[i + 1];

					from_to_record.push_back(from_to);
					thrds.push_back(std::thread(
						&SampleLON::generateNNetwork_ThreadTask<TSolution>,
						std::ref(group_info),
						std::cref(from_to),
						pro));
				}
				for (auto& thrd : thrds)
					thrd.join();


			}
	}


	template<typename TSolution>
	void SampleLON::sample_sols_NN_threads(
		std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
		std::vector<int>& belong_idx,
		Problem *pro, Random *rnd) {



		DivideGroupInfo<TSolution> group_info;
		group_info.initGroupInfo();
		group_info.addGroupInfo(mc_num_sample);


		int num_task = std::thread::hardware_concurrency();
		std::vector<int> temp_id_rnd(num_task);
		for (int idx(0); idx < num_task; ++idx) {
			temp_id_rnd[idx] = ofec::ADD_RND(double(idx+1) / double(num_task+1));
		}

	
		create_thread_generateSols<TSolution>(group_info, temp_id_rnd, pro);
		//auto from_to = group_info.updateGroupInfo(sols.size());
		group_info.addSols(sols);

		create_thread_updateNetworks<TSolution>(group_info, pro);



		std::vector<int> center_idxs;
		std::vector<int> sort_sidx;

	/*	bool updateFlag(true);
		while (updateFlag) {
			auto& samples(group_info.samples);
			auto& networks(group_info.networks);
			auto& udpateTime(group_info.updateTime);
			center_idxs.clear();
			updateFlag = false;
			sort_sidx.resize(samples.size());
			for (int idx(0); idx < sort_sidx.size(); ++idx) {
				sort_sidx[idx] = idx;
			}
			std::sort(sort_sidx.begin(), sort_sidx.end(), [&](int a, int b) {
				return networks[a].m_min_distance > networks[b].m_min_distance;
			});

			std::priority_queue<NetworkNearNodeInfo, std::vector<NetworkNearNodeInfo>, ComGreater> centers;

			int cur_idx(0);
			NetworkNearNodeInfo node(updateTime);
			while (centers.size() < mc_num_exploration || centers.top().m_value < networks[sort_sidx[cur_idx]].m_min_distance) {
				auto& cur_id(sort_sidx[cur_idx]);

				updateNNetworkInfo<TSolution>(samples, networks, cur_id, updateTime, pro);
				node.set(cur_id, networks[cur_id].m_min_distance);
				centers.push(node);
				if (centers.size() > mc_num_exploration) centers.pop();
				++cur_idx;
			}
			center_idxs.clear();
			while (!centers.empty()) {
				center_idxs.push_back(centers.top().m_id);
				centers.pop();
			}
			group_info.initGroupInfo();
			while (!center_idxs.empty()) {
				if (networks[center_idxs.back()].m_search_radiu > mc_num_exploration_radius) {
					group_info.addGroupInfo(mc_num_sample, center_idxs.back());					
				}
				center_idxs.pop_back();
			}

			create_thread_generateSols<TSolution>(group_info, temp_id_rnd, pro);
			create_thread_updateNetworks<TSolution>(group_info, pro);


			for (int idx(0); idx < group_info.samples.size(); ++idx) {
				updateNNetworkInfo<TSolution>(group_info.samples, group_info.networks,
					idx, group_info.updateTime, pro);
			}
		}*/

		for (int idx(0); idx < num_task; ++idx) {
			ofec::DEL_RND(temp_id_rnd[idx]);
		}

	//	getSampleResult<TSolution>(sols,belong_idx,)

		getSampleResult<TSolution>(sols, belong_idx, group_info.samples,group_info.networks);
		
		/*{
			auto& samples(group_info.samples);
			auto& networks(group_info.networks);
			sols.resize(samples.size());
			belong_idx.resize(samples.size());
			for (int idx(0); idx < sols.size(); ++idx) {
				if (networks[samples[idx]->id()].m_belong_idx == -1) {
					belong_idx[idx] = idx;
				}
				else {
					belong_idx[idx] = networks[samples[idx]->id()].m_belong_idx;
				}
				sols[idx] = std::move(samples[idx]);
			}
		}*/
		

	}

	



	template<typename TSolution>
	void SampleLON::sample_sols(std::vector<std::shared_ptr<ofec::SolutionBase>>& sols, Problem *pro, Random *rnd) {
		std::vector<std::shared_ptr<TSolution>> samples;
		std::shared_ptr<TSolution> center = nullptr;
		double radius = 1.0;
		int total_idx(0);
		{
			samples.resize(mc_num_sample);
			for (int idx(0); idx < samples.size(); ++idx) {
				generate_sols<TSolution>(samples[idx], total_idx++, pro, rnd, center, radius);
			}
		}


		sols.resize(samples.size());
		for (int idx(0); idx < sols.size(); ++idx) {
			sols[idx] = std::move(samples[idx]);
		}
	}



	template<typename TSolution >
	void SampleLON::sample_sols(
		std::vector<std::shared_ptr<ofec::SolutionBase>>& sols,
		std::vector<int>& belong_idx,
		Problem *pro, Random *rnd) {
		std::vector<std::shared_ptr<TSolution>> samples;
		std::vector<NNnetworkNode> networks;
		std::vector<int> sort_sidx;
		std::vector<int> center_idxs;
		int updateTime = 0;
		std::vector<int> idxs;
		std::shared_ptr<TSolution> center = nullptr;
		double radius = 1.0;
		{
			samples.resize(mc_num_sample);
			networks.resize(mc_num_sample);
			idxs.resize(samples.size());
			for (int idx(0); idx < samples.size(); ++idx) {
				generate_sols<TSolution>(samples[idx],idx, pro, rnd, center, radius);
				idxs[idx] = idx;
			}
		}

		generateNNetwork<TSolution>(samples, idxs, networks, pro, updateTime);
		getSampleResult(sols, belong_idx, samples,networks);

	}

	template<typename TSolution>
	void generatorSol(
		std::vector<std::shared_ptr<TSolution>> &samples,
		std::vector<int>& idxs,
		const std::pair<int,int>& from_to,
		Problem *pro, Random *rnd,const std::shared_ptr<TSolution>& center,
		double radius
	) {
		{
			samples.resize(mc_num_sample);
			networks.resize(mc_num_sample);
			idxs.resize(samples.size());
			for (int idx(0); idx < samples.size(); ++idx) {
				generate_sols<TSolution>(samples[idx], total_idx++, pro, rnd, center, radius);
				idxs[idx] = idx;
			}


		}
	}



	template<typename TSolution>
	void SampleLON::setFitness(ofec::SolutionBase& sol, Problem *pro) {
		if (sol.objective().size() == 1)
		{
			if (ofec::pro->optimizeMode()[0] == ofec::OptimizeMode::kMaximize) {
				sol.setFitness(sol.objective()[0]);
			}
			else {
				sol.setFitness(-sol.objective()[0]);
			}
		}
		else {
			auto& indi(dynamic_cast<TSolution&>(sol));
			sol.setFitness(-indi.fitness());
		}
	}


	template<typename TSolution>
	void SampleLON::generate_sols(std::shared_ptr<TSolution>& sol,
		int sol_id, Problem *pro, Random *rnd,
		const std::shared_ptr<TSolution>& center,
		double radius)
	{
		auto& pro(ofec::pro);
		sol.reset(new TSolution(pro.numberObjectives(), pro.numberConstraints(), pro.numberVariables()));
		sol->setId(sol_id);
		if (center != nullptr) {
			pro.initializeSolution(*sol, *center, radius, rnd);
		}
		else {
			//pro.initializeSolution()
			pro.initializeSolution(*sol, rnd);
		}
		sol->evaluate(pro, -1, false);
		setFitness<TSolution>(*sol, pro);
	}

	template<typename TSolution>
	void SampleLON::generate_sols_threadTask(std::shared_ptr<TSolution>& sol,
		int sol_id, Problem *pro, Random *rnd,
		const std::shared_ptr<TSolution>& center,
		double radius) {
		auto& pro(ofec::pro);
		sol.reset(new TSolution(pro.numberObjectives(), pro.numberConstraints(), pro.numberVariables()));
		sol->setId(sol_id);
		if (center != nullptr) {
			pro.initializeSolution(*sol, *center, radius, rnd);
		}
		else {
			//pro.initializeSolution()
			pro.initializeSolution(*sol, rnd);
		}
		sol->evaluate(pro, -1, false);
		setFitness<TSolution>(*sol, pro);
	}



	template<typename TSolution>
	void SampleLON::generate_sols(std::shared_ptr<TSolution>& sol,
		const std::vector<std::shared_ptr<TSolution>>& update_centers,
		const std::vector<double>& update_radius,
		const std::pair<int, int>& from_to,
		int from_idx,
		Problem *pro, Random *rnd) {
		std::pair<int, int> cur_idxs;
		for (int idx(from_to.first); idx < from_to.second; ++idx) {
			toIdx(idx-from_idx, cur_idxs);
			generate_sols<TSolution>(sols[idx], idx, pro, rnd, update_centers[cur_idxs.first],
				update_radius[cur_idxs.first]);
		}
	}


	template<typename TSolution>
	void SampleLON::generate_sols_threadTask_newly(
		DivideGroupInfo<TSolution>& group_info,
		const std::pair<int, int>& from_to,
		Problem *pro, Random *rnd) {

		auto& sols(group_info.samples);
		std::shared_ptr<TSolution> center = nullptr;
		double radius = 1.0;
		for (int idx(from_to.first); idx < from_to.second; ++idx) {
			int sol_idx = idx + group_info.from_idx;
			int center_idx = group_info.group_belong[idx];
			if (center_idx == -1) {
				generate_sols<TSolution>(sols[sol_idx], sol_idx, pro, rnd,
					center,
					radius);
			}
			else {
				generate_sols<TSolution>(sols[sol_idx], sol_idx, pro, rnd,
					group_info.samples[center_idx],
					group_info.networks[center_idx].m_search_radiu);
			}

		}
	}

	template<typename TSolution>
	void SampleLON::generateNNetwork(
		std::vector<std::shared_ptr<TSolution>>& sols,
		const std::vector<int>& idxs,
		std::vector<NNnetworkNode>& network,
		Problem *pro,
		int updateTime
	) {

		NetworkNearNodeInfo nodeInfo_x(updateTime);
		NetworkNearNodeInfo nodeInfo_y(updateTime);
		for (int idx(0); idx < idxs.size(); ++idx) {
			int sol_idx = idxs[idx];
			auto& node_x(network[sols[sol_idx]->id()]);
			auto& solx(sols[sol_idx]);
			//node_x.m_updateTime = updateTime;
			for (int idy(0); idy < idx; ++idy) {
				int sol_idy = idxs[idy];
				auto& soly(sols[sol_idy]);
				double dis = solx->variableDistance(*soly, pro);

				nodeInfo_x.set(sol_idy, dis);
				nodeInfo_y.set(sol_idx, dis);
				auto& node_y(network[sols[sol_idy]->id()]);
				node_x.m_nearest_network.pushNode(nodeInfo_x, updateTime);
				node_y.m_nearest_network.pushNode(nodeInfo_y, updateTime);
				if (solx->fitness() < soly->fitness()) {
					if (node_x.m_better_network.pushNode(nodeInfo_x, updateTime)) {
						
						node_x.updateInfo(nodeInfo_x);
	/*					if (node_x.m_min_distance > nodeInfo_x.m_value) {
							node_x.m_min_distance = nodeInfo_x.m_value;
							node_x.m_belong_idx = nodeInfo_x.m_id;
						}*/
					//	node_x.m_min_distance = std::min(nodeInfo_x.m_value, node_x.m_min_distance);
					}
				}
				else if (solx->fitness() > soly->fitness()) {
					if (node_y.m_better_network.pushNode(nodeInfo_y, updateTime)) {
						node_y.updateInfo(nodeInfo_y);
						
						//node_y.m_min_distance = std::min(nodeInfo_y.m_value, node_y.m_min_distance);
					}


				}

			}
		}

	}


	template<typename TSolution>
	void SampleLON::generateNNetwork(
		DivideGroupInfo<TSolution>& groupInfo,
		const std::pair<int, int>& from_to,
		Problem *pro
	) {
		NetworkNearNodeInfo nodeInfo_x(groupInfo.updateTime);
		for (int idx(from_to.first); idx < from_to.second; ++idx) {
			int sol_idx = groupInfo.from_idx + idx;
			auto& node_x(groupInfo.networks[groupInfo.samples[sol_idx]->id()]);
			auto& solx(sols[sol_idx]);

			int belong_idx(groupInfo.group_belong[idx]);
			auto& idy_from_to(groupInfo.group_idxs[belong_idx]);
			for (const auto& idy : idy_from_to) {
				int sol_idy = idy + groupInfo.from_idx;
				auto& soly(groupInfo.samples[sol_idy]);
				double dis = solx->variableDistance(*soly, pro);

				nodeInfo_x.set(sol_idy, dis);
				auto& node_y(groupInfo.networks[sols[sol_idy]->id()]);
				node_x.m_nearest_network.pushNode(nodeInfo_x, updateTime);

				if (solx->fitness() < soly->fitness()) {
					if (node_x.m_better_network.pushNode(nodeInfo_x, updateTime)) {

						node_x.updateInfo(nodeInfo_x);
						//node_x.m_min_distance = std::min(nodeInfo_x.m_value, node_x.m_min_distance);
					}
				}
			}

		}
	}



	template<typename TSolution >
	void SampleLON::generateNNetwork_ThreadTask(
		DivideGroupInfo<TSolution>& groupInfo,
		const std::pair<int, int>& from_to,
		Problem *pro
	) {
		auto& sols(groupInfo.samples);
		auto& updateTime(groupInfo.updateTime);
		NetworkNearNodeInfo nodeInfo_x(groupInfo.updateTime);
		for (int idx(from_to.first); idx < from_to.second; ++idx) {
			int sol_idx = groupInfo.from_idx + idx;
			NNnetworkNode node_x;
			node_x= (groupInfo.networks[groupInfo.samples[sol_idx]->id()]);
			auto& solx(sols[sol_idx]);

			int belong_idx(groupInfo.group_belong[idx]);
			auto& idy_from_to(groupInfo.group_idxs[belong_idx]);
			for (const auto& idy : idy_from_to) {
				int sol_idy = idy + groupInfo.from_idx;
				auto& soly(groupInfo.samples[sol_idy]);
				double dis = solx->variableDistance(*soly, pro);

				nodeInfo_x.set(sol_idy, dis);
			//	auto& node_y(groupInfo.networks[sols[sol_idy]->id()]);
				node_x.m_nearest_network.pushNode(nodeInfo_x, updateTime);

				if (solx->fitness() < soly->fitness()) {
					if (node_x.m_better_network.pushNode(nodeInfo_x, updateTime)) {

						node_x.updateInfo(nodeInfo_x);
						//node_x.m_min_distance = std::min(nodeInfo_x.m_value, node_x.m_min_distance);
					}
				}
			}

			groupInfo.networks[groupInfo.samples[sol_idx]->id()] = node_x;
		}
	}

	template<typename TSolution>
	void SampleLON::generateNNetworkThreadTask(
		std::vector<std::shared_ptr<TSolution>>& sols,
		const std::vector<int>& idxs,
		const std::pair<int, int>& from_to,
		std::vector<NNnetworkNode>& newtork,
		Problem *pro,
		int udpateTime
	) {
		NetworkNearNodeInfo nodeInfo_x(updateTime);
		for (int idx(from_to.first); idx < from_to.second; ++idx) {
			int sol_idx = idxs[idx];
			auto& node_x(network[sols[sol_idx]->id()]);
			auto& solx(sols[sol_idx]);
			for (auto& sol_idy : idxs) {
				auto& soly(sols[sol_idy]);
				double dis = solx->variableDistance(*soly, pro);

				nodeInfo_x.set(sol_idy, dis);
				auto& node_y(network[sols[sol_idy]->id()]);
				node_x.m_nearest_network.pushNode(nodeInfo_x, updateTime);

				if (solx->fitness() < soly->fitness()) {
					if (node_x.m_better_network.pushNode(nodeInfo_x, updateTime)) {
						
						node_x.updateInfo(nodeInfo_x);
						//node_x.m_min_distance = std::min(nodeInfo_x.m_value, node_x.m_min_distance);
					}
				}
			}
		}
	}

	template<typename TSolution>
	void SampleLON::generateNNetworkThreadTask(
		std::vector<std::shared_ptr<TSolution>>& sols,
		std::vector<NNnetworkNode>& newtork,
		const std::pair<int, int>& from_to,
		int from_idx, Problem *pro, int udpateTime
	) {
		std::pair<int, int> solx_idxs;
		std::pair<int, int> soly_idxs;
		for (int sol_idx(from_to.first); sol_idx < from_to.second; ++sol_idx) {
			auto& node_x(network[sols[sol_idx]->id()]);
			auto& solx(sols[sol_idx]);
			toIdx(sol_idx-from_idx, solx_idxs);
			soly_idxs.first = solx_idxs.first;
			int sol_idy(0);

			for (int idy(0); idy < mc_num_exploration_sample; ++idy) {
				soly_idxs.second = idy;
				toId(soly_idxs, sol_idy);
				sol_idy += from_idx;

				auto& soly(sols[sol_idy]);
				double dis = solx->variableDistance(*soly, pro);

				nodeInfo_x.set(sol_idy, dis);
				auto& node_y(network[sols[sol_idy]->id()]);
				node_x.m_nearest_network.pushNode(nodeInfo_x, updateTime);

				if (solx->fitness() < soly->fitness()) {
					if (node_x.m_better_network.pushNode(nodeInfo_x, updateTime)) {

						node_x.updateInfo(nodeInfo_x);
						//node_x.m_min_distance = std::min(nodeInfo_x.m_value, node_x.m_min_distance);
					}
				}
			}
		}
	}


	template<typename TSolution>
	void SampleLON::updateNNetworkInfo(
		std::vector<std::shared_ptr<TSolution>>& sols,
		std::vector<NNnetworkNode>& network,
		int sol_idx, int updateTime, Problem *pro
	) {
		std::queue<int> better_idxs;
		std::queue<int> nearest_idxs;

		std::vector<NetworkNearNodeInfo> info;
		std::queue<int> newly_bqueue;
		std::queue<int> newly_nqueue;
		auto& cur_sol(sols[sol_idx]);
		auto& cur_net(network[sol_idx]);
		NetworkNearNodeInfo nodeInfo(updateTime);
		cur_net.m_better_network.getAllInfo(info);
		for (auto& it : info) {
			if (it.m_updated_times < network[it.m_id].m_nearest_network.m_updateTime) {
				better_idxs.push(it.m_id);
			}
		}

		cur_net.m_nearest_network.getAllInfo(info);

		for (auto& it : info) {
			if (it.m_updated_times < network[it.m_id].m_nearest_network.m_updateTime) {
				nearest_idxs.push(it.m_id);
			}
		}



		while (!better_idxs.empty() && !nearest_idxs.empty()) {
			while (!better_idxs.empty()) {
				int nei_idx = better_idxs.front();
				better_idxs.pop();
				network[nei_idx].m_nearest_network.getAllInfo(info);
				for (auto& it : info) {
					if (cur_sol->fitness() < sols[it.m_id]->fitness()) {
						double dis = cur_sol->variableDistance(*sols[it.m_id], pro);
						nodeInfo.set(it.m_id, dis);
						if (cur_net.m_better_network.pushNode(nodeInfo, updateTime)) {
							newly_bqueue.push(it.m_id);
							cur_net.updateInfo(nodeInfo);
						}
					}
				}
			}
			while (!nearest_idxs.empty()) {
				int nei_idx = nearest_idxs.front();
				nearest_idxs.pop();
				network[nei_idx].m_nearest_network.getAllInfo(info);
				for (auto& it : info) {
					double dis = cur_sol->variableDistance(*sols[it.m_id], pro);
					nodeInfo.set(it.m_id, dis);
					if (cur_sol->fitness() < sols[it.m_id]->fitness()) {
						if (cur_net.m_better_network.pushNode(nodeInfo, updateTime)) {
							newly_bqueue.push(it.m_id);
							cur_net.updateInfo(nodeInfo);
						}
					}
					if (cur_net.m_nearest_network.pushNode(nodeInfo, updateTime)) {
						newly_nqueue.push(it.m_id);
					}

				}
			}
			swap(better_idxs, newly_bqueue);
			swap(nearest_idxs, newly_nqueue);
		}
	}


}


#endif 