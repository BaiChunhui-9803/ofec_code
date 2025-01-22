
#ifndef LOCAL_OPTIMA_NETWORK_H
#define LOCAL_OPTIMA_NETWORK_H
#include<set>
#include<vector>
#include<string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <sstream>
#include <thread>
#include <chrono>


class lon {
public:

	struct LonInfo {
		std::vector<int> m_nodeId;
		// node related, idx==-1 gray
		std::vector<int> m_node_funnel_idxs;
		std::vector<int> m_node_funnel_fit;
		std::vector<int> m_node_size;
		std::vector<int> m_node_fit;
		std::vector<std::vector<int>> m_graph;

		void output(const std::string& filepath, const std::string& filename) {

			std::ofstream out(filepath);
			out << filename << std::endl;
			out << m_graph.size() << std::endl;
			for (auto& it : m_graph) {

				out << it.front() << "\t" << it.back() << std::endl;
			}
			out << m_node_funnel_idxs.size() << std::endl;
			for (auto& it : m_node_funnel_idxs) {
				out << it << "\t";
			}
			out << std::endl;
			out << m_node_size.size() << std::endl;
			for (auto& it : m_node_size) {
				out << it << "\t";
			}
			out << std::endl;
			out << m_node_fit.size() << std::endl;
			for (auto& it : m_node_fit) {
				out << it << "\t";
			}
			out << std::endl;
			out.close();
		}


		void sortFit() {
			std::set<int> fitnessSet;
			for (int idx(0); idx < m_node_fit.size(); ++idx) {
				//std::set<int> fitnessSet;
				fitnessSet.insert(m_node_fit[idx]);
			}
			std::vector<int> sortedFit(fitnessSet.size());
			for (auto& it : fitnessSet) {
				sortedFit.push_back(it);
			}
			std::sort(sortedFit.begin(), sortedFit.end());
			int idFit = 1;
			std::map<int, int> fit2fitId;
			for (int idx(0); idx < sortedFit.size(); ++idx) {
				fit2fitId[sortedFit[idx]] = idx + 1;
			}

			for (auto& it :m_node_fit) {
				it = fit2fitId[it];
			}
		}
	};
	struct TraceEdge {
		int run;
		int iter;
		int id_from;
		int id_to;
		int num_kick;
		int count;
	};

	struct LonStatitic {
		//std::string m_tspname;
		//int LonType = 0;
		int m_numFunnel = 0;
		int m_numNode = 0;
		int m_numEdge = 0;

		void transfer(const LonInfo& lon) {
			m_numFunnel = 0;
			for (auto& it : lon.m_node_funnel_idxs) {
				m_numFunnel = std::max(it, m_numFunnel);
			}
			++m_numFunnel;
			m_numNode = lon.m_nodeId.size();
			m_numEdge = lon.m_graph.size();
		}



	};

	struct nodeInfo {
		int m_nodeId = 0;
		int m_sampleID = 0;
		int m_OptId = 0;
		long long cost = 0;
		bool is_exist = false;
		std::vector<int> m_sol;

		void outputNodeInfoHeader(std::ostream& out) {
			out << "Id\tSampleId\tCost\tOptId" << std::endl;
		}
		void outputNodeInfo(std::ostream& out) {
			out << m_nodeId << "\t" << m_sampleID << "\t" << cost << "\t";
			out << m_OptId << std::endl;
		}
		std::istream& inputNodeInfo(std::istream& in) {
			return in >> m_nodeId >> m_sampleID >> cost >> m_OptId;
		}
		static void inputHeader(std::istream& in) {
			std::string line;
			std::getline(in, line);
		}

		void outputSolHeader(std::ostream& out) {
			out << "Id\tSolution" << std::endl;
		}
		void outputSol(std::ostream& out) {
			out << m_nodeId;
			for (int idx(0); idx < m_sol.size(); ++idx) {
				out << "\t" << m_sol[idx];
			}
			out << std::endl;
		}
		std::istream& inputSol(std::istream& in) {
			std::string line;
			std::getline(in, line);
			std::stringstream ss(line);

			int number;
			m_sol.clear();
			ss >> m_nodeId;
			while (ss >> number) {
				m_sol.push_back(number);
			}
			return in;
		}
	};





	static void outputSolHeader(std::ostream& out) {
		out << "Id\tSolution" << std::endl;
	}
	static void outputSol(std::ostream& out,int nodeId, const std::vector<int>& sol) {
		out << nodeId;
		for (int idx(0); idx < sol.size(); ++idx) {
			out << "\t" << sol[idx];
		}
		out << std::endl;
	}

	static void filterNetworkRemoveLoopWishSink(
		std::set<int> sinkNodes,
		const std::vector<TraceEdge>& edges,
		const std::vector<double>& nodefit,
		LonInfo& lon);


	static void filterNetworkRemoveLoopWishSink(
		std::set<int> sinkNodes,
		const std::vector<TraceEdge>& edges,
		const std::vector<int>& nodefit,
		LonInfo& lon);

	static void filterLON(LonInfo& lon, LonInfo& lonfileter, double ratio) {
		std::vector<int> originId2newId;
		filterLON(lon, lonfileter, ratio, originId2newId);
	}
	static void filterLON(LonInfo& lon, LonInfo& lonfileter, double ratio, std::vector<int>& originId2newId);

	static void getFunnelIdx(std::set<int>& funnelSetIds,
		std::vector<TraceEdge>& trace
	);
	static void getSinkNode3(std::set<int> funnelSetIds, const std::vector<TraceEdge>& edges,
		const std::vector<double>& nodefit,
		std::set<int>& sinkNodes);

	static void getSinkNode3(std::set<int> funnelSetIds, const std::vector<TraceEdge>& edges,
		const std::vector<int>& nodefit,
		std::set<int>& sinkNodes);
	static void readTrace(const std::string& filepath, std::vector<TraceEdge>& trace);
	static void readNodeInfo(const std::string& filepath, std::vector<double>& nodefit);
	static void readEdge(const std::string& filepath, std::vector<TraceEdge>& edges);

	static void transfer(const std::vector<TraceEdge>& trace, std::vector<TraceEdge>& edges);
	////extern void transfer(const std::string& filepath, )

	static void outputTrace(const std::string& filepath, std::vector<TraceEdge>& trace);
	static void outputNodeInfo(const std::string& filepath, std::vector<int>& nodefit);
	//static void outputNodeInfoHeader(std::ofstream&out);
	//static void outputNodeInfoCur(std::ofstream&out, std::vector<int>& nodefit);
	static void outputEdge(const std::string& filepath, std::vector<TraceEdge>& edges);
	static void outputNodeSol(const std::string& filepath, std::vector<std::vector<int>>& nodeSol);
	static void outputLonNodeInfo(const std::string& filepath,const LonInfo& lon);
	static void outputLonEdgeInfo(const std::string& filepath, const LonInfo& lon);



	static void insertLonNodeInfo(const std::string& filepath, LonInfo& lon);
	static void insertLonEdgeInfo(const std::string& filepath, LonInfo& lon);



};

namespace name_tsp {
	extern void readSol(std::stringstream& buffer,
		std::vector<int>& solId,
		std::vector<std::vector<int>>& sols);
}
#endif
