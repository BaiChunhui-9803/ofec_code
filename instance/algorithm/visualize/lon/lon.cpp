#include"lon.h"
#include <queue>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
//using namespace lon;




void convertStr(const std::string& info, size_t from, size_t to,
	std::vector<int>& solIds,
	std::vector<std::vector<int>>& sols) {
	std::string line;
	int solId(0);
	int total(0);
	//auto startTime = std::chrono::system_clock::now();
	size_t oldFrom = from;
	std::vector<int> cuts;

	while (from < to) {
		cuts.push_back(from);
		size_t last = info.find_first_of('\n', from);
		line = info.substr(from, last - from);
		std::stringstream ss;
		ss << line;
		//ss >> solId;
		if (!(ss >> solId)) return;
		int number(0);
		solIds[solId - 1] = solId;
		sols[solId - 1].clear();
		while (ss >> number) {
			sols[solId - 1].push_back(number);
		}

		from = last + 1;
		++total;
		//if (total % 10000 == 0) {
		//	auto endTime = std::chrono::system_clock::now();
		//	std::chrono::duration<float> difference = endTime - startTime;
		//	double milliseconds = difference.count();
		//	std::cout << "total \t" << total << std::endl;
		//	std::cout << "read string time (s) \t" << milliseconds << std::endl;
		//}
	}
}

void convertStrMultithread(const std::vector<size_t>& cut, const std::string& info,
	std::vector<int>& solId, std::vector<std::vector<int>>& sols) {
	int num_task = std::min<int>(std::thread::hardware_concurrency(), cut.size() - 1);
	//num_task = 1;
	std::vector<std::thread> thrds;
	std::pair<int, int> from_to;
	for (size_t i = 0; i < num_task; ++i) {
		thrds.push_back(std::thread(
			&convertStr, std::cref(info), cut[i] + 1, cut[i + 1],
			std::ref(solId), std::ref(sols)));
	}
	//thrds.back().join();
	for (auto& thrd : thrds)
		thrd.join();
}



void name_tsp::readSol(std::stringstream& buffer, 
	std::vector<int>& solId,
	std::vector<std::vector<int>>& sols) {
	std::vector<int> sol;
	std::string info = buffer.str();
	buffer.flush();
	std::string line;
	size_t from = 0;
	size_t to = info.find_first_of('\n', from);
	line = info.substr(from, to);
	from = to + 1;
	//total.push_back(header1.substr(from + 1, to));

	int total(0);
	int bufferSize = std::thread::hardware_concurrency() * 10;

	auto startTime = std::chrono::system_clock::now();


	std::vector<size_t> cut;
	cut.push_back(to);
	int avg_cut = (info.size() - to) / std::thread::hardware_concurrency();
	from = avg_cut + to;
	to = info.find_first_of('\n', from);

	while (to <= info.size()) {
		cut.push_back(to);
		from = to + avg_cut;
		if (from >= info.size()) {
			cut.push_back(info.size());
			break;
		}
		to = info.find_first_of('\n', from);
		if (to > info.size()) {
			cut.push_back(info.size());
			break;
		}
	}

	//std::vector<char> cutVal;
	//for (auto& it : cut) {
	//	cutVal.push_back(info[it]);
	//}
	//int lineSize = info.find_first_of('\n',0);

	size_t lastIdx = info.find_last_of('\n');
	size_t lastNotIdx = lastIdx - 2;
	while (info[lastNotIdx] != '\n') {
		--lastNotIdx;
	}
	//std::cout << "lastnotIdx\t" << lastNotIdx << "\t" << "lastIdx\t" << lastIdx << std::endl;
	std::string curInfo = info.substr(lastNotIdx + 1, lastNotIdx - lastIdx);
	std::stringstream curStr;
	curStr << curInfo;
	int lastSolId(0);
	curStr >> lastSolId;

	std::cout << "total number of solutions\t" << lastSolId << std::endl;

	sols.resize(lastSolId);
	solId.resize(lastSolId);


	convertStrMultithread(cut, info,
		solId, sols);

	//convertStr(info, cut.front()+1, cut.back(),
	//	solId, sols);

	auto endTime = std::chrono::system_clock::now();
	auto difference = endTime - startTime;
	auto milliseconds = difference.count();
	std::cout << "read info time (s) \t" << milliseconds << std::endl;


	//std::vector<std::string> totallines;
	//while (from < info.size()) {
	//	to = info.find_first_of('\n', from);
	//	line = info.substr(from, to);
	//	if (info[to] != '\n') {
	//		totallines.push_back(line);
	//		//		totallines.push_back(line);
	//	}
	//	totallines.push_back(line);
	//	from = to + 1;
	//	if (totallines.size() > bufferSize) {
	//		auto endTime = std::chrono::system_clock::now();
	//		std::chrono::duration<float> difference = endTime - startTime;
	//		double milliseconds = difference.count();
	//		std::cout << "read string time (s) \t" << milliseconds << std::endl;
	//		int from = solId.size();
	//		int to = from + totallines.size();
	//		solId.resize(to);
	//		sols.resize(to);
	//		convertMultithread(from, to, totallines,
	//			solId, sols);
	//		totallines.clear();
	//		endTime = std::chrono::system_clock::now();
	//		difference = endTime - startTime;
	//		milliseconds = difference.count();
	//		std::cout << "cur size\t" << sols.size() << std::endl;
	//		std::cout << "read info time (s) \t" << milliseconds << std::endl;
	//	}
	//}
	//solId.resize(totallines.size());
	//sols.resize(totallines.size());
}

void initSet(std::vector<int>& belong) {
	for (int idx(0); idx < belong.size(); ++idx) {
		belong[idx] = idx;
	}
}

int unionSet(std::vector<int>& belong, int id) {
	if (belong[id] == id)return id;
	else {
		return belong[id] = unionSet(belong, belong[id]);
	}
}
int mergeSet(std::vector<int>& belong, int a, int b) {
	if (unionSet(belong, a) != unionSet(belong, b)) {
		return belong[belong[b]] = belong[a];
	}
	return belong[a];
}


void lon::filterNetworkRemoveLoopWishSink(
	std::set<int> sinkNodes,
	const std::vector<TraceEdge>& edges,
	const std::vector<double>& nodefit,
	LonInfo& lon) {
	struct NetworkEdgeRSL {
		int id = 0;
		int m_from = 0;
		int m_to = 0;
		bool m_effective = false;
	};
	std::vector<int> sortedId(nodefit.size() - 1);
	for (int idx(1); idx < nodefit.size(); ++idx) {
		sortedId[idx - 1] = idx;
	}
	std::sort(sortedId.begin(), sortedId.end(), [&](int a, int b) {
		return nodefit[a] < nodefit[b];
	});

	std::vector<int> numParensts(nodefit.size(), 0);
	std::vector<int> numSons(nodefit.size(), 0);

	//std::vector<int> belongFitIds(nodefit.size(),-1);
	std::vector<int> belongFitVals(nodefit.size(), 0);


	std::vector<bool> isSinkNode(nodefit.size(), false);
	for (auto& it : sinkNodes)isSinkNode[it] = true;

	std::vector<std::vector<NetworkEdgeRSL*>> sons(nodefit.size());
	std::vector<NetworkEdgeRSL> networkEdges(edges.size());
	for (int idx(0); idx < edges.size(); ++idx) {
		if (edges[idx].id_to != -1 && edges[idx].id_to != edges[idx].id_from) {
			networkEdges[idx].m_from = edges[idx].id_from;
			networkEdges[idx].m_to = edges[idx].id_to;
			networkEdges[idx].id = idx;
			networkEdges[idx].m_effective = true;


			if (nodefit[edges[idx].id_to] < nodefit[edges[idx].id_from]) {
				sons[networkEdges[idx].m_to].push_back(&networkEdges[idx]);
				numParensts[networkEdges[idx].m_from]++;
				numSons[networkEdges[idx].m_to]++;
			}

			else if (nodefit[edges[idx].id_to] == nodefit[edges[idx].id_from]) {
				if (!isSinkNode[edges[idx].id_from]) {
					sons[networkEdges[idx].m_to].push_back(&networkEdges[idx]);
					numParensts[networkEdges[idx].m_from]++;
					numSons[networkEdges[idx].m_to]++;
				}
			}
			//else {
			//	int stop = -1;
			//}
		//	else networkEdges[idx].m_effective = true;

		}
		else {
			networkEdges[idx].m_from = edges[idx].id_from;
			networkEdges[idx].m_to = edges[idx].id_to;
			networkEdges[idx].id = idx;
			networkEdges[idx].m_effective = false;
			//		sons[networkEdges[idx].m_to].push_back(&networkEdges[idx]);
		}

	}

	//for (auto& edge : networkEdges) {
	//	if (edge.m_effective) {
	//	
	//	}
	//}


	std::vector<int> visited(nodefit.size(), false);
	std::queue<int> que_frontNode;

	// for test
	std::set<int> fitValues;

	for (auto& idx : sinkNodes) {

		que_frontNode.push(idx);
		belongFitVals[idx] = nodefit[idx];
		fitValues.insert(nodefit[idx]);
	}
	//for (int idx(1); idx < nodefit.size(); ++idx) {
	//	if (numParensts[idx] == 0 && numSons[idx] != 0&&sinkNodes.find(idx)!=sinkNodes.end()) {
	//		que_frontNode.push(idx);
	//		belongFitVals[idx] = nodefit[idx];
	//		fitValues.insert(nodefit[idx]);
	//	}
	//}

	std::vector<int> fitValuesVec;
	for (auto& it : fitValues) {
		fitValuesVec.push_back(it);
	}


	std::sort(fitValuesVec.begin(), fitValuesVec.end(), [](int a, int b) {
		return a < b;
	});

	std::map<int, int> fitValToFitIdx;

	for (int idx(0); idx < fitValuesVec.size(); ++idx) {
		fitValToFitIdx[fitValuesVec[idx]] = idx;
	}

	int firstIdx(0);
	while (true) {
		while (!que_frontNode.empty()) {
			int curId = que_frontNode.front();
			que_frontNode.pop();
			visited[curId] = true;
			for (auto& edge : sons[curId]) {
				if (visited[edge->m_from]) {
					edge->m_effective = false;
				}
				else {
					int nextId = edge->m_from;
					if (belongFitVals[nextId] == 0) {
						belongFitVals[nextId] = belongFitVals[curId];
					}
					else if (belongFitVals[nextId] != belongFitVals[curId]) {
						belongFitVals[nextId] = -1;
					}
					if (--numParensts[nextId] == 0) {
						que_frontNode.push(nextId);
					}
				}
			}
		}

		for (; firstIdx < sortedId.size(); ++firstIdx) {
			int idx = sortedId[firstIdx];
			if (!visited[idx]) {
				que_frontNode.push(idx);
				if (belongFitVals[idx] == 0) {
					belongFitVals[idx] = nodefit[idx];
				}
				//else if (belongFitVals[idx] != nodefit[]) {
				//	belongFitVals[nextId] = -1;
				//}
				break;
			}
		}
		if (firstIdx >= sortedId.size())break;
	}

	for (auto& curEdge : networkEdges) {
		if (curEdge.m_effective) {
			if (belongFitVals[curEdge.m_from] != belongFitVals[curEdge.m_to]) {
				curEdge.m_effective = false;
			}
		}
	}


	for (int idx(0); idx < edges.size(); ++idx) {
		if (edges[idx].id_to != -1 && edges[idx].id_to != edges[idx].id_from) {
			networkEdges[idx].m_from = edges[idx].id_from;
			networkEdges[idx].m_to = edges[idx].id_to;
			networkEdges[idx].id = idx;
			networkEdges[idx].m_effective = true;
			if (nodefit[edges[idx].id_to] < nodefit[edges[idx].id_from]) {
				sons[networkEdges[idx].m_to].push_back(&networkEdges[idx]);
			}

		}
	}


	lon.m_node_fit.resize(nodefit.size() - 1);
	lon.m_node_funnel_idxs.resize(nodefit.size() - 1);
	lon.m_node_funnel_fit.resize(nodefit.size() - 1);
	lon.m_nodeId.resize(nodefit.size() - 1);
	//std::set<int> colorFunnels;
	for (int idx(1); idx < nodefit.size(); ++idx) {
		lon.m_node_fit[idx - 1] = nodefit[idx];
		if (belongFitVals[idx] < 0) {
			lon.m_node_funnel_idxs[idx - 1] = -1;
		}
		else
			lon.m_node_funnel_idxs[idx - 1] = fitValToFitIdx[belongFitVals[idx]];
		lon.m_node_funnel_fit[idx - 1] = belongFitVals[idx];
		//colorFunnels.insert(belongFitVals[idx]);
		lon.m_nodeId[idx - 1] = idx;
	}

	lon.m_node_size.resize(nodefit.size() - 1, 0);
	std::fill(lon.m_node_size.begin(), lon.m_node_size.end(), 0);

	std::vector<int> curEdge(2, 0);
	for (auto& edge : networkEdges) {
		if (edge.m_effective) {
			curEdge.front() = edge.m_from;
			curEdge.back() = edge.m_to;
			lon.m_node_size[edge.m_to - 1]++;
			lon.m_graph.push_back(curEdge);
		}
	}
}

void lon::filterNetworkRemoveLoopWishSink(
	std::set<int> sinkNodes,
	const std::vector<TraceEdge>& edges,
	const std::vector<int>& nodefit,
	LonInfo& lon) {
	struct NetworkEdgeRSL {
		int id = 0;
		int m_from = 0;
		int m_to = 0;
		bool m_effective = false;
	};
	std::vector<int> sortedId(nodefit.size() - 1);
	for (int idx(1); idx < nodefit.size(); ++idx) {
		sortedId[idx - 1] = idx;
	}
	std::sort(sortedId.begin(), sortedId.end(), [&](int a, int b) {
		return nodefit[a] < nodefit[b];
		});

	std::vector<int> numParensts(nodefit.size(), 0);
	std::vector<int> numSons(nodefit.size(), 0);

	//std::vector<int> belongFitIds(nodefit.size(),-1);
	std::vector<int> belongFitVals(nodefit.size(), 0);


	std::vector<bool> isSinkNode(nodefit.size(), false);
	for (auto& it : sinkNodes)isSinkNode[it] = true;

	std::vector<std::vector<NetworkEdgeRSL*>> sons(nodefit.size());
	std::vector<NetworkEdgeRSL> networkEdges(edges.size());
	for (int idx(0); idx < edges.size(); ++idx) {
		if (edges[idx].id_to != -1 && edges[idx].id_to != edges[idx].id_from) {
			networkEdges[idx].m_from = edges[idx].id_from;
			networkEdges[idx].m_to = edges[idx].id_to;
			networkEdges[idx].id = idx;
			networkEdges[idx].m_effective = true;


			if (nodefit[edges[idx].id_to] < nodefit[edges[idx].id_from]) {
				sons[networkEdges[idx].m_to].push_back(&networkEdges[idx]);
				numParensts[networkEdges[idx].m_from]++;
				numSons[networkEdges[idx].m_to]++;
			}

			else if (nodefit[edges[idx].id_to] == nodefit[edges[idx].id_from]) {
				if (!isSinkNode[edges[idx].id_from]) {
					sons[networkEdges[idx].m_to].push_back(&networkEdges[idx]);
					numParensts[networkEdges[idx].m_from]++;
					numSons[networkEdges[idx].m_to]++;
				}
			}
			//else {
			//	int stop = -1;
			//}
		//	else networkEdges[idx].m_effective = true;

		}
		else {
			networkEdges[idx].m_from = edges[idx].id_from;
			networkEdges[idx].m_to = edges[idx].id_to;
			networkEdges[idx].id = idx;
			networkEdges[idx].m_effective = false;
			//		sons[networkEdges[idx].m_to].push_back(&networkEdges[idx]);
		}

	}

	//for (auto& edge : networkEdges) {
	//	if (edge.m_effective) {
	//	
	//	}
	//}


	std::vector<int> visited(nodefit.size(), false);
	std::queue<int> que_frontNode;

	// for test
	std::set<int> fitValues;

	for (auto& idx : sinkNodes) {

		que_frontNode.push(idx);
		belongFitVals[idx] = nodefit[idx];
		fitValues.insert(nodefit[idx]);
	}
	//for (int idx(1); idx < nodefit.size(); ++idx) {
	//	if (numParensts[idx] == 0 && numSons[idx] != 0&&sinkNodes.find(idx)!=sinkNodes.end()) {
	//		que_frontNode.push(idx);
	//		belongFitVals[idx] = nodefit[idx];
	//		fitValues.insert(nodefit[idx]);
	//	}
	//}

	std::vector<int> fitValuesVec;
	for (auto& it : fitValues) {
		fitValuesVec.push_back(it);
	}


	std::sort(fitValuesVec.begin(), fitValuesVec.end(), [](int a, int b) {
		return a < b;
		});

	std::map<int, int> fitValToFitIdx;

	for (int idx(0); idx < fitValuesVec.size(); ++idx) {
		fitValToFitIdx[fitValuesVec[idx]] = idx;
	}

	int firstIdx(0);
	while (true) {
		while (!que_frontNode.empty()) {
			int curId = que_frontNode.front();
			que_frontNode.pop();
			visited[curId] = true;
			for (auto& edge : sons[curId]) {
				if (visited[edge->m_from]) {
					edge->m_effective = false;
				}
				else {
					int nextId = edge->m_from;
					if (belongFitVals[nextId] == 0) {
						belongFitVals[nextId] = belongFitVals[curId];
					}
					else if (belongFitVals[nextId] != belongFitVals[curId]) {
						belongFitVals[nextId] = -1;
					}
					if (--numParensts[nextId] == 0) {
						que_frontNode.push(nextId);
					}
				}
			}
		}

		for (; firstIdx < sortedId.size(); ++firstIdx) {
			int idx = sortedId[firstIdx];
			if (!visited[idx]) {
				que_frontNode.push(idx);
				if (belongFitVals[idx] == 0) {
					belongFitVals[idx] = nodefit[idx];
				}
				//else if (belongFitVals[idx] != nodefit[]) {
				//	belongFitVals[nextId] = -1;
				//}
				break;
			}
		}
		if (firstIdx >= sortedId.size())break;
	}

	for (auto& curEdge : networkEdges) {
		if (curEdge.m_effective) {
			if (belongFitVals[curEdge.m_from] != belongFitVals[curEdge.m_to]) {
				curEdge.m_effective = false;
			}
		}
	}


	for (int idx(0); idx < edges.size(); ++idx) {
		if (edges[idx].id_to != -1 && edges[idx].id_to != edges[idx].id_from) {
			networkEdges[idx].m_from = edges[idx].id_from;
			networkEdges[idx].m_to = edges[idx].id_to;
			networkEdges[idx].id = idx;
			networkEdges[idx].m_effective = true;
			if (nodefit[edges[idx].id_to] < nodefit[edges[idx].id_from]) {
				sons[networkEdges[idx].m_to].push_back(&networkEdges[idx]);
			}

		}
	}


	lon.m_node_fit.resize(nodefit.size() - 1);
	lon.m_node_funnel_idxs.resize(nodefit.size() - 1);
	lon.m_node_funnel_fit.resize(nodefit.size() - 1);
	lon.m_nodeId.resize(nodefit.size() - 1);
	//std::set<int> colorFunnels;
	for (int idx(1); idx < nodefit.size(); ++idx) {
		lon.m_node_fit[idx - 1] = nodefit[idx];
		if (belongFitVals[idx] < 0) {
			lon.m_node_funnel_idxs[idx - 1] = -1;
		}
		else
			lon.m_node_funnel_idxs[idx - 1] = fitValToFitIdx[belongFitVals[idx]];
		lon.m_node_funnel_fit[idx - 1] = belongFitVals[idx];
		//colorFunnels.insert(belongFitVals[idx]);
		lon.m_nodeId[idx - 1] = idx;
	}

	lon.m_node_size.resize(nodefit.size() - 1, 0);
	std::fill(lon.m_node_size.begin(), lon.m_node_size.end(), 0);

	std::vector<int> curEdge(2, 0);
	for (auto& edge : networkEdges) {
		if (edge.m_effective) {
			curEdge.front() = edge.m_from;
			curEdge.back() = edge.m_to;
			lon.m_node_size[edge.m_to - 1]++;
			lon.m_graph.push_back(curEdge);
		}
	}
}

void lon::getSinkNode3(std::set<int> funnelSetIds, const std::vector<TraceEdge>& edges,
	const std::vector<double>& nodefit,
	std::set<int>& sinkNodes) {
	std::vector<int> belong(nodefit.size());
	initSet(belong);
	for (auto& edge : edges) {
		if (edge.id_to != -1) {
			if (nodefit[edge.id_from] == nodefit[edge.id_to]) {
				mergeSet(belong, edge.id_from, edge.id_to);
			}
		}
	}
	std::vector<int> numOutEdge(nodefit.size(), 0);
	for (auto& edge : edges) {
		if (edge.id_to != -1) {
			if (nodefit[edge.id_from] > nodefit[edge.id_to]) {
				int parent = unionSet(belong, edge.id_from);
				++numOutEdge[parent];
			}
		}
	}

	for (auto& it : funnelSetIds) {
		int parent = unionSet(belong, it);
		if (numOutEdge[parent] <= 1) {
			sinkNodes.insert(it);
		//	std::cout << "nodefit\t" << nodefit[it] << "node Id\t" << it << std::endl;
		}
	}

}



void lon::getSinkNode3(std::set<int> funnelSetIds, const std::vector<TraceEdge>& edges,
	const std::vector<int>& nodefit,
	std::set<int>& sinkNodes) {
	std::vector<int> belong(nodefit.size());
	initSet(belong);
	for (auto& edge : edges) {
		if (edge.id_to != -1) {
			if (nodefit[edge.id_from] == nodefit[edge.id_to]) {
				mergeSet(belong, edge.id_from, edge.id_to);
			}
		}
	}
	std::vector<int> numOutEdge(nodefit.size(), 0);
	for (auto& edge : edges) {
		if (edge.id_to != -1) {
			if (nodefit[edge.id_from] > nodefit[edge.id_to]) {
				int parent = unionSet(belong, edge.id_from);
				++numOutEdge[parent];
			}
		}
	}

	for (auto& it : funnelSetIds) {
		int parent = unionSet(belong, it);
		if (numOutEdge[parent] <= 1) {
			sinkNodes.insert(it);
			//	std::cout << "nodefit\t" << nodefit[it] << "node Id\t" << it << std::endl;
		}
	}

}


void lon::getFunnelIdx(std::set<int>& funnelSetIds,
	std::vector<TraceEdge>& trace//, 
	//const std::vector<TraceEdge>& edges,
	//const std::vector<int>& nodefit
) {
	//std::set<int> funnelSetIds;
	std::vector<int> totalRun(1e6, -1);
	std::vector<TraceEdge> finalTrace(1e6);
	for (auto& it : trace) {
		if (it.id_to == -1) {
			totalRun[it.run] = it.id_from;
			finalTrace[it.run] = it;
		}
	}
	for (auto& it : totalRun) {
		if (it != -1) {
			funnelSetIds.insert(it);
		}
	}
}


void lon::filterLON(LonInfo& lon, LonInfo& lonfileter, double ratio, std::vector<int>& originId2newId) {
	originId2newId.resize(lon.m_node_fit.size());
	std::fill(originId2newId.begin(), originId2newId.end(), -1);
	int minValue(std::numeric_limits<int>::max());

	for (auto& it : lon.m_node_fit) minValue = std::min(minValue, it);
	int minValueCut = minValue * (1.0 + 0.0001);



	int newIdNum(0);
	for (int idx(0); idx < lon.m_node_fit.size(); ++idx) {

		if ((lon.m_node_fit[idx] - minValue) <= ratio * minValue) {
			originId2newId[idx] = newIdNum++;
		}
	}

	std::vector<int> newId2originId(newIdNum);
	for (int idx(0); idx < originId2newId.size(); ++idx) {
		if (originId2newId[idx] >= 0) {
			newId2originId[originId2newId[idx]] = idx;
		}
	}

	lonfileter.m_node_fit.resize(newIdNum);
	lonfileter.m_node_funnel_idxs.resize(newIdNum);
	lonfileter.m_node_funnel_fit.resize(newIdNum);
	lonfileter.m_node_size.resize(newIdNum);
	lonfileter.m_nodeId.resize(newIdNum);
	for (int idx(0); idx < newIdNum; ++idx) {
		int originId = newId2originId[idx];
		lonfileter.m_node_fit[idx] = lon.m_node_fit[originId];
		lonfileter.m_node_funnel_idxs[idx] = lon.m_node_funnel_idxs[originId];
		lonfileter.m_node_funnel_fit[idx] = lon.m_node_funnel_fit[originId];
		lonfileter.m_node_size[idx] = lon.m_node_size[originId];
		lonfileter.m_nodeId[idx] = lon.m_nodeId[originId];

	}

	lonfileter.m_graph.clear();
	std::vector<int> curEdge(2);

	for (auto& it : lon.m_graph) {
		curEdge.front() = originId2newId[it.front() - 1] + 1;
		curEdge.back() = originId2newId[it.back() - 1] + 1;
		if (curEdge.front() >= 1 && curEdge.back() >= 1) {
			lonfileter.m_graph.push_back(curEdge);
		}
		//else if (curEdge.front() >= 1) {
		//	int stop = -1;
		//}

	}

	//int stop = -1;
}

void lon::readTrace(const std::string& filepath, std::vector<TraceEdge>& trace) {
	std::string line;
	std::ifstream infile(filepath);;
	std::getline(infile, line);
	trace.clear();
	TraceEdge curInfo;
	while (infile >> curInfo.run >> curInfo.iter >> curInfo.id_from >> curInfo.id_to >> curInfo.num_kick) {
		trace.push_back(curInfo);
	}
	infile.close();
}

void lon::readNodeInfo(const std::string& filepath, std::vector<double>& nodefit) {
	int id(0), fit(0);
	std::ifstream infile(filepath);
	std::string line;
	std::getline(infile, line);
	while (infile >> id >> fit) {
		if (nodefit.size() <= id) {
			nodefit.resize(id + 1);
		}
		nodefit[id] = fit;
	}
	infile.close();
}
void lon::readEdge(const std::string& filepath, std::vector<TraceEdge>& edges) {
	TraceEdge curEdge;

	std::string line;
	std::ifstream infile(filepath);;
	std::getline(infile, line);
	edges.clear();
	while (infile >> curEdge.id_from >> curEdge.id_to >> curEdge.count) {
		edges.push_back(curEdge);
	}
	infile.close();
}


void lon::transfer(const std::vector<TraceEdge>& trace, std::vector<TraceEdge>& edges) {
	std::map<std::pair<int, int>, int> edgeInfo2EdgeId;
	std::vector<int> edgeVisited;
	int maxId = 0;
	int curId = 0;
	std::pair<int, int> edgeInfo;
	std::vector<std::pair<int, int>> edgeId2EdgeInfo;
	for (auto& it : trace) {
		edgeInfo.first = it.id_from;
		edgeInfo.second = it.id_to;

		if (edgeInfo2EdgeId.find(edgeInfo) == edgeInfo2EdgeId.end()) {
			curId = edgeInfo2EdgeId[edgeInfo] = maxId++;
			edgeVisited.push_back(0);
			edgeId2EdgeInfo.push_back(edgeInfo);
		}
		else {
			curId = edgeInfo2EdgeId[edgeInfo];
		}
		++edgeVisited[curId];
		//int curId= edgeInfo2EdgeId[]
	}

	edges.resize(edgeVisited.size());
	for (int idx(0); idx < edgeInfo2EdgeId.size(); ++idx) {
		edges[idx].id_from = edgeId2EdgeInfo[idx].first;
		edges[idx].id_to = edgeId2EdgeInfo[idx].second;
		edges[idx].count = edgeVisited[idx];

		//edgeOut << edgeId2EdgeInfo[idx].first << "\t" << edgeId2EdgeInfo[idx].second << "\t";
		//edgeOut << edgeVisited[idx] << std::endl;
	}
}
//

void lon::outputTrace(const std::string& filepath, std::vector<TraceEdge>& trace) {

	std::ofstream edgeHistoryOut(filepath + ".edge_history");
	edgeHistoryOut << "RUN\tITER\tID_START\tID_END\tNUM_KICKS" << std::endl;

	for (auto& it : trace) {
		edgeHistoryOut << it.run << "\t" << it.iter << "\t" << it.id_from << "\t" << it.id_to << "\t" << it.num_kick << std::endl;
	}
	edgeHistoryOut.close();
}
void lon::outputNodeInfo(const std::string& filepath, std::vector<int>& nodefit) {
	std::ofstream nodeFitOut(filepath + ".nodes");
	nodeFitOut << "ID\tFITNESS" << std::endl;
	for (int idx(1); idx < nodefit.size(); ++idx) {
		nodeFitOut << idx << "\t" << nodefit[idx] << std::endl;
	}
	nodeFitOut.close();
}
void lon::outputEdge(const std::string& filepath, std::vector<TraceEdge>& edges) {
	std::ofstream edgeOut(filepath+ ".edges");
	edgeOut << "ID_START\tID_END\tCOUNT" << std::endl;
	for (int idx(0); idx < edges.size(); ++idx) {
		edgeOut << edges[idx].id_from << "\t" << edges[idx].id_to<< "\t";
		edgeOut << edges[idx].count<< std::endl;

	}
	edgeOut.close();
}


void lon::outputNodeSol(const std::string& filepath, std::vector<std::vector<int>>& nodeSol) {
	//int solId(0);

	std::ofstream nodeSolOut(filepath + ".sols");
	nodeSolOut << "ID\tSOLUTION" << std::endl;
	for (int solId(1); solId < nodeSol.size(); ++solId) {
		nodeSolOut << solId;
		auto& sol = nodeSol[solId];
		for (int idx(0); idx < sol.size(); ++idx) {
			nodeSolOut << "\t" << sol[idx];
		}
		nodeSolOut << std::endl;
	}
	nodeSolOut << std::endl;

	
}

void lon::outputLonNodeInfo(const std::string& filepath,const LonInfo& lon) {
	std::ofstream lonOut(filepath + ".lonNodeInfo");
	lonOut << "ID\tFunnelIdx\tFunnelFit\tnodeSize\tnodeFit\t" << std::endl;
	for (int idx(0); idx < lon.m_nodeId.size(); ++idx) {
		lonOut << lon.m_nodeId[idx] << "\t";
		lonOut << lon.m_node_funnel_idxs[idx] << "\t";
		lonOut << lon.m_node_funnel_fit[idx] << "\t";
		lonOut << lon.m_node_size[idx] << "\t";
		lonOut << lon.m_node_fit[idx] << "\t";
		lonOut << std::endl;
	}
	lonOut.close();
}



void lon::outputLonEdgeInfo(const std::string& filepath, const LonInfo& lon) {
	std::ofstream lonOut(filepath + ".lonEdgeInfo");
	lonOut << "EdgeFrom\tEdgeTo" << std::endl;
	for (int idx(0); idx < lon.m_graph.size(); ++idx) {
		lonOut << lon.m_graph[idx].front() << "\t" << lon.m_graph[idx].back();
		lonOut << std::endl;
	}
	lonOut.close();
}


void lon::insertLonNodeInfo(const std::string& filepath, LonInfo& lon) {
	std::ifstream lonIn(filepath + ".lonNodeInfo");
	std::string line;
	std::getline(lonIn, line);
	int lonId(0), funelIdx(0), funnelFit(0), nodeSize(0), nodeFit(0);
	lon.m_nodeId.clear();
	lon.m_node_funnel_idxs.clear();
	lon.m_node_funnel_fit.clear();
	lon.m_node_fit.clear();
	lon.m_node_size.clear();
	while (lonIn >> lonId >> funelIdx >> funnelFit >> nodeSize >> nodeFit) {
		lon.m_nodeId.push_back(lonId);
		lon.m_node_funnel_idxs.push_back(funelIdx);
		lon.m_node_funnel_fit.push_back(funnelFit);
		lon.m_node_size.push_back(nodeSize);
		lon.m_node_fit.push_back(nodeFit);
	}
	lonIn.close();
}
void lon::insertLonEdgeInfo(const std::string& filepath, LonInfo& lon) {
	std::ifstream lonIn(filepath + ".lonEdgeInfo");
	std::string line;
	std::getline(lonIn, line);
	std::vector<int> edge(2);
	lon.m_graph.clear();
	while (lonIn >> edge.front()>>edge.back()) {
		lon.m_graph.push_back(edge);
	}

	lonIn.close();
}