#include "nbn_fla_tsp.h"
#include "../../../../instance/problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../../../core/parameter/variants_to_stream.h"

namespace ofec {
	namespace nbn {
		namespace tsp {

			void getNarrowGap(const NBNinfo& nbnInfo, std::vector<int>& bridgeIds, double bestFitThreadhold, double basinSisethread) {
				std::vector<std::vector<int>> sons;
				getDirectBasin(nbnInfo, sons);
				std::vector<int> basinSize;
				getBasinSize(nbnInfo, sons, basinSize);



				auto& fitness = nbnInfo.m_vFitness;
				auto& belong = nbnInfo.m_belong;
				auto& dis2parent = nbnInfo.m_dis2parent;


				auto norFit = fitness;
				ofec::dataNormalize(norFit);

				int numBasinSize = nbnInfo.m_belong.size() * basinSisethread;

				for (int idx(0); idx < nbnInfo.m_vFitness.size(); ++idx) {
					if (basinSize[idx] >= numBasinSize && norFit[idx] >= bestFitThreadhold) {
						bridgeIds.push_back(idx);
					}
				}

				std::vector<bool> visited(fitness.size(), false);
				std::queue<int> queIds;
				for (auto& it : bridgeIds) {
					visited[it] = true;
					queIds.push(it);
				}
				while (!queIds.empty()) {
					int cur = queIds.front();
					queIds.pop();
					if (belong[cur] != -1 && !visited[belong[cur]]) {
						queIds.push(belong[cur]);
						visited[belong[cur]] = true;
					}
				}


			}

			void calculateNarrowGapInfo(const NBNinfo& nbnInfo, const std::vector<int>& bridgeIds, GapInfo& info) {

				std::vector<std::vector<int>> sons;
				getDirectBasin(nbnInfo, sons);
				std::vector<int> basinSize;
				getBasinSize(nbnInfo, sons, basinSize);




				std::vector<bool> visited(nbnInfo.m_vFitness.size(), false);
				for (auto& it : bridgeIds) visited[it] = true;
				std::vector<int> leafBridgeIds(nbnInfo.m_vFitness.size(), 0);
				for (int idx(0); idx < bridgeIds.size(); ++idx) {
					for (auto& sonIter : sons[bridgeIds[idx]]) {
						if (visited[sonIter]) {
							++leafBridgeIds[bridgeIds[idx]];
						}
					}
				}

				
				auto bridgeIdsCopy = bridgeIds;


				std::sort(bridgeIdsCopy.begin(), bridgeIdsCopy.end(), [&](int a,int b) {
					if (nbnInfo.m_vFitness[a] == nbnInfo.m_vFitness[b]) {
						return leafBridgeIds[a] < leafBridgeIds[b];
					}
					else
					return nbnInfo.m_vFitness[a] < nbnInfo.m_vFitness[b];
					});




			

				
				struct Info {
					std::vector<int> m_path;
					double m_maxL = 0;
					
				};

				std::vector<Info> paths(bridgeIdsCopy.size());
	
		/*		for (auto& it : bridgeIds) {
					visited[it] = true;
				}*/
				std::vector<int> bridgeOriId(nbnInfo.m_vFitness.size(), -1);
				for (int idx(0); idx < bridgeIdsCopy.size(); ++idx) {
					bridgeOriId[bridgeIdsCopy[idx]] = idx;
				}
			
				for (int idx(0); idx < bridgeIdsCopy.size(); ++idx) {
					auto& curPath = paths[idx];
					auto& cur = bridgeIdsCopy[idx];
					if (nbnInfo.m_dis2parent[cur] > 1e9) {
						curPath.m_maxL += 0;
					}
					else 		curPath.m_maxL += nbnInfo.m_dis2parent[cur];
					curPath.m_path.push_back(cur);
					if (nbnInfo.m_belong[cur] != -1) {
						auto parentId = nbnInfo.m_belong[cur];

						int parentIdx = bridgeOriId[parentId];
						if (parentIdx >=0) {

							auto& parentPath = paths[parentIdx];
							if (curPath.m_maxL > parentPath.m_maxL) {
								parentPath = curPath;
							}
						}
					}
				}

				std::vector<GapInfo> maxPaths(paths.size());
				for (int idx(0); idx < maxPaths.size(); ++idx) {
					auto& info = maxPaths[idx];
					auto& maxPath = paths[idx];
					info.m_maxlength = maxPath.m_maxL;
					info.m_numNodes = maxPath.m_path.size();
					double sumNode(0);

					std::vector<int> nodeBasinSize;

					for (auto& it : maxPath.m_path) {
						sumNode += basinSize[it];
						nodeBasinSize.push_back(basinSize[it]);
					}
					info.m_avgBasinSize = sumNode / info.m_numNodes;
					ofec::calMin(nodeBasinSize, info.m_avgBasinSize);
				}


				for (auto& it : maxPaths) {
					if (info.m_maxlength < it.m_maxlength) {
						info = it;
					}
					else if (info.m_maxlength == it.m_maxlength) {
						if (info.m_avgBasinSize > it.m_avgBasinSize) {
							info = it;
						}
					}
				}



				// calculate direct sons
				{
					info.m_minDirectionSons = std::numeric_limits<int>::max();
					for (auto& it : bridgeIds) {
						info.m_minDirectionSons = std::min<int>(info.m_minDirectionSons, sons[it].size());
					}
				}

				// calculate max NBD

				{
					info.maxNBD = 0;
					for (auto& it : bridgeIds) {
						if(nbnInfo.m_dis2parent[it]<1)
						info.maxNBD = std::max<int>(info.maxNBD, nbnInfo.m_dis2parent[it]);
					}
				}


			}



			void TspSolutionsToNBNinfo(NBNinfo& nbnInfo, const std::string& saveDir, const std::string& filename, ofec::Environment* env) {

				nbnInfo.clear();
				auto& m_sols = nbnInfo.m_sols;
				ofec::ParameterVariantStream paramsStream;
				variants_stream::inputFromFile(paramsStream, saveDir + filename + "_solVariables.txt");
				size_t solSize(0);

				paramsStream >> solSize;
				m_sols.resize(solSize);
				auto pro = env->problem();

				for (auto& it : m_sols) {
					it.reset(pro->createSolution());
					auto& cursol = dynamic_cast<ofec::TravellingSalesman::SolutionType&>(*it);
					paramsStream >> cursol.variable().vect();
				}

				auto& solbases = nbnInfo.m_solbases;
				solbases.clear();
				for (auto& it : m_sols) {
					solbases.push_back(it.get());
				}

			}
		}
	}
}