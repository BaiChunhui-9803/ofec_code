
#ifndef OFEC_NBN_LEARNING2_H
#define OFEC_NBN_LEARNING2_H

#include "../../../../../core/problem/solution.h"
#include <memory>
#include <iomanip>
#include <set>
#include <map>
#include <vector>
#include <stack>
#include <queue>

#include "../../../../../utility/random/newran.h"

namespace ofec {
	namespace NBN_learn2 {


		class NetworkNode;
		class RelationNode;
		class RelationMapNode;
		class RelationListNode;

		class RelationShip;
		class RelationMap;
		class RelationList;



		//using RelationMap = std::map<unsigned long long, std::shared_ptr<RelationMapNode>>;
		class IndividualInfo {
		//	int m_lastUpdateTime = 0;
		//	int m_calculatedTimes = 0;
			double m_minDis = 0;

		public:			
	//		static unsigned long long  m_totalId;
	//		unsigned long long m_id = 0;



			IndividualInfo() = default;

			virtual double distance(const IndividualInfo& otherInfo)const {
				return 0;
			}
			virtual double fitness()const {
				return 0;
			}

			double minDis()const {
				return m_minDis;
			}


			virtual void setMinDis(double dis) {
				m_minDis = dis;
			}
			// -1: farer from, 0 equal to minDis, 1 smaller than minDis
			virtual int judgeNearerTo(const IndividualInfo& otherInfo,
				double& curdis) {
				curdis = distance(otherInfo);
				if (curdis > m_minDis) {
					return -1;
				}
				else if (curdis == m_minDis) {
					return 0;
				}
				else return 1;
			}

			virtual void initialize() {
			//	m_lastUpdateTime = 0;
			//	m_calculatedTimes = 0;
				m_minDis = std::numeric_limits<double>::max();
			}
		};



		

		class RelationNode : public EnableSharedPtr<RelationNode> {
			friend class RelationShip;
			friend class RelationList;
		protected:

			bool m_active = false;
			std::weak_ptr<RelationNode> m_parentSon;

			bool m_visited = false;
			
		public:


			RelationNode() = default;
			//~RelationNode() {
			//	clearMemory();
			//}
			virtual NetworkNode* getCurInfo() = 0;
			virtual RelationShip* getRelationShip() = 0;


			virtual std::shared_ptr<RelationNode>* getNext() {
				return nullptr;
			}
			bool isActive()const {
				return m_active;
			}

			virtual void setAfterNode(std::shared_ptr<RelationNode>& afterNode) {}
			virtual void setBeforeNode(RelationNode* beforeNode) {}


			virtual void initialize() {
				m_active = false;
				m_parentSon.reset();
			}
			virtual void removeFromStruct() = 0;
			virtual void removeNode() {

				removePSfromStruct();
				this->removeFromStruct();

			}
			virtual void removePSfromStruct() {
				{
					/*if (m_parentSon)*/ {
						if (!m_parentSon.expired()) {
							auto p = m_parentSon.lock();
							if (p) p->removeFromStruct();
						}
						m_parentSon.reset();
					}
				}
			}
			void clearMemory() {
				removeNode();
				m_active = false;
			}

			virtual void set(
				std::shared_ptr<NetworkNode>& cur,
				std::shared_ptr<RelationNode>& parentSon,
				RelationShip* curStruct
			) {
				m_active = true;
				m_parentSon = parentSon;
				m_visited = false;
			}
		};


		class RelationMapNode : public RelationNode{

			friend class RelationMap;
		protected:

		public:

			NetworkNode* m_cur = nullptr;
			RelationMap* m_pMap = nullptr;

			virtual NetworkNode* getCurInfo() override {
				return m_cur;
			}

			virtual RelationShip* getRelationShip()override;
			virtual void removeFromStruct() override;
			

			virtual void set(
				std::shared_ptr<NetworkNode>& cur,
				std::shared_ptr<RelationNode>& parentSon,
				RelationShip* curStruct
			)override;

		};
		

		class RelationListNode : public RelationNode{
			friend class RelationList;
		protected:

		public:
			std::shared_ptr<NetworkNode> m_cur = nullptr;
			RelationNode* m_before = nullptr;
			std::shared_ptr<RelationNode> m_after = nullptr;
			RelationList* m_list = nullptr;

			virtual void set(
				std::shared_ptr<NetworkNode>& cur,
				std::shared_ptr<RelationNode>& parentSon,
				RelationShip* curStruct
			)override;

			virtual NetworkNode* getCurInfo() override {
				return m_cur.get();
			}
			virtual RelationShip* getRelationShip() override;


			virtual std::shared_ptr<RelationNode>* getNext() override {
				return &m_after;
			}
			virtual void setAfterNode(std::shared_ptr<RelationNode>& afterNode) override;
			virtual void setBeforeNode(RelationNode* beforeNode) override;

			virtual void removeFromStruct() override;

		};


		class RelationShip {

		protected:
			int m_visited = 0;
		public:
			virtual void initialize() {
				m_visited = 0;
			}
			virtual void clear() = 0;
			virtual size_t size() const = 0;
		};

		class RelationMap : public RelationShip {
			friend class RelationMapNode;
			friend class NetworkNode;
		protected:


		public:

			std::map<NetworkNode*, std::shared_ptr<RelationNode>> m_map;

			virtual void remove(NetworkNode* cur) {
				m_map.erase(cur);
			}
			void insert(RelationMapNode* cur) {
				m_map[cur->m_cur] = cur->getSharedPtr();
				cur->m_pMap = this;
			}

			virtual void initialize() override {
				RelationShip::initialize();
				m_map.clear();
			}

			virtual void clear() {
				for (auto& it : m_map) {
				//	auto pit = it.second->getSharedPtr();
					it.second->removePSfromStruct();
				}
				m_map.clear();
			}

			bool judgeInside(NetworkNode* cur) const {
				return m_map.find(cur) != m_map.end();
			}

			virtual size_t size() const{
				return m_map.size();
			}


		};

		struct RelationList : public RelationShip {
		//	friend class RelationListNode;
		//	friend class NetworkNode;
		//	
		//protected:

			std::shared_ptr<RelationListNode> m_head;
			int m_size = 0;
			void insert(RelationListNode* cur);
		//public:

			virtual void initialize() override {

				RelationShip::initialize();
				m_size = 0;
				m_head.reset(new RelationListNode);
				std::shared_ptr<RelationListNode> end(new RelationListNode);
				m_head->m_after = end;
			//	m_head->m_list = this;
				end->m_before = m_head.get();
			//	end->m_list = this;
			}

			virtual void clear() {
				while (m_head->m_after->m_active) {
					auto p = m_head->m_after;
					p->removeNode();
				}
			}

			std::shared_ptr<RelationListNode>* getHead() {
				return &m_head;
			}


			virtual size_t size() const{
				return m_size;
			}
		};

		class NetworkNode : public EnableSharedPtr<NetworkNode> {

			std::shared_ptr<IndividualInfo> m_sol_info;
			RelationMap m_parents;
			RelationList m_sons;

			NetworkNode* m_direct_parent = nullptr;
			//RelationList m_sons;
			int m_basinSize = 1;

			int m_updateStamp = 0;
			bool m_lazy_update = false;

		public:


			~NetworkNode() {
				clear();
			}


			NetworkNode* directParent() {
				return m_direct_parent;
			}
			std::shared_ptr<IndividualInfo>& getSharedSolInfo() {
				return m_sol_info;
			}
			IndividualInfo* getSolInfo() {
				return m_sol_info.get();
			}
			RelationList& sonList() {
				return m_sons;
			}

			inline double minDis()const {
				return m_sol_info->minDis();
			}

			inline double fitness()const {
				return m_sol_info->fitness();
			}

			inline double distanceTo(NetworkNode* other) const {
				return m_sol_info->distance(*other->m_sol_info);
			}

			int calFitRel(NetworkNode* other) {
				if (m_sol_info->fitness() < other->fitness()) {
					return 1;
				}
				else if (m_sol_info->fitness() == other->fitness()) {
					return 0;
				}
				else return -1;
			}
			
			int calDisRel(NetworkNode* other, double& curdis) {
				return m_sol_info->judgeNearerTo(*other->m_sol_info,curdis);
			}

			

			void updateSons(std::shared_ptr<NetworkNode>& newNode, ofec::Random* rnd) {
				
				auto sonIter = m_sons.m_head->getSharedPtr();
				int reFlag(0);
				while ((*sonIter->getNext())->isActive()) {
					auto& cur = (*sonIter->getNext());
					double curdis(0);
					reFlag = cur->getCurInfo()->calFitRel(newNode.get());
					if (reFlag == 1) {
						reFlag = cur->getCurInfo()->calDisRel(newNode.get(), curdis);
					}
					if (reFlag == 1) {
						std::shared_ptr<RelationNode> cur_ptr = cur;
						cur_ptr->removeNode();
						auto curNetworkNode = cur_ptr->getCurInfo()->getSharedPtr();
						NetworkNode::bindRelationship(newNode,curNetworkNode);
						NetworkNode::bindDirectRelationship(newNode, curNetworkNode,curdis);
					}
					else if (reFlag == 0) {

						auto curNetworkNode = cur->getCurInfo()->getSharedPtr();
						NetworkNode::bindRelationship(newNode, curNetworkNode);
						if (cur->getRelationShip()->size() * rnd->uniform.next() < 1.0) {
							NetworkNode::bindDirectRelationship(newNode, curNetworkNode, curdis);
						}
						sonIter = (*sonIter->getNext());
					}
					else {
						sonIter = (*sonIter->getNext());
					}
				}
			}



			static void bindRelationship(
				std::shared_ptr<NetworkNode>& parent,
				std::shared_ptr<NetworkNode>& son) {
				
				//m_direct_parent.get();

				if (!son->m_parents.judgeInside(parent.get())) {
					std::shared_ptr<RelationNode> r1(new RelationMapNode);
					std::shared_ptr<RelationNode> r2(new RelationListNode);
					r1->set(parent, r2, &son->m_parents);
					r2->set(son, r1, &parent->m_sons);
				}

			}

			static void bindDirectRelationship(
				std::shared_ptr<NetworkNode>& parent,
				std::shared_ptr<NetworkNode>& son, double curdis) {

			
				{
						son->getSolInfo()->setMinDis(curdis);
						auto p = son->m_direct_parent;
						if(p)
						p->m_basinSize -= son->m_basinSize;
						p = son->m_direct_parent =  parent.get();
						p->m_basinSize += son->m_basinSize;

				}
			}


			void clear() {
				m_sol_info = nullptr;

				m_basinSize = 0;
				m_parents.clear();
				m_sons.clear();
				m_direct_parent = nullptr;
				m_lazy_update = false;
			}

			void initialize(const std::shared_ptr<IndividualInfo>& curinfo) {
				m_updateStamp = 0;
				m_lazy_update = false;
				m_sol_info = curinfo;
				m_parents.initialize();
				m_sons.initialize();
				m_direct_parent = nullptr;
				m_basinSize = 1;
			}
			void clearParents() {
				m_parents.clear();
			}
			void clearSons() {
				m_sons.clear();
			}

			size_t parentSize()const {
				return m_parents.size();
			}
			size_t sonSize()const {
				return m_sons.size();
			}


			void updateSonsons(ofec::Random* rnd) {
				auto sonIter = m_sons.m_head->getSharedPtr();
				int reFlag(0);
				std::shared_ptr<NetworkNode> curnode = getSharedPtr();
				while ((*sonIter->getNext())->isActive()) {
					sonIter = (*sonIter->getNext());
					sonIter->getCurInfo()->updateSons(curnode, rnd);
	
				}
			}
		};


		struct Network {

			//std::list<NetworkNode> m_network;
			//std::vector<std::shared_ptr<NetworkNode>> m_heads;


			std::vector<std::shared_ptr<NetworkNode>> m_found_peaks;
			//std::vector<std::shared_ptr<NetworkNode>> m_promising_nodes;
			std::vector<std::shared_ptr<NetworkNode>> m_promising_areas;
			std::shared_ptr<NetworkNode> m_root;
			std::shared_ptr<ofec::Random> m_rnd;
			//std::shared_ptr<ofec::Problem> m_pro;


			void updateParent(std::shared_ptr<NetworkNode>& cur, std::shared_ptr<NetworkNode>& other, ofec::Random * rnd) {
				other->updateSons(cur, rnd);
				double curdis(0);
				auto pOther = other.get();
				while (pOther != nullptr) {
					
				}
			//	reFlag = cur->getCurInfo()->calFitRel(newNode.get());

			}



			void initialize(const std::shared_ptr<ofec::Random>& rnd) {
				m_rnd = rnd;
			}


			void getGraph(
				const std::shared_ptr<NetworkNode>& root,
				std::vector<int>& belong, 
				std::vector<double>& dis2parent,
				std::vector<double>& fitness,
				std::vector<std::shared_ptr<IndividualInfo>>& sols
			) {
			
				int id = 0;
				std::map<NetworkNode*, int> nodeToId;
				std::queue<NetworkNode*> que;
				que.push(root.get());
				nodeToId[root.get()] = id++;	
				while (!que.empty()) {

					auto cur = que.front();
					que.pop();
					std::shared_ptr<RelationNode> sonIter = (*cur->sonList().getHead())->getSharedPtr();
					while ((*sonIter->getNext())->isActive()) {
						sonIter = (*sonIter->getNext());
						auto newnode = sonIter->getCurInfo();
						if (nodeToId.find(newnode) == nodeToId.end()) {
							nodeToId[newnode] = id++;
							que.push(newnode);
						}

					}
				}

				belong.resize(nodeToId.size());
				dis2parent.resize(nodeToId.size());
				fitness.resize(nodeToId.size());
				sols.resize(nodeToId.size());
				for (auto& it : nodeToId) {
					
					belong[it.second] = nodeToId[it.first->directParent()];
					dis2parent[it.second] = it.first->minDis();
					fitness[it.second] = it.first->fitness();
					sols[it.second] = it.first->getSharedSolInfo();
				}
				
			}


			int updateReltionShip(
				std::shared_ptr<NetworkNode>& better,
				std::shared_ptr<NetworkNode>& worse
			) {
				double curdis(0);
				int updateFlag = worse->calDisRel(better.get(), curdis);
				if (updateFlag != -1) {
					if (updateFlag == 1) {

						worse->clearParents();

						NetworkNode::bindRelationship(better, worse);
						NetworkNode::bindDirectRelationship(better, worse,curdis);
					}
					else {
						NetworkNode::bindRelationship(better, worse);
						if (worse->parentSize() * m_rnd->uniform.next() < 1.0) {
							NetworkNode::bindDirectRelationship(better, worse,curdis);
						}
					}

				}
				return updateFlag;
			}



			int updateReltionShip(
				std::vector<std::shared_ptr<NetworkNode>*>& better,
				std::shared_ptr<NetworkNode>& worse
			) {
				if (better.empty())return -1;

				int updateFlag = -1;
				double curdis(0);
				{
					auto& betterOne = *better.front();

					updateFlag = worse->calDisRel(betterOne.get(),curdis);

				}

				if (updateFlag != -1) {
					if (updateFlag == 1) {
						worse->clearParents();
						//NetworkNode::clearRelationshipMap(worse->m_parents);
						//worse->clearParents();

						

						for (auto& betterOneIter : better) {
							auto& betterOne = *betterOneIter;
							NetworkNode::bindRelationship(betterOne, worse);
							//worse->updateSons(betterOne);
						}
						auto& betterOne = *better[m_rnd->uniform.nextNonStd<int>(0, better.size())];
						{
							NetworkNode::bindDirectRelationship(betterOne, worse,curdis);
						}
					}
					else {
						for (auto& betterOneIter : better) {
							auto& betterOne = *betterOneIter;
							NetworkNode::bindRelationship(betterOne, worse);

						}
						//	m_rnd->uniform.next()* worse->m_parents.size();

						if (worse->parentSize() * m_rnd->uniform.next() < better.size()) {
							auto& betterOne = *better[m_rnd->uniform.nextNonStd<int>(0, better.size())];
							{
								NetworkNode::bindDirectRelationship(betterOne, worse,curdis);
							}
						}
					}

				}
				return updateFlag;

			}




			void updateTwoRelationShip(
				std::shared_ptr<NetworkNode>& node1,
				std::shared_ptr<NetworkNode>& node2
			) {
				if (node1->fitness() < node2->fitness()) {
					swap(node1, node2);
				}
				updateReltionShip(node1, node2);
			}

			void updateBest(
				std::vector<std::shared_ptr<NetworkNode>>& bestNodes,
				std::shared_ptr<NetworkNode>& cur) {

				if (bestNodes.empty()) {
					bestNodes.push_back(cur);
				}
				else {
					if (bestNodes.front()->fitness() == cur->fitness()) {
						bestNodes.push_back(cur);
					}
					else if (bestNodes.front()->fitness() < cur->fitness()) {
						bestNodes.clear();
						bestNodes.push_back(cur);
					}
				}
			}
			

			void updateNode(
				std::shared_ptr<NetworkNode>& parent1,
				std::shared_ptr<NetworkNode>& parent2,
				std::shared_ptr<NetworkNode>& betterSon,
				std::vector<std::shared_ptr<NetworkNode>>& bestNodes) {


				parent1->updateSons(betterSon, m_rnd.get());
				parent2->updateSons(betterSon, m_rnd.get());

				updateTwoRelationShip(parent1, betterSon);
				updateTwoRelationShip(parent2, betterSon);

				bestNodes.clear();


				updateBest(bestNodes, parent1);
				updateBest(bestNodes, parent2);
				updateBest(bestNodes, betterSon);

			}


			


			void updateNeighbor(std::vector<std::shared_ptr<NetworkNode>*>& totalSols, 
				std::vector<std::shared_ptr<NetworkNode>>& bestNodes) {
				bestNodes.clear();

				std::vector<std::vector<std::shared_ptr<NetworkNode>*>> parentIds(totalSols.size());
				std::vector<double> minDis(totalSols.size(), std::numeric_limits<double>::max());

				for (int idx(0); idx < minDis.size(); ++idx) {
					minDis[idx] = (*totalSols[idx])->minDis();
				}

				std::sort(totalSols.begin(), totalSols.end(), [&](
					std::shared_ptr<NetworkNode>* a, std::shared_ptr<NetworkNode>* b) {
					return (*a)->fitness() < (*b)->fitness();
				});

				for (int idx(0); idx < totalSols.size(); ++idx) {
					auto& vPId = parentIds[idx];
					auto& cursol = *totalSols[idx];
					auto& curMinDis = minDis[idx];
					for (int idy(idx + 1); idy < totalSols.size(); ++idy) {
						double curDis = cursol->distanceTo(totalSols[idy]->get());
						if (curDis < curMinDis) {
							curMinDis = curDis;
							vPId.clear();
							vPId.push_back(totalSols[idy]);
						}
						else if (curDis == curMinDis) {
							vPId.push_back(totalSols[idy]);
						}
					}
				}


				for (int idx(0); idx < totalSols.size(); ++idx) {
					updateReltionShip(parentIds[idx], *totalSols[idx]);
				}

				bestNodes.push_back(*totalSols.back());
				for (auto iter = totalSols.rbegin() + 1; iter != totalSols.rend(); ++iter) {
					if (bestNodes.front()->fitness() == (**iter)->fitness()) {
						bestNodes.push_back(**iter);
					}
					else break;
				}
			}


			void updateNeighbor(
				std::vector<std::shared_ptr<NetworkNode>>& parents,
				std::vector<std::shared_ptr<NetworkNode>>& sons,
				std::vector<std::shared_ptr<NetworkNode>>& bestNodes
			) {


				std::vector<std::shared_ptr<NetworkNode>*> totalSols;
				for (auto& it : parents) {
					totalSols.push_back(&it);
				}
				for (auto& it : sons) {
					totalSols.push_back(&it);
				}

				updateNeighbor(totalSols,bestNodes);
			}


			void updateNeighbor(
				std::vector<std::shared_ptr<NetworkNode>>& sons,
				std::vector<std::shared_ptr<NetworkNode>>& bestNodes
			) {

				std::vector<std::shared_ptr<NetworkNode>*> totalSols;
				for (auto& it : sons) {
					totalSols.push_back(&it);
				}
				updateNeighbor(totalSols, bestNodes);
			}
		};


	};
}

#endif