
#ifndef OFEC_NBN_LEARNING_H
#define OFEC_NBN_LEARNING_H

#include "../../../../../core/problem/solution.h"
#include <memory>
#include <iomanip>
namespace ofec {
	class NBN_learn {
	protected:
		struct RelationList;
		struct NetworkNode;
		
		class IndividualInfo {
		public:

			static unsigned long long  m_totalId;

			unsigned long long m_id = 0;

			int m_lastUpdateTime = 0;
			int m_calculatedTimes = 0;
			double m_minDis = 0;


			IndividualInfo() {
				m_id = m_totalId++;
			}

			virtual double distance(const IndividualInfo& otherInfo) {
				return 0;
			}
			virtual double fitness() {
				return 0;
			}

			void initialize() {
				m_lastUpdateTime = 0;
				m_calculatedTimes = 0;
				m_minDis = std::numeric_limits<double>::max();
			}
		};


		struct RelationNode : public EnableSharedPtr<RelationNode> {
			bool m_active = false;
			//NetworkNode* m_ptr_cur = nullptr;
			std::shared_ptr<NetworkNode> m_cur = nullptr;
			
			

			std::weak_ptr<RelationNode> m_parentSon;
			RelationList* m_list = nullptr;
			std::weak_ptr<RelationNode> m_Lbefore;
			std::shared_ptr<RelationNode> m_Lafter = nullptr;


			RelationNode() = default;

			//RelationNode(
			//	std::shared_ptr<NetworkNode>& cur,
			//	std::shared_ptr<RelationNode>& parentSon,
			//	RelationList& list) {
			//	set(cur, parentSon, list);
			//}


			~RelationNode() {
				clearMemory();
			}
			void initialize() {
				m_active = false;
				m_cur = nullptr;
				m_list = nullptr;
				m_parentSon.reset();
				m_Lbefore.reset();
				m_Lafter.reset();
			}


			void removeFromLink() {
				if (m_Lafter != nullptr) {
					{
						auto p = m_Lbefore.lock();
						if (p) p->m_Lafter = m_Lafter;
					}
					//m_Lbefore->m_Lafter = m_Lafter;
					m_Lafter->m_Lbefore = m_Lbefore;
					--m_list->m_numNode;
					auto id = m_cur->m_sol_info->m_id;
					m_list->m_nodeIds.erase(id);

				}
				m_Lbefore.reset();
				m_Lafter.reset();
				m_list = nullptr;
			}

			void removeNode() {

				{
					auto p = m_parentSon.lock();
					if (p) p->removeFromLink();
					m_parentSon.reset();
					this->removeFromLink();
				}
	
			}
			void clearMemory() {
				removeNode();
				{
					auto p = m_parentSon.lock();
					if(p) p->removeNode();
				}
				m_active = false;
				m_cur.reset();
				m_list = nullptr;
				m_parentSon.reset();
				m_Lbefore.reset();
				m_Lafter.reset();

			}

			

			bool set(
				const std::shared_ptr<NetworkNode>& cur,
				const std::shared_ptr<RelationNode>& parentSon,
				RelationList& list) {

				auto id = cur->m_sol_info->m_id;
				if (list.m_nodeIds.find(id) == list.m_nodeIds.end()) {
					m_cur = cur;
					m_parentSon = parentSon;
					m_list = &list;
					auto& head = list.m_head;
					m_Lafter = head->m_Lafter;
					m_Lbefore = head;
					{
						auto p = m_Lbefore.lock();
						if (p)
							p->m_Lafter = this->getSharedPtr();
					}
					//m_Lbefore->m_Lafter = this->getSharedPtr();
					m_Lafter->m_Lbefore = this->getSharedPtr();
					++m_list->m_numNode;
					m_list->m_nodeIds.insert(id);
					return true;
				}
				else return false;
			}
		};


		struct RelationList {
			std::shared_ptr<RelationNode> m_head;
			std::shared_ptr<RelationNode> m_tail;
			int m_numNode = 0;
			std::set<unsigned long long> m_nodeIds;

			
			~RelationList() {
				clear();
			}



			void initialize() {
				m_head.reset(new RelationNode);
				m_tail.reset(new RelationNode);
				m_head->m_Lafter = m_tail;
				m_tail->m_Lbefore = m_head;
				m_numNode = 0;
			}
			void clearList() {
				//m_numNode = 0;
				auto iter = m_head.get();

				while (iter->m_Lafter->m_active) {
					auto cur = iter->m_Lafter;
				//	cur->m_parentSon->removeNode();
					cur->removeNode();
					cur = nullptr;
				}
			}

			void clear() {
				clearList();
				m_head->m_Lafter.reset();
				m_tail->m_Lbefore.reset();
				m_tail = nullptr;
				m_head = nullptr;
			}


		};


		struct NetworkNode : public EnableSharedPtr<RelationNode> {
	
			std::shared_ptr<IndividualInfo> m_sol_info;
		
			//std::vector<ListNode> m_parents;
			//std::shared_ptr<NetworkNode> m_directParent;

			RelationList m_parents;
			RelationList m_sons;

			RelationList m_direct_parent;
			RelationList m_direct_sons;
			int m_basinSize = 0;

			int m_updateStamp = 0;
			bool m_lazy_update = false;


			~NetworkNode() {
				clear();
			}

			static void bindRelationship(
				 std::shared_ptr<NetworkNode>& parent, 
				 std::shared_ptr<NetworkNode>& son) {
				
				std::shared_ptr<RelationNode> r1(new RelationNode);
				std::shared_ptr<RelationNode> r2(new RelationNode);
				if (r2->set(parent, r1, son->m_parents)) {
					r1->set(son, r2, parent->m_sons);
				}

			}

			static void bindDirectRelationship(
				 std::shared_ptr<NetworkNode>& parent,
				 std::shared_ptr<NetworkNode>& son) {

				son->removeDirectParent();

				std::shared_ptr<RelationNode> r1(new RelationNode);
				std::shared_ptr<RelationNode> r2(new RelationNode);
				r1->set(son, r2, parent->m_direct_sons);
				r2->set(parent, r1, son->m_direct_parent);

				parent->m_basinSize += son->m_basinSize;
			}


			void clear() {
				m_sol_info = nullptr;

				m_basinSize = 0;
				m_parents.clear();
				m_sons.clear();
				m_direct_parent.clear();
				m_direct_sons.clear();
				m_lazy_update = false;
			}

			void initialize(const std::shared_ptr<IndividualInfo>& curinfo) {
				m_updateStamp = 0;
				m_lazy_update = false;
				m_sol_info = curinfo;
				m_parents.initialize();
				m_sons.initialize();
			}

			void clearParents() {
				m_parents.clearList();
			}
			void removeDirectParent() {
				auto& parent = m_direct_parent.m_head->m_Lafter;
				if (parent->m_active) {
					parent->m_cur->m_basinSize -= m_basinSize;
					m_direct_parent.clearList();
				}

			}

		};


		struct Network {
			
			//std::list<NetworkNode> m_network;
			//std::vector<std::shared_ptr<NetworkNode>> m_heads;
			

			std::vector<std::shared_ptr<NetworkNode>> m_roots;
			std::vector<std::shared_ptr<NetworkNode>> m_promising_nodes;
			std::vector<std::shared_ptr<NetworkNode>> m_history_peak;

			std::shared_ptr<ofec::Random> m_rnd;



			void updateSons(
				std::shared_ptr<NetworkNode>& near, 
				std::shared_ptr<NetworkNode>& newNode) {
				auto iter = near->m_sons.m_head.get();

				int updateFlag = -1;
				double curDis = 0;
				while (iter->m_Lafter->m_active) {
					updateFlag = -1;
					{
						auto& cur = iter->m_Lafter;
						if (cur->m_cur->m_sol_info->fitness() < newNode->m_sol_info->fitness()) {
							updateFlag = updateReltionShip(newNode, cur->m_cur);
						}
					}

					if (updateFlag == -1) {
						iter = iter->m_Lafter.get();
					}
				}

			}


			int updateReltionShip(
				std::shared_ptr<NetworkNode>& better,
				std::shared_ptr<NetworkNode>& worse
				) {
				int updateFlag = -1;
				
				double curDis = better->m_sol_info->distance(*worse->m_sol_info);
				if (curDis < worse->m_sol_info->m_minDis) {
					updateFlag = 1;
				}
				else if (curDis == worse->m_sol_info->m_minDis) {
					updateFlag = 2;
				}

				if (updateFlag != -1) {
					if (updateFlag == 1) {
						worse->clearParents();

						NetworkNode::bindRelationship(better, worse);
						NetworkNode::bindDirectRelationship(better, worse);
					}
					else {
						NetworkNode::bindRelationship(better, worse);
						if (worse->m_parents.m_numNode * m_rnd->uniform.next() < 1.0) {
							NetworkNode::bindDirectRelationship(better, worse);
						}
					}

				}
				return updateFlag;
			}


			void updateTwoRelationShip(
				std::shared_ptr<NetworkNode>& node1,
				std::shared_ptr<NetworkNode>& node2
			) {
				if (node1->m_sol_info->fitness() < node2->m_sol_info->fitness()) {
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
					if (bestNodes.front()->m_sol_info->fitness() == cur->m_sol_info->fitness()) {
						bestNodes.push_back(cur);
					}
					else if (bestNodes.front()->m_sol_info->fitness() < cur->m_sol_info->fitness()) {
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

				updateSons(parent1, betterSon);
				updateSons(parent2, betterSon);

				//std::vector<
				

				//std::shared_ptr<NetworkNode> head = 

				updateTwoRelationShip(parent1, betterSon);
				updateTwoRelationShip(parent2, betterSon);
				
				bestNodes.clear();
				

				updateBest(bestNodes, parent1);
				updateBest(bestNodes, parent2);
				updateBest(bestNodes, betterSon);

			}
			
			
			void updateNeighbor(
				std::vector<std::shared_ptr<NetworkNode>>& parents,
				std::vector<std::shared_ptr<NetworkNode>>& sons
			) {
				
			}
		};
		
		std::shared_ptr<ofec::Problem> m_pro;
	};
}

#endif