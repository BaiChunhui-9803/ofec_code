/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 31 Jan 2024 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// �����ӿռ������������ͼ

#ifndef OFEC_SPACEGRAPH_H
#define OFEC_SPACEGRAPH_H

namespace ofec {

	class SpaceNode {
	public:
		SpaceNode(size_t inx) : m_space_inx(inx) {}

		// ���ǰ���ھӽڵ�
		void addFrontNeighbor(const SpaceNode& neighbor) {
			m_front_neighbors.emplace_back(std::make_shared<SpaceNode>(neighbor));
		}

		// ��Ӻ����ھӽڵ�
		void addBehindNeighbor(const SpaceNode& neighbor) {
			m_behind_neighbors.emplace_back(std::make_shared<SpaceNode>(neighbor));
		}

		// ��ȡǰ���ھӽڵ�
		std::vector<std::shared_ptr<SpaceNode>> getFrontNeighbors() {
			return m_front_neighbors;
		}

		// ��ȡ�����ھӽڵ�
		std::vector<std::shared_ptr<SpaceNode>> getBehindNeighbors() {
			return m_behind_neighbors;
		}

		// ��ȡ�����ӿռ�����
		std::vector<Real> getBehindNeighborDists() {
			return m_spacenode_dists;
		}

		size_t getSpaceInx() {
			return m_space_inx;
		}

	private:
		size_t m_space_inx;//�ӿռ�����
		std::vector<std::shared_ptr<SpaceNode>> m_front_neighbors;
		std::vector<std::shared_ptr<SpaceNode>> m_behind_neighbors;
		std::vector<Real> m_spacenode_dists;//�����ھӵļ����С����
	};

	class SpaceEdge {
	public:
		//SpaceEdge(size_t inx) : m_space_inx(inx) {}

		// ���ǰ��ڵ�
		void addFrontNode(const SpaceNode& neighbor) {
			m_front_node = std::make_shared<SpaceNode>(neighbor);
		}

		// ��Ӻ���ڵ�
		void addBehindNode(const SpaceNode& neighbor) {
			m_behind_node = std::make_shared<SpaceNode>(neighbor);
		}

		// ��ȡǰ���ھӽڵ�
		std::shared_ptr<SpaceNode> getFrontNode() {
			return m_front_node;
		}

		// ��ȡ�����ھӽڵ�
		std::shared_ptr<SpaceNode> getBehindNode() {
			return m_behind_node;
		}

		// ��ȡ�ߵ�ֵ
		Real getEdgeValue() {
			return m_edge_value;
		}

		// ���ñߵ�ֵ
		void setEdgeValue(Real v) {
			m_edge_value = v;
		}

	private:
		size_t m_space_inx;//�ӿռ�����
		std::shared_ptr<SpaceNode> m_front_node;
		std::shared_ptr<SpaceNode> m_behind_node;
		Real m_edge_value;//�ߵ�ֵ
	};

	// ����ͼ��
	class SpaceDirectedGraph {
	public:
		// ��ӽڵ�
		void addNodes(size_t space_inx, size_t space_inx2) {
			SpaceNode node(space_inx);
			SpaceNode behind_node(space_inx2);
			node.addBehindNeighbor(behind_node);
			m_graph_nodes.insert(std::make_pair<>(space_inx, std::make_shared<SpaceNode>(node)));
		}

		// ��ӽڵ�
		void addNode(size_t space_inx) {
			SpaceNode node(space_inx);
			//SpaceNode behind_node(space_inx2);
			//node.addBehindNeighbor(behind_node);
			m_graph_nodes.insert(std::make_pair<>(space_inx, std::make_shared<SpaceNode>(node)));
		}

		// ��������
		void addEdge(SpaceNode& from_node, const SpaceNode& to_node, Real v) {
			from_node.addBehindNeighbor(to_node);
			SpaceEdge edge;
			edge.addFrontNode(from_node);
			edge.addBehindNode(to_node);
			edge.setEdgeValue(v);
			m_graph_edges.emplace_back(std::make_shared<SpaceEdge>(edge));
		}

		// ��������
		void addEdge(size_t space_inx1, size_t space_inx2, Real v) {
			//SpaceNode from_node(space_inx1);
			SpaceNode to_node(space_inx2);
			m_graph_nodes[space_inx1]->addBehindNeighbor(to_node);
			SpaceEdge edge;
			edge.addFrontNode(*m_graph_nodes[space_inx1]);
			edge.addBehindNode(to_node);
			edge.setEdgeValue(v);
			m_graph_edges.emplace_back(std::make_shared<SpaceEdge>(edge));
		}

		// ��ӡͼ���ڽ��б�
		void printGraph() {
			for (const auto& entry : m_graph_nodes) {
				std::shared_ptr<SpaceNode> current = entry.second;
				std::cout << "Node " << current->getSpaceInx() << " is connected to:";
				for (auto neighbor : current->getBehindNeighbors()) {
					std::cout << " " << neighbor->getSpaceInx();
				}
				std::cout << std::endl;
			}
		}

		std::map<size_t, std::shared_ptr<SpaceNode>> getGraphNode(){return m_graph_nodes;}
		std::vector<std::shared_ptr<SpaceEdge>> getGraphEdge() { return m_graph_edges; }

	private:
		std::map<size_t, std::shared_ptr<SpaceNode>> m_graph_nodes;  // ���ڴ洢ͼ�еĽڵ�
		std::vector<std::shared_ptr<SpaceEdge>> m_graph_edges;  // ���ڴ洢ͼ�еı�
	};
}

#endif  // !OFEC_SPACEGRAPH_H