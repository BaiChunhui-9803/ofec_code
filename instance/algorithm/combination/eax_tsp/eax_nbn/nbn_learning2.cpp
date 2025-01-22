#include "nbn_learning2.h"


//unsigned long long  ofec::NBN_learn2::IndividualInfo::m_totalId = 0;


namespace ofec {



	void NBN_learn2::RelationMapNode::set(
		std::shared_ptr<NetworkNode>& cur,
		std::shared_ptr<RelationNode>& parentSon,
		RelationShip* curStruct
	) {
		RelationNode::set(cur, parentSon, curStruct);
		m_cur = cur.get();
		m_pMap = dynamic_cast<RelationMap*> (curStruct);
		m_pMap->insert(this);
	}

	void NBN_learn2::RelationMapNode::removeFromStruct()
	{
		if (m_pMap != nullptr) {
			m_pMap->remove(m_cur);
			m_pMap = nullptr;
		}
	}

	NBN_learn2::RelationShip* NBN_learn2::RelationMapNode::getRelationShip() {
		return m_pMap;
	}
	void NBN_learn2::RelationListNode::set(
		std::shared_ptr<NetworkNode>& cur,
		std::shared_ptr<RelationNode>& parentSon,
		RelationShip* curStruct
	) {
		RelationNode::set(cur, parentSon, curStruct);
		m_cur = cur;
		m_list = dynamic_cast<RelationList*>(curStruct);
		m_list->insert(this);
		//m_pMap = dynamic_cast<RelationMap*> (curStruct);
	}



	NBN_learn2::RelationShip* NBN_learn2::RelationListNode::getRelationShip() {
		return m_list;
	}

	void NBN_learn2::RelationListNode::setAfterNode(std::shared_ptr<RelationNode>& afterNode)  {
		m_after = afterNode;
	}
	void NBN_learn2::RelationListNode::setBeforeNode(RelationNode* beforeNode) {
		m_before = beforeNode;
	}

	void NBN_learn2::RelationListNode::removeFromStruct()
	{
		this->m_after->setBeforeNode(this->m_before);
		this->m_before->setAfterNode(this->m_after);
		this->m_after.reset();
		this->m_before = nullptr;
		this->m_list->m_size--;
	}

	void NBN_learn2::RelationList::insert(RelationListNode* cur)
	{
		auto cur_ptr = cur->getSharedPtr();
		cur_ptr->setAfterNode(m_head->m_after);
		cur_ptr->setBeforeNode(m_head.get());
		//cur_ptr->m_after = m_head->m_after;
		//cur_ptr->m_before = m_head.get();
		m_head->m_after->setBeforeNode(cur_ptr.get());
		m_head->m_after = cur_ptr;
		++m_size;
	}

}