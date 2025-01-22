#include "nbn_eax_pop_v5.h"


#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../../../../utility/nbn_visualization/nbn_fla.h"
#include "../../../../../utility/typevar/typevar.h"



ofec::PopNBN_EAX_V5::Indi* ofec::PopNBN_EAX_V5::RelationNode::getCur() {
	return m_cur.m_cur;
}
ofec::PopNBN_EAX_V5::Indi* ofec::PopNBN_EAX_V5::RelationNode::getParent() {
	return m_cur.m_other;
}
void ofec::PopNBN_EAX_V5::RelationNode::removeFromLink() {

	if(m_list !=nullptr){
		m_before->m_after = m_after;
		m_after->m_before = m_before;
		m_after = m_before = nullptr;
		--m_list->m_num;
	}
}




void ofec::PopNBN_EAX_V5::RelationNode::insertToLink(RelationList* list) {

	list->insertNode(this);
	//m_list = list;
}


void ofec::PopNBN_EAX_V5::RelationList::insertNode(RelationNode* cur) {

	cur->removeFromLink();
	cur->m_after = m_head->m_after;
	cur->m_before = m_head.get();
	cur->m_list = this;
	m_head->m_after->m_before = cur;
	m_head->m_after = cur;
	++m_num;
	//m_sumDis2parent += cur.m_cur->m_dis2parent;
}

void ofec::PopNBN_EAX_V5::RelationList::removeNode(RelationNode* cur)
{
	cur->removeFromLink();
}


void ofec::PopNBN_EAX_V5::RelationList::clearLink() {
	auto head = m_head.get();
	while (head->m_after->isActive()) {
		head->m_after->removeFromLink();
	}

}
