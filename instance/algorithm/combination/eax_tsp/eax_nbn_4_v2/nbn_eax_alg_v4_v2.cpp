#include "nbn_eax_alg_v4_v2.h"



#ifdef OFEC_DEMO
#include <ui/custom/buffer/algorithm/others/nbn_algorithm/buffer_nbn_algorithm_data.h>
#include <core/global_ui.h>
//#include "../../../problem/combination/travelling_salesman/travelling_salesman.h"

#endif
void ofec::NBN_EAX_Alg_V4_V2::run_()
{
	while (!terminating()) {
		m_pop->evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
		updateBuffer();
#endif
	}
}

void ofec::NBN_EAX_Alg_V4_V2::initialize_()
{
	Algorithm::initialize_();
	m_pop.reset(new PopNBN_EAX_V4_V2);
	m_pop->initialize(m_problem.get(), m_random.get());
}

#ifdef OFEC_DEMO
void ofec::NBN_EAX_Alg_V4_V2::updateBuffer()
{

	auto& info = ofec_demo::BufferNBNAlgData::m_communicate_info;
	m_pop->calNBN(info.nbn_data, info.nbn_data_best, info.m_curPop, info.curBasin, info.m_curSearchArea, info.m_bestIds, m_problem.get(),m_random.get());
//	m_pop->getPopIds(info.m_curPop);
	ofec_demo::g_buffer->appendAlgBuffer(this);
}
#endif