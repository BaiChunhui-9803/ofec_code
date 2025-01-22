#include "nbn_eax_alg_v4.h"



#ifdef OFEC_DEMO
#include <ui/custom/buffer/algorithm/combination/buffer_nbn_sols.h>
#include <ui/custom/buffer/algorithm/combination/density_matrix/buffer_density_matrix.h>
#include <core/global_ui.h>
//#include "../../../problem/combination/travelling_salesman/travelling_salesman.h"

#endif
void ofec::NBN_EAX_Alg_V4::run_()
{
	while (!terminating()) {
		m_pop->evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
		updateBuffer();
#endif
	}
}

void ofec::NBN_EAX_Alg_V4::initialize_()
{
	Algorithm::initialize_();
	m_pop.reset(new PopNBN_EAX_V4);
	m_pop->initialize(m_problem.get(), m_random.get());
}

#ifdef OFEC_DEMO
void ofec::NBN_EAX_Alg_V4::updateBuffer()
{

	auto& info = ofec_demo::BufferComNBNAlg::m_communicate_info;
	m_pop->calNBN(info.m_belong, info.m_dis2parent, info.m_fitness, info.m_bestIds, info.m_bestFit, m_random.get());
	m_pop->getPopIds(info.m_curPop);
	ofec_demo::g_buffer->appendAlgBuffer(this);
}
#endif