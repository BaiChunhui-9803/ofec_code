#include "nbn_alg_com.h"

#ifdef OFEC_DEMO
#include <ui/custom/buffer/algorithm/combination/buffer_nbn_sols.h>
#include <ui/custom/buffer/algorithm/combination/density_matrix/buffer_density_matrix.h>
#include <core/global_ui.h>
#include "../../../problem/combination/travelling_salesman/travelling_salesman.h"

#endif

void ofec::NBN_COM_ALG::run_()
{

	m_pop.reset(new PopGL_NBN_COM_ALG(m_pop_size, m_problem.get()));
	m_pop->initialize(m_problem.get(), m_random.get());
	while (!terminating()) {
		m_pop->evolve(m_problem.get(), this, m_random.get());
#ifdef OFEC_DEMO
		updateBuffer();
#endif
	}
}

void ofec::NBN_COM_ALG::initialize_()
{

	Algorithm::initialize_();
	m_pop_size = 1e3;
	m_pop.reset();
}

#ifdef OFEC_DEMO
void ofec::NBN_COM_ALG::updateBuffer()
{

	auto& info = ofec_demo::BufferComNBNAlg::m_communicate_info;
	m_pop->getSolInfo(info.m_belong, info.m_dis2parent, info.m_fitness);
	std::vector<std::array<int, 2>> optEdge;
	if (CAST_TSP(m_problem.get())->optima()->isSolutionGiven()) {
		const auto& optSol = CAST_TSP(m_problem.get())->optima()->solution(0).variable();
		CAST_TSP(m_problem.get())->transferEdgeSol(optSol, optEdge);
	}
	
	auto& densityMat = ofec_demo::BufferDensityMatrix::ms_density_matrix.front();
	auto& densityMatRange = ofec_demo::BufferDensityMatrix::ms_density_range.front();
	auto& frameMat = ofec_demo::BufferDensityMatrix::ms_frame_matrix.front();
	
	auto& mat = m_pop->getProMat();
	std::vector<int> sortIds(m_problem->numVariables());
	for (int idx(0); idx < sortIds.size(); ++idx) {
		sortIds[idx] = idx;
	}
	std::vector<int> originIds2Id(sortIds);
	densityMat.resize(mat.size());
	//densityMatRange.resize(mat.size());
	frameMat.resize(mat.size());
	for (int idx(0);idx<mat.size();++idx) {
		densityMat[idx].resize(mat[idx].size());
		frameMat[idx].resize(mat[idx].size());
	}
	densityMatRange.first = std::numeric_limits<double>::max();
	densityMatRange.second = -std::numeric_limits<double>::max();

	std::vector<bool> optEdgeInfo;
	for (int idx(0); idx < mat.size(); ++idx) {
		auto& curpro = mat[idx];
		optEdgeInfo.resize(mat[idx].size());
		std::fill(optEdgeInfo.begin(), optEdgeInfo.end(), false);
	
		if (!optEdge.empty()) {
			for (auto& curEdge : optEdge[idx]) {
				optEdgeInfo[curEdge] = true;
			}
		}
		sort(sortIds.begin(), sortIds.end(), [&](int a,int b) {
			return curpro[a] > curpro[b];
		});

		auto& showPro = densityMat[idx];
		auto& showFrame = frameMat[idx];
		for (int idy(0); idy < sortIds.size(); ++idy) {
			showPro[idy] = curpro[sortIds[idy]];
			showFrame[idy] = optEdgeInfo[sortIds[idy]];
		}
		double maxPro(0), minPro(0);
		calMax(curpro, maxPro);
		calMin(curpro, minPro);
		densityMatRange.first = std::min<double>(minPro, densityMatRange.first);
		densityMatRange.second = std::max<double>(maxPro, densityMatRange.second);

	



	}

	ofec_demo::g_buffer->appendAlgBuffer(this);


}
#endif