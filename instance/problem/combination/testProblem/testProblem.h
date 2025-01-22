/********* Begin Register Information **********
{
	"name": "ComOP_TP",
	"identifier": "ComOP_TP",
	"problem tags": [ "SOP" ]
}
*********** End Register Information **********/

#ifndef OFEC_TEST_PROBLEM_H
#define OFEC_TEST_PROBLEM_H

#include"../../../../core/problem/problem.h"
#include "../../../../core/definition.h"
#include"../selectionProblem/ttFun.h"

namespace ofec {
#define CAST_TESTPRO(pro) dynamic_cast<TestProblem*>(pr))

	class TestProblem : public Problem {
	protected:
	//	SP::Tt_fun m_funs;
	//	SP::Tt_fun_para m_par;
	//	SP::ParamTtFun m_randomPar;
		std::array<std::unique_ptr<SP::FunTt>,2> m_funs;
		std::array<std::vector<FUN_GEN::elementaryFunctionParameters>, 2> m_pars;
		int mc_T = 1000;
		SP::ParFunTt m_randomPar;
//		SP::Function3DTt m_3Dfun;

		double m_sampleRatio = 0.1;

		std::pair<double, double> m_range = { 1.0,2.0 };
		std::vector<double> m_centers;
		int m_wallNum = 20;
//		std::vector<SP::Function3DTt> m_3Dfuns;

		SP::FunTt3D m_dynamicMesh;
	public:

		//const std::pair<double, double>& x1_range() {
		//	return m_par.x_from.front();
		//}

		//const std::pair<double, double>& x2_range() {
		//	return m_pars.back().m_from_x_range;
		//}

		


		//SP::Function3DTt& get3DFun() {
		//	return m_3Dfun;
		//}
		//void generateSamples(std::vector<std::array<ofec::Real, 3>>& samples, double T, int div) const {
		//	//samples.resize(div * div);
		//	//double offset = 1.0 / double(div);
		//	//for (int idx(0); idx < div; ++idx) {
		//	//	for (int idy(0); idy < div; ++idy) {
		//	//	//	samples[idx * div + idy] = { Real(div * idx * offset),Real(div * idy * offset),Real(0) };
		//	//		samples[idx * div + idy] = { Real(idx* offset),Real(idy* offset),Real(m_3Dfun.getValue(T,  idx * offset, idy * offset)) };
		//	//	}
		//	//}  
		//}



		void generateSamples(std::vector<std::vector<std::array<ofec::Real, 3>>>& totalSamples, double T, int div) const {
			double offset = 1.0 / double(div);
			int totalDim = 20;
			totalSamples.resize(totalDim);
			double Toffset = 1.0 / double(totalDim);
			for (int wallId(0); wallId < totalDim; ++wallId) {
				auto& sample(totalSamples[wallId]);
				sample.resize(div * div);

			//	qDebug() << Toffset * double(wallId) <<"\t" << m_dynamicMesh.m_funSeed.getValAdd(1.0, Toffset * double(wallId));
				for (int idx(0); idx < div; ++idx) {
					for (int idy(0); idy < div; ++idy) {
						sample[idx * div + idy] = 
						{ Real(idx * offset),Real(idy * offset),
							Real(m_dynamicMesh.getValue(1.0,Toffset * double(wallId), idx * offset, idy * offset) + m_centers[wallId]) };
						//	samples[idx * div + idy] = { Real(div * idx * offset),Real(div * idy * offset),Real(0) };
					//	sample[idx * div + idy] = { Real(idx * offset),Real(idy * offset),Real(m_dynamicMesh.getValue(1.0,Toffset*double(wallId), idx * offset, idy * offset) + m_centers[wallId] )};
					}
				}
			}
			//totalSamples.resize(20);
			//for (int wallId(0); wallId < m_wallNum; ++wallId) {
			//	auto& sample(totalSamples[wallId]);
			//	sample.resize(div * div);
			//	for (int idx(0); idx < div; ++idx) {
			//		for (int idy(0); idy < div; ++idy) {
			//			//	samples[idx * div + idy] = { Real(div * idx * offset),Real(div * idy * offset),Real(0) };
			//			sample[idx * div + idy] = { Real(idx * offset),Real(idy * offset),Real(m_3Dfuns[wallId].getValueDivT(wallId,  idx * offset, idy * offset)) + m_centers[wallId] };
			//		}
			//	}
			//}
			//totalSamples.resize(m_wallNum);
			//for (int wallId(0); wallId < m_wallNum; ++wallId) {
			//	auto& sample(totalSamples[wallId]);
			//	sample.resize(div * div);
			//	for (int idx(0); idx < div; ++idx) {
			//		for (int idy(0); idy < div; ++idy) {
			//			//	samples[idx * div + idy] = { Real(div * idx * offset),Real(div * idy * offset),Real(0) };
			//			sample[idx * div + idy] = { Real(idx * offset),Real(idy * offset),Real(m_3Dfuns[wallId].getValue(T,  idx * offset, idy * offset))+ m_centers[wallId] };
			//		}
			//	}
			//}
		}

		void generateSamplesDivT(std::vector<std::vector<std::array<ofec::Real, 3>>>& totalSamples, double T, int div) const {
			double offset = 1.0 / double(div);
			int totalDim = 20;
			totalSamples.resize(totalDim);
			for (int wallId(0); wallId < totalDim; ++wallId) {
				auto& sample(totalSamples[wallId]);
				sample.resize(div * div);
				for (int idx(0); idx < div; ++idx) {
					for (int idy(0); idy < div; ++idy) {
						sample[idx * div + idy] =
						{ Real(idx * offset),Real(idy * offset),
							Real(m_dynamicMesh.getValue(1.0+wallId, 0.1, idx * offset, idy * offset) + m_centers[wallId]) };				}
				}
			}
		}


		virtual void initialize_()override {
			Problem::initialize_();
			m_randomPar.initialize(mc_T);
			int idx = 0;

			m_randomPar.initEdgeVarPar();
			m_funs[idx].reset(new SP::FunTt());
			m_randomPar.initPar(m_pars[idx], -1);
			m_randomPar.initFun(m_pars[idx], *m_funs[idx], m_random.get());

			m_randomPar.initEdgeNoisePar();
			idx = 1;
			m_funs[idx].reset(new SP::FunTt());
			m_randomPar.initPar(m_pars[idx], -1);
			m_randomPar.initFun(m_pars[idx], *m_funs[idx], m_random.get());

			//m_3Dfun.initialize(0, mc_T, m_sampleRatio, m_random.get());

			//m_3Dfuns.resize(m_wallNum);
			//for (auto& it : m_3Dfuns) {
			//	it.initialize(0, mc_T, m_sampleRatio, m_random.get());
			//}
			m_centers.resize(m_wallNum);
			m_centers.front() = 0;
			for (int idx(1); idx < m_centers.size(); ++idx) {
				m_centers[idx] = m_centers[idx - 1] + m_random->uniform.nextNonStd<double>(m_range.first, m_range.second);
			}

			m_randomPar.initMeshSeedFun();
			m_randomPar.initPar(m_pars[idx], -1);
			m_randomPar.initFun(m_pars[idx], m_dynamicMesh.m_funSeed, m_random.get());

			
		}


		std::array<std::unique_ptr<SP::FunTt>, 2>& get_fun(){
			return m_funs;
		}
		int getT() {
			return mc_T;
		}
		void initSolution(SolBase& s, Random *rnd) const {}
		bool same(const SolBase& s1, const SolBase& s2) const {
			return true;
		}
		Real variableDistance(const SolBase& s1, const SolBase& s2) const {
			return 0;
		}
		Real variableDistance(const VarBase& s1, const VarBase& s2) const {
			return 0;
		}

		virtual void evaluate_(SolBase& s, bool effective) {}
		
	};

	using ComOP_TP = TestProblem;

}


#endif //  TEST_PROBLEM_H