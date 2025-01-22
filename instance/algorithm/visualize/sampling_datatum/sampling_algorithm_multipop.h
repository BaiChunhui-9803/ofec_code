/*************************************************************************
* Project:Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Yiya Diao
* Email: diaoyiyacug@gmail.com
* Language: C++
*-------------------------------------------------------------------------------
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*-------------------------------------------------------------------------------
*
*  Modified by Diao Yiya on 2024 09 16
*
*-------------------------------------------------------------------------------
*************************************************************************/

#ifndef OFEC_SAMPLING_ALGORITHM_MULTIPOP_H
#define OFEC_SAMPLING_ALGORITHM_MULTIPOP_H

#include "../../../../core/algorithm/algorithm.h"
#include "../../../../core/environment/environment.h"
#include "../../../../core/parameter/parameter_variant.h"

#include "../../../../datum/algorithm/multi_pop.h"
#include "../../../../instance/algorithm/continuous/single_objective/global/cma_es/cma_es.h"
namespace ofec {

#define CAST_SAMPLE_DATA(alg) dynamic_cast<SamplingData*>(alg)

	class SamplingData : virtual public Algorithm {
		OFEC_ABSTRACT_INSTANCE(SamplingData)
	public:

		virtual ~SamplingData() = default;

		struct SolutionInfo:public EnableSharedPtr<SolutionInfo>{
			int m_runId = 0;
			int m_iter = 0;
			int m_popId = 0;
			int m_indiId = 0;

			virtual ~SolutionInfo() = default;

			// for output to file 
			virtual void toParameterVariants(ParameterVariantStream& out)const {
				out << m_runId << m_iter << m_popId << m_indiId;
			}
			virtual void fromParameterVariants(ParameterVariantStream& in) {
				in >> m_runId >> m_iter >> m_popId >> m_indiId;
			}

		};


	protected:
		void clear() {
			m_sols.clear();
			m_solInfos.clear();
		}
	protected:
		void addInputParameters() {}

	public:
		const std::vector<std::shared_ptr<SolutionBase>>& getSols() {
			return m_sols;
		}
		const std::vector<std::shared_ptr<SolutionInfo>>& getSolInfos() {
			return m_solInfos;
		}

	protected:
		std::vector<std::shared_ptr<SolutionBase>> m_sols;
		std::vector<std::shared_ptr<SolutionInfo>> m_solInfos;


	};


	template<typename T>
	class SamplingAlgorithm : virtual public SamplingData, virtual public T {
		OFEC_CONCRETE_INSTANCE(SamplingAlgorithm) 
	protected:
		int m_iteration = 0;
	protected:
		void addInputParameters() {}

		virtual void initialize_(Environment* env) override {
			SamplingData::clear();
			T::initialize_(env);
			m_iteration = 0;
		}

		virtual void run_(Environment* env) override {
			T::run_(env);
		}

		virtual void handleDatumUpdated(Environment* env) override {
			if (g_multi_pop.updateFlag()) {
				auto pro = env->problem();
				++m_iteration;
				SolutionInfo curInfo;
				curInfo.m_iter = m_iteration;
				curInfo.m_popId = 0;
				curInfo.m_indiId = 0;
				auto& pops = g_multi_pop.pops;
				for (int idpop(0); idpop < pops.size(); ++idpop) {
					curInfo.m_popId = idpop;
					for (int idIndi(0); idIndi < pops[idpop].size(); ++idIndi) {
						curInfo.m_indiId = idIndi;
						m_sols.emplace_back(pro->createSolution(*pops[idpop][idIndi]));
						m_solInfos.emplace_back(new SolutionInfo(curInfo));
					}

				}
			}
		}

	public:

		virtual void reset() override {
			T::reset();
		}
		virtual bool terminating()override {
			return T::terminating();
		}


		virtual void datumUpdated(Environment* env, DatumBase& datum) override {
			T::datumUpdated(env, datum);
		}


	};
}

#endif