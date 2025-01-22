#ifndef OFEC_NBN_EAX_POP_H
#define OFEC_NBN_EAX_POP_H

#include "../../../../../core/algorithm/population.h"
#include "../../../../../core/problem/solution.h"


#include "../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../environment.h"
#include "nbn_learning2.h"
#include "../../nbn_alg_com/gl_calculator.h"

namespace ofec {





	class IndiEAX_info: public NBN_learn2::IndividualInfo {
	protected:
		eax_tsp::TIndi m_indi;
		
		double m_radius = 0;
		unsigned long long m_createdTime = 0;

		double m_pos = 1.0;
		
	public:

		eax_tsp::TIndi& indi() {
			return m_indi;
		}

		virtual double distance(const IndividualInfo& otherInfo)const override {
			return m_indi.distanceTo(dynamic_cast<const IndiEAX_info& >(otherInfo).m_indi);
		}
		virtual double fitness()const {
			return m_indi.fEvaluationValue* m_pos;
		}
		void setIndi(const eax_tsp::TIndi& indi) {
			m_indi = indi;
		}
		unsigned long long  getTime()const {
			return m_createdTime;
		}
		void setTime(unsigned long long time) {
			m_createdTime = time;
		}

		virtual void initialize(ofec::Problem* pro) {
			NBN_learn2::IndividualInfo::initialize();
			//	m_lastUpdateTime = 0;
			//	m_calculatedTimes = 0;
			m_pos= (pro->optMode(0) == ofec::OptMode::kMaximize) ? 1 : -1;
			m_radius = std::numeric_limits<double>::max();
		}
	};

	class PopNBN_EAX : public Population<Solution<VarVec<int>>>, public eax_tsp::TEnvironment {
	public:
		PopNBN_EAX() : Population(), TEnvironment() {}
		PopNBN_EAX(size_t size_pop, Problem* pro) : Population(size_pop, pro), TEnvironment(){};
	
		virtual ~PopNBN_EAX() = default;
		int evolve(Problem* pro, Algorithm* alg, Random* rnd) override;
		virtual void initialize(Problem* pro, Random* rnd) override;
	//	int evaluate(Problem* pro, Algorithm* alg) override;

		

		

	protected:
		NBN_learn2::Network m_net;
		unsigned long long m_maxNumStagation = 20;
		unsigned long long m_curIter = 0;
		
		GL_calculator m_gl_calculator;
		
		double m_minDis = 5;


		std::vector<eax_tsp::TIndi> m_bestSols;
		double m_bestObj = 0;
	public:

		void calNBN(
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<double>& fitness,
			std::vector<int>& optNBNid,
			std::vector<double>& optFit,
			ofec::Random* rnd
		);
		
	};
}

#endif // !OFEC_PopGL_NBN_COM_ALG_H