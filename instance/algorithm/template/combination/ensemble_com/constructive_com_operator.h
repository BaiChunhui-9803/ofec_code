#ifndef CONSTRUCTIVE_COM_POP_H
#define CONSTRUCTIVE_COM_POP_H

#include<memory>
#include<vector>
#include"../../../../../core/definition.h"

namespace ofec {

	template<typename TInterpreter>
	class ConstructiveComOp {

	protected:
		bool m_flag_init_pop = false;
		bool m_flag_step_init = false;
		bool m_flag_step_next = false;
		bool m_flag_step_final = false;
		bool m_flag_update_final = false;


		Real m_init_pop_pro = 1.0;
		Real m_step_init_pro = 1.0;
		Real m_step_next_pro = 1.0;
		Real m_step_final_pro = 1.0;
		Real m_update_final_pro = 1.0;


	public:
		using IndividualType = typename TInterpreter::SolutionType;

		ConstructiveComOp() = default();
		~ConstructiveComOp() = default();

		bool flagInitPop() { return m_flag_init_pop; }
		bool flagStepInit() { return m_flag_step_init; }
		bool flagStepNext() { return m_flag_step_next; }
		bool flagStepFinal() { return m_flag_step_final; }
		bool flagUpdateFinal() { return m_flag_update_final; }



		Real getInitPopPro() { return m_init_pop_pro; }
		Real getStepInitPro() { return m_step_init_pro; }
		Real getStepNextPro() { return m_step_next_pro; }
		Real getStepFinalPro() { return m_step_final_pro; }
		Real getUpdateFinalPro() { return m_update_final_pro; }


		void setInitPopPro(Real pro) { m_init_pop_pro = pro; }
		void setStepInitPro(Real pro) { m_step_init_pro = pro; }
		void setStepNextPro(Real pro) { m_step_next_pro = pro; }
		void setStepFinalPro(Real pro) { m_step_final_pro = pro; }
		void setUpdateFinalPro(Real pro) { m_update_final_pro = pro; }



		virtual void initialize() {}
		virtual bool initializePop(std::vector<std::unique_ptr<IndividualType>>& pops,const TInterpreter& interperter, Algorithm *alg,Problem *pro, Random *rnd);
		virtual bool stepInit(IndividualType& sol,const std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter, Algorithm *alg, Problem *pro, Random *rnd) {return false;}
		virtual void stepBack(const IndividualType& sol, const std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter, Algorithm *alg, Problem *pro, Random *rnd) {}
		virtual bool stepNext(IndividualType& sol, const std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter, Algorithm *alg, Problem *pro, Random *rnd) { return false; }
		virtual bool stepFinal(IndividualType& sol, const std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter, Algorithm *alg, Problem *pro, Random *rnd) {return false;}
		virtual bool updateFinal(std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter, Algorithm *alg, Problem *pro, Random *rnd) { return false; }
		virtual void constructiveFinal(const std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter, Algorithm *alg, Problem *pro, Random *rnd) {}
	};
	template<typename TInterpreter>
	inline bool ConstructiveComOp<TInterpreter>::initializePop(std::vector<std::unique_ptr<IndividualType>>& pops, const TInterpreter& interperter, Algorithm *alg, Problem *pro, Random *rnd)
	{
		for (int i = 0; i < m_individuals.size(); ++i) {
			m_individuals[i]->initialize(i, pro, rnd);
		}
		return true;
	}

	
	
}






#endif