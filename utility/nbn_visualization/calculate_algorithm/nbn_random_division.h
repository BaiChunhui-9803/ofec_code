

#ifndef NBN_RANDOM_SAMPLE_H
#define NBN_RANDOM_SAMPLE_H

#include <vector>
#include <memory>
#include <list>
#include "../../core/problem/solution.h"
#include "../../utility/typevar/typevar.h"
#include "nbn_division_base.h"


namespace ofec {
	class NBN_RandomDivision : public NBN_DivisionBase{
		struct node {
			int m_id = -1;
			std::shared_ptr<SolBase> m_sol;
			double m_fitness;
			int m_belong_id = -1;
			double m_dis2parent = 1e9;
			bool m_flag_opt = false;
			int m_popIter = -1;
			int m_popSolId = -1;
			int m_algId = -1;

			node() = default;
			node(int id, const std::shared_ptr<SolBase>& sol, double fitness,
				int belong_id, double dis2parent,
				bool flag_opt=false,
				int popIter = -1, int popSolId = -1, int algId=-1 )
				: m_id(id), m_sol(sol), m_fitness(fitness),
				m_belong_id(belong_id), m_dis2parent(dis2parent),
				m_flag_opt(flag_opt)
				, m_popIter(popIter), m_popSolId(popSolId),m_algId(algId){}
			//std::unique_ptr<node> m_next = nullptr;
		};
		std::list<node> m_list;

		//void addNode(int id, 
		//	const std::shared_ptr<SolBase>& sol,
		//	double fitness, int belong_id, double dis2parent, bool flag_opt=false
		//	);


	protected:
		virtual void initialize_(bool flag_grid_sample = true)override;
		void addRandomSols(int num_sample);
	public:

		NBN_RandomDivision() = default;
		~NBN_RandomDivision() = default;

		virtual size_t size()const  override{
			return m_list.size();
		}
		void addRandomSol(const SolBase& sol,int & id, bool flag_opt = false, 
			int popIter = -1, int popSolId = -1, int algId=-1);

		virtual void addSol(const SolBase& new_sol, int& belong_id, bool flag_opt = false,
			int popIter = -1, int popSolId = -1, int algId = -1) override {
			addRandomSol(new_sol, belong_id, flag_opt, popIter, popSolId, algId);
		}

		void getSharedData(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong)const;


		void getSharedData(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<int> popIters,
			std::vector<int> popSolIds,
			std::vector<int> algIds
			)const;



		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<bool>& flagOpt
		)const override;


		virtual void getSharedNBN(
			std::vector<std::shared_ptr<SolBase>>& sols,
			std::vector<double>& fitness,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			std::vector<int> popIters,
			std::vector<int> popSolIds,
			std::vector<int> algIds
		)const override;
	};
}



#endif