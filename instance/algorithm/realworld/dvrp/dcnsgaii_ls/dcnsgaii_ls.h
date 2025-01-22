/********* Begin Register Information **********
{
	"name": "DCNSGAII-LS",
	"identifier": "DCNSGAII_LS",
	"problem tags": [ "DVRP" ]
}
*********** End Register Information **********/

#ifndef OFEC_DCNSGAII_LS_H
#define OFEC_DCNSGAII_LS_H

#include "../../../../../core/algorithm/population.h"
#include "../int_dvrp_pop/int_dvrp_pop.h"
#include "ga_dvrp.h"

namespace ofec {
	class DCNSGAII_LS_pop : public Population<DCMOEA_ind<Solution<DVRP::routes>>>, 
							public DCMOEA<DCMOEA_ind<Solution<DVRP::routes>>> {
	public:
		using Population<DCMOEA_ind<Solution<DVRP::routes>>>::IndividualType;
	public:
		DCNSGAII_LS_pop(size_t size_pop, Problem *pro);
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
		void initialize(Problem *pro, Random *rnd) override;
		void initializeAfterEvaluation();
		void sort();
		const size_t get_num_vio_obj() { return m_num_vio_obj; }
		void Online(Problem *pro, Algorithm *alg, Random *rnd);
		void set_bestInd(const DCMOEA_ind<Solution<DVRP::routes>> &best_ind) { m_best_ind_Offline = best_ind; }
		const DCMOEA_ind<Solution<DVRP::routes>> &get_bestInd() { return m_best_ind_Online; }
		
		/*void set_m_cnt_fail() { m_cnt_fail++; }
		void clear_m_cnt_fail() { m_cnt_fail = 0; }
		void set_m_convergence() { m_convergence = true; }
		int &get_m_cnt_fail() { return m_cnt_fail; }*/
	protected:
		void select_next_parent_population(std::vector<IndividualType*>& pop, Problem *pro);
		//void reduce_radius();
		//void caculate_nichecount(std::vector<IndividualType*>& pop);
		Dominance e_Pareto_compare(IndividualType* const&s1, IndividualType* const&s2);
		Dominance Pareto_compare(IndividualType* const&s1, IndividualType* const&s2);
		void generate_offspring_LS(Problem *pro, Algorithm *alg, Random *rnd);
		void generate_offspring_ALS(Problem *pro, Algorithm *alg, Random *rnd);
		void calculate_cp();
		void tsp_nearest(std::vector<size_t>& member_index, Real presentTime, Problem *pro, Random *rnd);
		void get_cost(size_t s, size_t e, Real &present_time, Problem *pro, Random *rnd);
		void ls_Rswap2(IndividualType &offspring, Problem *pro, Random *rnd);
		void ls_tsp(IndividualType &offspring, Problem *pro, Random *rnd);
		void ls_maxDT(IndividualType &offspring, Problem *pro);
		void ls_r1(IndividualType &offspring, Problem *pro, Random *rnd);
		void ls_rMulti(IndividualType &offspring, Problem *pro, Random *rnd);
		void ls_2opt(IndividualType &offspring, Problem *pro, Random *rnd);
		void ls_maxWTlength(IndividualType &offspring, Problem *pro);
		void ls_maxLength(IndividualType &offspring, Problem *pro, Random *rnd);
		bool check_order(IndividualType &offspring, Problem *pro);

	protected:
		INIT_DVRP_POP m_ahc_tsp;
		//GA_DVRP m_GA;
		std::vector<std::unique_ptr<IndividualType>> m_offspring;
		std::vector<std::unique_ptr<IndividualType>> m_temp_pop;
		//DE::recombine_strategy m_recombine_strategy;
		std::vector<int> m_rand_seq;
		
		bool m_evo_success = false;
		size_t m_num_vio_obj = 1;
		int m_cnt_fail = 0;
		std::vector<int> m_score;
		std::vector<Real> m_p_select;
		std::vector<Real> m_cp;//cumulative probability
		bool m_convergence = false;
		int m_num_operators = 8;
		
		
		//Online
		IndividualType m_best_ind_Online;
		IndividualType m_best_ind_Offline;
		std::vector<std::unique_ptr<IndividualType>> m_best_inds;
		bool m_insert_success = false;
		bool m_evo_success_Online = false;
		int m_cnt_fail_Online = 0;
		int m_OnlineIters = 2000;
		std::vector<int> m_score_Online;
		std::vector<Real> m_p_select_Online;
		std::vector<Real> m_cp_Online;//cumulative probability
		bool m_convergence_Online = false;
		int m_num_operators_Online = 3;
		void construct_path_Online(IndividualType & ind, size_t i, Problem *pro);
		void updateSeqOnline(IndividualType & ind, size_t i, Problem *pro, Algorithm *alg, Random *rnd);
		void updateSeqOnline_inserted(IndividualType & ind, size_t i, std::pair<int,int>, Problem *pro, Algorithm *alg, Random *rnd);
		void Online_ALS(IndividualType & ind, size_t i, Problem *pro, Algorithm *alg, Random *rnd);
		void ls_Rswap2_Online(IndividualType &ind, size_t i, Problem *pro, Random *rnd);
		void ls_maxDT_Online(IndividualType &ind, size_t i, Problem *pro);
		void ls_2opt_Online(IndividualType &ind, size_t i, Problem *pro, Random *rnd);
		//void remove_nOrder(IndividualType & ind, size_t i, size_t j);
	public:
		std::string m_algName_ = "ALSDCMOEA";
		std::string m_proName_ = "VRPRTC";
		std::string m_initialName_ = "AHC_initial";
	};

	class DCNSGAII_LS : public Algorithm {
	protected:
		std::unique_ptr<DCNSGAII_LS_pop> m_pop;
		size_t m_pop_size;
		int m_local_num = 100;
		bool m_isOnline = false;

		void initialize_() override;
		void run_() override;
	public:
		DCNSGAII_LS(param_map& v);
		void initialize();
		void run_();
		void record();
#ifdef OFEC_DEMO
		void updateBuffer() override{}
#endif
	};
}



#endif // !OFEC_DCNSGAII_LS_H
