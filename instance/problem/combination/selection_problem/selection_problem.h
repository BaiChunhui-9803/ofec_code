/********* Begin Register Information **********
{
	"name": "ComOP_DSP",
	"identifier": "ComOP_DSP",
	"problem tags": [ "SP", "ComOP", "DOP", "NoisyOP" ]
}
*********** End Register Information **********/


#ifndef SELECTION_PROBLEM_H
#define SELECTION_PROBLEM_H

#include "../../../../core/global.h"
#include "../../../../core/problem/problem.h"
#include "../../../../core/problem/optima.h"
#include "../../../../core/problem/domain.h"
#include "../../../../core/problem/solution.h"
#include "../../../../utility/memory_record/memory_record.h"
#include "../../../../utility/linear_algebra/vector.h"
#include "../../../../utility/functional.h"
#include "../../../../core/problem/uncertainty/dynamic.h"
#include "../../../../core/problem/uncertainty/noisy.h"
#include "../../../../instance/record/robust_sp/rcr_selection_problem_dynamic.h"
#include "ttFun.h"
#include "ttFunPar.h"
#include "sp_map.h"
#include <string>
#include <map>
#include <array>
#include"../metric_com_uncertianty.h"
//#include "../../instance/algorithm/visualize/graph_calculator.h"

namespace ofec {
#define CAST_DYN_SP(pro) dynamic_cast<sp::SelectionProblem*>(pro)

	namespace sp {
		class SelectionProblem : virtual public Noisy, virtual public Dynamic {
		private:
			
			static thread_local SolutionSPInfo m_temp_solinfo;
		public:
			using solution_type = MapInfo::solutionType;
			using variable_type = MapInfo::variableType;
			//	static const Real m_eps ;
		protected:
			// flag 
			bool m_flag_points_correlation = false;
			// for memory
			MemeryRecord m_memory_record;
			//	void init_record();
			//	void record_solution(const SolutionBase& sol);
			ParaInfo m_para;
			int m_number_variables;
			MapInfo m_map;

			//			int m_change_fre = 1000
			int m_time_span = 7;
			Real m_survival_time = 0.4;

			std::vector<std::pair<std::string, Real>> objectives;
			std::unique_ptr<SolutionBase> m_curX = nullptr;
			std::vector<std::pair<std::string, Real>> beforeObjectives;
			std::unique_ptr<SolutionBase> m_beforeX = nullptr;
			std::vector<std::string> m_heads;
			std::ofstream out;


			const int m_sample_times = 50;
			const int m_future_eval_times = 100;
			std::vector<Real> m_info_sol_objectives;
			std::vector<Real> m_info_sol_future_objective;
			std::vector<std::string> m_info_heads;
			std::string m_info_file_name = "total_info";
			std::ofstream m_info_out;

			std::string m_future_effective_file_name = "future_effective";
			std::ofstream m_future_effective_out;
			std::vector<Real> m_record_effective_vals;
			std::vector<int> m_decision_dim;
			std::vector<int> m_working_sol;

			std::vector<std::vector<std::vector<std::pair<int, double>>>> m_opt_nodes;
			Real m_dijstra_objecitive = 0;
			solution_type m_dijstra_sol;
			std::vector<std::vector<std::vector<SolutionSPInfo>>> m_opt_nodes_info;
			std::vector<double> m_dim_min_obj;
			Real m_near_opt_epsion = 0;

			//std::vector<std::vector<std::vector<itn>>

			// for test
			std::map<int, int> m_same_sols;
			std::vector<std::vector<int>> m_each_best_sols;
			std::vector<double> m_each_best_errs;
			// near opt
			std::vector<solution_type> m_near_opt_sols;



			// for multi opts
			// parameters
			double m_near_opt_accu = 0.005;
			//double m_max_near_opt_dis = 0;
			std::vector<std::pair<std::vector<int>, double>> m_near_opts;

			std::vector<int> m_neighghbor_size;
			bool m_variableMemoryChange = false;

			int m_test_idRand = -1;

			

			
			// cal multi optima
			const int m_max_opts = 1e2;
			std::vector<int> m_CodingMapToRealIdxs;
			long long m_multi_optima_running_time;
			bool m_flag_multi_optima = true;

		    
			std::string m_local_optima_filename = "local_optima";
			
			std::string m_file_type = ".sp_local_opt";
			std::string m_file_running_time = "local_optima_running_time.csv";


		//	int m_test_T = 0;
		protected:


			void setShowFitness(SolutionBase& sol) const{
				sol.setFitness(-sol.objective()[0]);
			}
			
			
			void outputMultiOpt();
			void updateOptimas();
			
			void insertOptmFileInside(int T, std::vector<std::shared_ptr<solution_type>>& sols)const;
			void insertOptimasFile(int T, std::vector<std::shared_ptr<solution_type>>& record_sols);
			void insertOptimasFile(int T);
			void insertOptimasFiles(std::pair<int,int> from_to_T);
			void insertOptimasFileThreads();

			void setOptimalStruct(int T, std::vector<std::shared_ptr<SolutionBase>>& sols);

			void analyseDynOptDis();

			std::string getDirName() const {
				std::string filename = g_working_dir;
				filename += "/instance/problem/combination/selection_problem/";
				return std::move(filename);
			}

			


			void calInfo(const SolutionBase& sol);
			void setInfoHead();
			void opepInfoFile();
			void insertInfoFile();
			void updateNeighborSize() {
				auto& station(m_map.getStations());
				m_neighghbor_size.resize(m_decision_dim.size());
				for (int idx(0); idx < m_decision_dim.size(); ++idx) {
					m_neighghbor_size[idx] = station[m_decision_dim[idx]].size();
				}
			}

			bool changeDecisionDims() {
				std::vector<int> temp_decision_dim;
				for (int idx(0); idx < m_map.numWall(); ++idx) {
					if (!m_map.getScenesDim(getT(m_T))[idx]) {
						temp_decision_dim.push_back(idx);
					}
				}
				if (temp_decision_dim != m_decision_dim) {
					m_decision_dim = temp_decision_dim;
					m_number_variables = m_decision_dim.size();
					updateNeighborSize();
					return true;
				}
				else {
					m_number_variables = m_decision_dim.size();
					return false;
				}

			}
			void setT(int T) {
				m_T = T;

			}
			double getT(int T)const {
				return double(T) / double(m_para.mc_T);
			}
			//Real calValue(int intT, const SolutionBase& s, Random *rnd);
			Real getPriceValue(int intT, const SolutionBase& s, Random *rnd)const;
			//Real getMeanPriceValue(int T , const SolutionBase& s);
			// mu& sigma : f(x)+\delta ,  f(x,t+delta)
			void evaluate_soltuion(std::vector<Real>& mus, std::vector<Real>& sigmas,
				double T, const SolutionBase& s,
				bool mean = false, int sample_num = 1, int future_t = 1);
			// mu& sigma :  f(x)+delta ,  f(x,t+delta),   f(x+\delta,t)
			void evaluate_solution_effVal(
				std::vector<Real>& mus,
				std::vector<Real>& sigmas,
				double T, const SolutionBase& s,
				bool mean = false, int sample_num = 1, int future_t = 1);

			void addSpInfoPos(SolutionSPInfo& info,
				const std::vector<int>& x, int to_dim, int intT, Random *rnd)const;
			void addSpInfoPosNew(SolutionSPInfo& info,
				const std::vector<int>& x, int to_dim, int intT, Random *rnd)const;


			// for test
			void testStatic(const std::vector<int>& sol,
				int pos, int T, Random *rnd)const;
			void testStatic(int T, Random *rnd)const;
			void testSubSolution(
				SolutionSPInfo& temp,
				std::array<double, 3>& dis_info,
				const std::vector<int>& sol,
				int pos, int T, Random *rnd)const;
			double testSubSolution(
				const std::vector<int>& sol,
				int from,
				int total_test,
				int T, Random *rnd
			)const;
		public:


			virtual void updateLocalOptimaInfo()const override;


			virtual void mapToVD01(const SolutionBase& sol, std::vector<double>& vals) {
				auto& x = dynamic_cast<const solution_type&>(sol).variable().vect();
				vals.resize(x.size());
				for (int idx(0); idx < x.size(); ++idx) {
					vals[idx] = double(x[idx]) / double(m_map.numCandidates(idx));
				}
			}

			


			virtual double getWorstObj(int numSample);
			//virtual void getDynamicObjective(const SolutionBase& s, std::vector<Real>& objs, int t = 0) {
			//	objs.resize(1);
			//	objs[0] = getPriceValue(t, s);
			//}


			virtual const std::vector<int>& getNeighborSize() {
				return m_neighghbor_size;
			}

			const std::vector<std::pair<std::vector<int>, double>>& getNearOpts() {
				return m_near_opts;

			}
			//const std::vector<solution_type>& getNearOpts() {
			//	return m_near_opt_sols;
			//}
			void recordSol(const SolutionBase& sol);

			const MapInfo& getMap()const {
				return m_map;
			}
			ParaInfo& getPara() {
				return m_para;
			}
			// T,t range[0,1], x,y range[-2,2]
			void getMeshSeed(int wall_id, 
				double T, double t,
				std::array<double, 2>& vals
			) {
				double mu(0), sigma(0);
				m_map.getWall()[wall_id].
					m_mesh_price.getDistribution(
						vals[0], vals[1], T, t, vals[0], vals[1]);
			}
			// T,t range[0,1], x,y range[-2,2]
			void getMeshVal(int wall_id, double T, double t,
				double x,double y,
				std::array<double,2>& vals)const {
				double mu(0), sigma(0);
				m_map.getWall()[wall_id].
					m_mesh_price.getDistribution(
						x,y, T, t, vals[0], vals[1]);
			}
			// wall id pos id
			size_t numPoints(int dim) const {
				return m_map.numCandidates(m_decision_dim[dim]);
			}
			virtual size_t numberVariables() const override {
				return m_number_variables; 
			}
			//double getT(int T) const {
			//	return double(T) / double(m_para.m_dynamic_info.m_curT);
			//}
			virtual int updateEvaluationTag(SolutionBase& s, Algorithm *alg, bool effective_eval) override;
			virtual void change() override {
				if (m_flag_change) {
					
					Dynamic::change(); 
					updateOptimas();
					if (m_opt_base.number() == 0) {
						int stop = -1;
					}
					++m_T;
					m_variableMemoryChange = changeDecisionDims();
					//updateNeighborSize();
					//std::cout << "change" <<"\t"<<m_T<< std::endl;
					//std::cout << "change fre" << "\t" << m_frequency << std::endl;
					//std::cout << "evaluation " << GET_ALG(this).evaluations() << std::endl;
					m_record_effective_vals.clear();



	
				//	m_info_out <<"change" << std::endl;
					//dijstra(m_T, -1);
					//if (m_curX != nullptr) {
					//	
					//	auto& originS = dynamic_cast<const solution_type&>(*m_curX);
					//	m_beforeX.reset(new solution_type(originS));
					//	beforeObjectives = objectives;
					//	//for (auto& it : beforeObjectives) {
					//	//	out << it.first << ",";
					//	//}
					//	for (auto& it : beforeObjectives) {
					//		out << it.second << ",";
					//	}
					//	out << std::endl;
					//	//calInfo(*m_beforeX);
					//	//insertInfoFile();
					//}
				}

			}

			virtual void updateCandidates(const SolutionBase& sol, std::list<std::unique_ptr<SolutionBase>>& candidates) const override {
				if (candidates.empty()) {
					candidates.emplace_back(new solution_type(dynamic_cast<const solution_type&>(sol)));
				}
				else if (sol.objective().front() == candidates.front()->objective()[0]) {
					candidates.emplace_back(new solution_type(dynamic_cast<const solution_type&>(sol)));
				}
				else if (sol.objective().front() < candidates.front()->objective().front()) {
					candidates.clear();
					candidates.emplace_back(new solution_type(dynamic_cast<const solution_type&>(sol)));
				}
			}

#ifdef OFEC_DEMO
			virtual void appendBuffer(Algorithm *alg, bool effective_eval)const override;
#endif

			// 通过 problem 继承
			virtual bool same(const SolutionBase& s1, const SolutionBase& s2) const override;
			virtual Real variableDistance(const SolutionBase& s1, const SolutionBase& s2) const override;
			virtual Real variableDistance(const VariableBase& s1, const VariableBase& s2) const override;
			virtual void initializeSolution(SolutionBase& s, Random *rnd) const override;
			virtual void initializeSolution(SolutionBase& s, Random *rnd) const override;

			virtual void setCurSolution(const SolutionBase& sol)override;

			// init solution around center within radius

			virtual bool initializeSolution(SolutionBase& s, const SolutionBase& center, double radius, Random *rnd)const;
			Real getHeuristicInfo(int dim, int from, int to) {
				return m_map.get_heuristic_info(m_T, dim, from, to);
			}
			// sample solutions
			void generate_level_samples(std::vector<solution_type> & samples, Random *rnd);
			template<typename SolutionType>
			void generate_HLS_samples(const SolutionBase& s,std::vector<std::unique_ptr<SolutionBase>>& samples, Random *rnd)
			{
				auto& x = dynamic_cast<const solution_type&>(s).variable();
				auto& obj = dynamic_cast<const solution_type&>(s).objective();
				int sample_id(0);
				std::vector<int> dim_sample(m_number_variables, 0);
				generate_samples_numbers(dim_sample, samples.size(), rnd);
				std::vector<int> pos_sample;
				for (int dim(0); dim < m_number_variables; ++dim) {
					pos_sample.resize(numPoints(dim));
					generate_samples_numbers(pos_sample, dim_sample[dim], rnd);
					for (int pos_idx(0); pos_idx < pos_sample.size(); ++pos_idx) {
						while (pos_sample[pos_idx]--) {
							samples[sample_id].reset(new SolutionType());
							SolutionType& cur(dynamic_cast<SolutionType&>(*samples[sample_id]));
							cur.variable() = x;
							cur.objective().resize(1);
							cur.variable()[dim] = pos_idx;
							++sample_id;
						}
					}
				}
			}

			void evaluate_inside(SolutionBase& s,int T) const;
			virtual void evaluate_(SolutionBase& s, bool effective) override;
			// metric 
			Real evaluate_mean_value(const SolutionBase& s) ;
			Real evaluate_effective_mean_value(const SolutionBase& s);
			Real evaluate_mean_value_sample(const SolutionBase& s,int sample_num=1000);
			Real evaluate_effective_mean_value_sample(const SolutionBase& s,int sample_num=1000);

			virtual void setHeads();
			virtual void showInfomations(Algorithm *alg)override;
			virtual void printfSolution(Algorithm *alg,const SolutionBase& sol);
			// solution 
			virtual void outSolution(std::ostream& out, const std::unique_ptr<SolutionBase>& sol) ;
			virtual bool inSolution(std::istream& in, std::unique_ptr<SolutionBase>& sol) ;
			// test problem
			void testProblemCharacer();

		protected:
			void dijstra(int T, Random *rnd);
			void dijstraEachNodeOpt(int T, Random *rnd);
			void dijstraEachNodeOptTest(int T, Random *rnd);
			void dijstraNearOpt(int T, Random *rnd)const;
			void dijstraEachNodeOpt(
				std::vector<std::vector<double>>& min_dis_toStart,
				std::vector<std::vector<double>>& min_dis_toEnd,
				std::vector<std::vector<std::vector<int>>>& best_sols,
				int T, Random *rnd)const;
			void dijstraNearOpt(
				std::vector<std::pair<std::vector<int>,double>>& sols,
				double near_opt_accu,
				int T, Random *rnd)const;


			void dijstraMultiOpt(
				std::vector<std::unique_ptr<solution_type>>& sols
			)const;

			void dijstraMultiOptThreads(
				std::vector<std::unique_ptr<solution_type>> &sols
			)const;

			void outputMultiOptThreads(const std::string& filename)const;


			void updateSolType(
				int idDim,
				const std::vector<std::vector<std::unique_ptr<solution_type>>>& before_opts,
				const std::vector<std::vector<SolutionSPInfo>>& before_opts_info,
				std::vector<std::vector<std::unique_ptr<solution_type>>>& cur_opts,
				std::vector<std::vector<SolutionSPInfo>>& cur_opts_info,
				const std::vector<std::array<int,2>> & sol_idxs,
				std::pair<int, int> from_to
			)const ;
			
			void updateSolSolsBestDis(
				std::vector<std::vector<std::unique_ptr<solution_type>>>& cur_opts,
				std::vector<std::vector<SolutionSPInfo>>& cur_opts_info,
				const std::vector<int> &sortedC,
				std::pair<int, int> from_to
			)const;

			

		protected:
			virtual void initialize_()override;
			virtual void copy(const Problem& rP)override;
		
		protected:
			// generate random sample times 
			std::vector<int> generate_samples_numbers(int size, int sample_num, Random *rnd);
			void generate_samples_numbers(std::vector<int>& idx_samples, int sample_num, Random *rnd);

			// 通过 Noisy 继承
			virtual void getNoisyObjective(const SolutionBase& s, std::vector<Real>& objs, Random *rnd = -1) override;

			// 通过 Dynamic 继承
			virtual void getDynamicObjective(const SolutionBase& s, std::vector<Real>& objs, int t = 0) override;
};
	}
	using ComOP_DSP = sp::SelectionProblem;
}



#endif