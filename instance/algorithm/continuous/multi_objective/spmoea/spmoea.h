/********* Begin Register Information **********
{
	"name": "SPMOEA",
	"identifier": "SPMOEA",
	"problem tags": [ "OOMOP", "ConOP", "MOP" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yong Xia
* Email: changhe.lw@google.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.

*************************************************************************/
// Created: 20 Feb 2023 by Qingshan Tan (email:qingshan.t@cug.edu.cn)
// Last modified: 

#ifndef OFEC_SPMOEA_H
#define OFEC_SPMOEA_H

#include "../../../../../core/algorithm/multi_population.h"
#include "../../../../../core/problem/continuous/continuous.h"

//#include "../../../template/selection/multi_objective/nsgaii.h"
#include "../../../template/classic/ga/sbx_pop.h"
//#include "../../../template/classic/de/population.h"
#include "../../../../algorithm/continuous/multi_objective/moea_de_pop.h"
#include "../../../template/classic/de/individual.h"

#include "../../../../record/multi_objective/rcr_vec_real_moea.h"
#include "../../../../record/rcr_vec_real.h"
#include "../metrics.h"

#include "../../../../../utility/kd-tree/kdtree_space.h"
#include "../../../../../utility/nondominated_sorting/fast_sort.h"
#include "../../../../../utility/environment_selection/selection_methods.h"
#include "../../../../../utility/functional.h"
#include "../../../../../utility/linear_algebra/matrix.h"
#include "../../../../../utility/metricsMOP/IGD.h"

#include "../mo_hlc.h"

#include<fstream>
#include<tuple>
#include<Eigen/Core>
#include<Eigen/Dense>

#ifdef OFEC_DEMO
#include <core/global_ui.h>
#endif

namespace ofec {

	class Subspace {
	private:
		std::vector <std::pair<Real, Real>> m_bound;
		std::vector<std::shared_ptr<Solution<>>> m_history_inds;//子空间内的历史解
		std::vector<size_t> m_subspace_front_inx;  //子空间前沿解在历史解中的索引

	public:
		Subspace(const std::vector<std::pair<Real, Real>>& bound) :m_bound(bound) {}
		std::vector<std::pair<Real, Real>> getBox() { return m_bound; }
		void setBoxBound(const std::vector<std::pair<Real, Real>>& bound) { m_bound = bound; }
		std::vector<std::shared_ptr<Solution<>>>& getSubspaceHisSols() { return m_history_inds; }
		std::vector<size_t>& getSubspaceFrontInx() { return m_subspace_front_inx; }
		void updateSubspaceInfo(Population<Solution<>>& new_pop, Problem* pro);
		void clearInfos() {
			m_history_inds.clear();
			m_subspace_front_inx.clear();
		}
	};

	class SPMOEA_pop :public Population<Solution<>> {
	public:
		explicit SPMOEA_pop(size_t size_pop, Problem *pro);
		void initialize(Problem *pro, Random *rnd) override;
		int evolve(Problem *pro, Algorithm *alg, Random *rnd) override;
		void envirSelection(Problem *pro, Random *rnd, bool b);
		void envirSelection(std::vector<std::pair<Real,Real>> &bound, Problem *pro, Random *rnd);
		void envirSelection(Population<Solution<>>& pop, Population<Solution<>>& compound_pop,Problem* pro);
		void NDSort(Population<Solution<>>& pop,Problem* pro);
		void evalDense(Population<Solution<>>& pop, Population<Solution<>>& pop_combined,Problem* pro);
		std::vector<size_t> denseSelect(Population<Solution<>>& pop, Population<Solution<>>& pop_combined, Problem* pro);
		Population<Solution<>>& getCombinPop() { return m_pop_combined; }
		Population<Solution<>>& getOffspring() { return m_offspring; }
		void setSearchRange(std::vector<std::pair<Real, Real>>& v) { m_pop_search_range = v; }
		std::vector<std::pair<Real, Real>> getSearchRange() { return m_pop_search_range; }
		std::string getPopState() { return m_pop_state; }
		void setPopState(std::string ss) { m_pop_state=ss; }
		void setLocateSubspace(std::vector<size_t> v) { m_space_inx = v; }
		std::vector<size_t> getLocateSubspace() { return m_space_inx; }
		//void addPopDistribute(std::vector<std::vector<Real>>& s) { m_pop_distribute.emplace_back(s); }
		std::vector<std::vector<std::vector<Real>>>& getPopDistribute() { return m_pop_distribute; }
		void updatePopDistribute(Problem *pro);
		void updatePopHisInfo(Population<Solution<>>& new_pop, Problem *pro);
		void addPopdist();
		std::vector<std::vector<Real>>& getPopdist() { return m_pop_dist; }
		bool popConverged();
		std::vector<std::pair<Real, Real>> getPopSearchBox() { return m_search_box; }
		void setSearchBox(std::vector<std::pair<Real, Real>>& box) { m_search_box = box; }
		std::vector<std::shared_ptr<Solution<>>>& getPopHisFrontSols() { return m_pop_history_front_sols; }
		std::vector<std::shared_ptr<Population<Solution<>>>>& getPopGenFrontSols() { return m_pop_gen_front_sols; }

	private:
		std::string m_pop_state = "active";//子种群状态
		std::vector<size_t> m_space_inx;//子种群所在子空间索引
		size_t m_stag_gen = 3;//停滞代数
		Population<Solution<>> m_pop_combined;  // combination of parent and children
		Population<Solution<>> m_offspring;  // 
		std::vector<std::pair<Real, Real>> m_pop_search_range;//单个种群的搜索范围
		std::vector<std::pair<Real, Real>> m_pop_cur_obj_range;//单个种群当前的目标值范围
		std::vector<std::pair<Real, Real>> m_pop_his_front_range;//单个种群历史前沿的目标值范围
		std::vector <std::shared_ptr<Solution<>>> m_pop_history_front_sols;//种群的历史非支配解
		std::vector <std::shared_ptr<Solution<>>> m_pop_archive;//从历史非支配解中选择的部分代表解
		std::vector <std::shared_ptr<Solution<>>> m_pop_history_sols;//种群的历史评估解
		std::vector <std::shared_ptr<Population<Solution<>>>> m_pop_gen_front_sols;//种群中每一代的新的前沿解
		std::vector<std::vector<std::vector<Real>>> m_pop_distribute;//记录种群在搜索空间和目标空间的分布状态
		std::vector<std::vector<Real>> m_pop_dist;//记录种群与前一代的相似度

		std::vector<std::pair<Real, Real>> m_search_box;
	};

	typedef std::vector<std::vector<std::map<size_t, std::vector<size_t>>>> pop_attach;//目标和搜索空间的依附
	typedef std::vector<std::map<size_t, std::vector<size_t>>> space_attach;//种群在某一空间的依附，前排和非前排子空间
	using KDTree = nanoflann::KDTreeSpace<Real>;

	class SPMOEA :public Algorithm {
	public:
		void initialize_() override;
		void initPop(Problem *pro, Algorithm *alg, Random *rnd);
		void initiVarSpace(Problem *pro);
		void record() override;

		int evolve(Problem *pro, Algorithm *alg, Random *rnd);
		virtual void generateOffspring(Problem *pro, Random *rnd);
		void repairSol(std::vector<Real>& sol, std::vector<std::pair<Real, Real>>& sub_bound, Random *rnd);
		void repairSol(std::vector<Real>& sol, std::vector<size_t> space_inx, Random *rnd);
		bool normalSol(std::vector<Real>& sol, std::vector<std::pair<Real, Real>>& sub_bound);

		std::vector<std::vector<Real>> sampleByDE(SPMOEA_pop& pop, std::vector<std::pair<Real, Real>>& bound, size_t off_num, size_t k, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleByDE(SPMOEA_pop& pop, size_t base_inx,std::vector<std::pair<Real, Real>>& bound, size_t off_num, size_t k, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleByGA(SPMOEA_pop& pop, std::vector<std::pair<Real, Real>>& bound, size_t off_num, Problem *pro, Random *rnd);
		std::vector<std::vector<Real>> sampleByExtend(size_t space1, size_t space2, size_t off_num, Problem *pro, Random *rnd);
		std::vector<std::vector<Real>> sampleByNeighSpace(size_t space1, size_t off_num, Problem *pro, Algorithm* alg, Random *rnd);
		std::vector<std::vector<Real>> sampleByFrontNeighSpace(size_t space1, std::vector<size_t>& front_link_spaces, std::vector<std::vector<size_t>>& manifold_dist, size_t sample_num, Problem *pro, Algorithm* alg, Random *rnd);
		std::vector<std::vector<Real>> sampleByFrontNeighSpace(size_t space1,size_t ind_inx, std::vector<size_t>& front_link_spaces, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg,Random *rnd);
		std::vector<std::vector<Real>> sampleByNeighSpaceFront(size_t space1, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t off_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, std::vector<Real>& sol, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		virtual std::vector<std::vector<Real>> sampleInFrontNeighSpace(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInFrontNeighSpace2(std::vector<size_t>& front_link_spaces, size_t ind_inx, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleByLinkSpaces(std::vector<size_t>& front_spaces, std::vector<std::vector<size_t>>& manifold_dist, size_t sample_num, Problem *pro, Algorithm* alg, Random *rnd);
		std::vector<std::vector<Real>> sampleByManifoldSpaces(std::vector<size_t>& front_link_spaces, std::vector<std::vector<size_t>> &manifold_dist, size_t sample_num, Problem *pro, Algorithm* alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInLinkSpaces(std::vector<size_t>& front_spaces, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleBetweenSpaces(std::vector<size_t> &front_link_space1, std::vector<size_t>& front_link_space2, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg,Random *rnd);
		std::vector<std::vector<Real>> sampleInSpace(size_t parent_space, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro, Algorithm *alg,Random *rnd);
		std::vector<std::vector<Real>> sampleInSpace(size_t parent_space, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInFrontSpaces(std::vector<size_t>& front_spaces, std::vector<std::pair<Real, Real>>& bound, size_t sample_num, Problem *pro,Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleBySolution(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInSolution(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleByImprove(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		virtual std::vector<std::vector<Real>> samplePushorExt(Solution<>& sol, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		virtual std::vector<std::vector<Real>> samplePush(Solution<>& sol, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleInSubspace(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleBySubspace(Solution<>& sol, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleInSubspace(size_t space1, size_t ind_inx, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInSubspace(size_t space_inx, size_t sample_num, Problem* pro, Algorithm* alg, Random* rnd);
		std::vector<std::vector<Real>> sampleInFrontSubspace(size_t space1, size_t ind_inx, std::vector<size_t>& front_link_spaces, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInRange(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleInRange(Solution<>& sol, size_t sample_num, const std::vector<std::pair<Real, Real>>& bound, Problem *pro, Algorithm *alg, Random *rnd);
		std::vector<std::vector<Real>> sampleByNeighGradient(Solution<>& sol, size_t sample_num, Problem *pro, Algorithm *alg, Random *rnd);

		void manifoldSpread(std::vector<size_t>& front_link_spaces,Problem *pro,Random *rnd);

		virtual void updateInteractiveSols();

		void setPopInitRange(SPMOEA_pop &pop, const std::vector<std::pair<Real, Real>>& pre_bound,const std::vector<std::pair<Real,Real>> &new_bound,Problem *pro);

		std::vector<Real> vectorOperator(std::vector<Real>& sol1, std::vector<Real>& sol2, Real offset, Problem* pro, Algorithm* alg, Random *rnd);
		std::vector<Real> vectorPredict(std::vector<Real>& sol1, std::vector<Real>& sol2, std::vector<std::pair<Real, Real>>& bound, Random *rnd);

		virtual void PopResourceAssign(std::vector<size_t>& assign_pop_resource, size_t switch_period, Problem* pro);

		void updateRegionLinkSpace(size_t inx);
	    virtual void updateFrontRegionLinkSpace(Problem* pro, Random* rnd);
		std::vector<std::tuple<size_t, size_t, Real>> findCloseSpaces(std::vector<std::vector<size_t>>& clusters, Random *rnd);
		std::vector<std::vector<size_t>> calSpaceManifoldDist(std::vector<size_t>& spaces);

		bool ifLinearLinkSpace(size_t inx, std::vector<size_t> clustered_inx1, size_t inx2,Problem* pro, Random* rnd);
		
		bool ifLinearLinkSpace(size_t inx1, size_t inx2);

		Real ifContinuousSpace(size_t inx1, size_t inx2);

		std::map<size_t,Real> spaceLinkProbability(size_t inx1, std::vector<size_t> neis);
		std::map<size_t, Real> frontSpaceLinkProbability(size_t inx1);

		Real ifSingleSegSpace(size_t inx1);
		std::pair<bool,std::pair<size_t,Real>> ifMultiSegSpace(size_t inx1);

		std::pair<int, Real> findSplitDimAndPos(int inx, Problem* pro);

		std::vector<std::vector<std::vector<size_t>>>& getRegionLinkSpace() { return m_region_link_space; }
		std::vector<std::vector<size_t>>& getFrontRegionLinkSpace() { return m_front_region_link_space; }
		std::vector<std::vector<size_t>> clusterRegionFrontSpace(std::vector<size_t>& frontspace);
		virtual std::vector<std::vector<size_t>> linearClusterFrontSpace(std::vector<size_t>& frontspace,Problem* pro, Random* rnd);
		std::vector<size_t> getNonlapNeighSpaces(std::vector<size_t>& spaces, std::vector<size_t>& cur_inx, std::vector<size_t>& dealt_spaces);

		void updateVarSpaceInfo(Population<Solution<>>& pop, Problem *pro, Random *rnd);
		void updateSubspaceFrontSol(Population<Solution<>>& pop, Problem *pro, Random *rnd);
		void updateVarSpaceRank(Problem* pro, Random* rnd);
		void updateVarSpaceDominanceRank(Problem* pro, Random* rnd);

		void updateHisVarSpaceInfo(std::vector<std::shared_ptr<Solution<>>>& pop, Problem *pro, Random *rnd);

		std::vector<size_t> getFrontSubspace();
		std::vector<size_t> getBehindSubspace();

		/*
		线性回归
		*/
		std::vector<Real> linearRegression(const std::vector<std::vector<Real>> &input_data, Real learn_ratio, int iterate_times);
		void gradient_descent(const std::vector<std::vector<Real>>& training_set, int feature_num, std::vector<Real>& w, Real learn_ratio, int iterator_time);
		Real Theta(const std::vector<std::vector<Real>>& training_set, int featue_num, const std::vector<Real>& w);
		Real predict(const std::vector<Real>& weigh, const std::vector<Real>& data, int feature_num);



		//void updateVarSpaceInfo(Solution<>& pop);
		//void updateHistoryFrontSol(Population<Solution<>>& pop, Problem *pro);
		//void updateHistorySols(Population<Solution<>>& pop);
		void updateHistoryInfo(Population<Solution<>>& new_pop,Problem *pro);
		void updateArchive(size_t num, Problem *pro);
		void updateObjRange(Population<Solution<>>& pop,Problem *pro);
		void updateObjSpace();
		void updatePopDistribute(size_t inx, Problem *pro);
		void updatePopdist(size_t inx);
		void updateNewPop(Population<Solution<>>& new_pop);
		const Population<Solution<>>& getNewPops() const { return m_new_pops; }

		void generatePop(const std::vector<std::vector<size_t>>& space_clusters, Problem *pro, Algorithm *alg, Random *rnd);

		virtual int subspaceSeperable(size_t inx,Problem *pro);//检测子空间是否具有可分性
		void divideSubspace(size_t inx,size_t num);
		void splitSubspace(size_t inx, int dim,Real pos);

		bool updateCluster();
		virtual void clusterSubspace();
		bool popConverged(size_t inx);
		void updateFrontSpace();

		void NDSort(Population<Solution<>>& pop);
		size_t layerNDSort(Population<Solution<>>& pop);
		void recordMetrics(Problem *pro, Algorithm *alg);
		//void recordHisFront(Problem *pro);

		space_attach spaceAttach(const Population<Solution<>>& pop, size_t num_spaces,Problem *pro, bool b);//b选择目标空间还是搜索空间
		space_attach spaceAttach(const std::vector <std::shared_ptr<Solution<>>>& pop,size_t num_spaces, Problem *pro,bool b);//b选择目标空间还是搜索空间
		pop_attach popAttach(const Population<Solution<>>& pop, size_t num_spaces, Problem *pro);
		pop_attach popAttach(const std::vector <std::shared_ptr<Solution<>>>& pop, size_t num_spaces, Problem *pro);
		
		/*环境选择*/
		std::vector<size_t> selectIndiFromFront(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd);
		std::vector<size_t> selectIndiFromFront(const std::vector <std::shared_ptr<Solution<>>>& pop, size_t select_num, Problem *pro);
		std::vector<size_t> selectMaxMinFromFront(const Population<Solution<>>& pop, size_t select_num);
		std::vector<size_t> selectMaxMinFromFront(const std::vector <std::shared_ptr<Solution<>>>& pop, size_t select_num);
		std::vector<size_t> selectMaxMinFromFront(const std::vector <std::vector<Real>>& data, size_t select_num);
		std::vector<size_t> selectMaxMinFromFront(const std::vector<std::vector<Real>>& data, std::vector<size_t> ranks, size_t select_num);
		std::vector<size_t> selectByNorm(const std::vector<std::vector<Real>>& data, size_t select_num,Problem *pro,Random *rnd);

		std::vector<size_t> selectIndiFromSpace(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd);
		std::vector<size_t> selectFromSpace(const Population<Solution<>>& pop, size_t select_num, Problem *pro, Random *rnd);
		std::vector<size_t> selectFromSpace(const std::vector<std::vector<Real>>& pop, std::vector<size_t> rank, size_t select_num, Problem *pro, Random *rnd);
		virtual std::vector<size_t> selectIndiFromSpace(const std::vector <std::shared_ptr<Solution<>>>& pop, size_t select_num, Problem *pro,Random *rnd);

		std::vector<size_t> multiReferSelection();
		std::vector<size_t> neighLocalSelection(size_t neighbor_num,Problem* pro);

		void multiObjSelection(Problem* pro);
		void sparseSelection(Problem* pro);
		void localSelection(Problem* pro);
		void localCrowdSelection(Problem* pro, Random* rnd);
		void localCrowdSelection2(Problem* pro, Random* rnd);
		void localCrowdSelection3(Problem* pro, Random* rnd);
		void localCrowdSelection4(Problem* pro, Random* rnd);
		void localCrowdSelection5(size_t select_num, Problem* pro, Random* rnd);
		void localIndSelection(Problem* pro, Random* rnd);
		void doubleCrowdSelection(Problem* pro);
		void indicatorSelection(Problem* pro, Random* rnd);
		void subspaceSelection(Problem* pro, Random* rnd);
		void ensembleSelection(size_t select_num, Problem* pro, Random* rnd);
		virtual void SolutionSelection(Problem* pro, Random* rnd);
		void knnSelection(Problem* pro, Random* rnd);

		void visitCount(const std::vector<std::vector<Real>>& sols, std::vector<Real>& middle_point);

		Real spaceCoverRatio(size_t inx,size_t subspace_num);
		bool ifFrontChanged(std::vector<size_t> &pre_front_spaces);


		void setVarPartitionNum(const size_t v) { m_num_region_var = v; }
		void setObjPartitionNum(const size_t v) { m_num_region_obj = v; }
		size_t numVarRegion() const { return m_num_region_var; }
		size_t numObjRegion() const { return m_num_region_obj; }
		std::vector<std::pair<Real, Real>> getPopRange() { return m_pop_obj_range; }
		std::vector<std::pair<Real, Real>> getFrontObjRange() { return m_front_obj_range; }
		std::vector<std::pair<Real, Real>> getHisObjRange() { return m_his_obj_range; }
		std::vector<std::vector<Real>> getEE() { return m_explore_exploit_ratio; }
		Real getVarSpaceVolume() { return m_var_space_volume; }
		std::vector<std::vector<Real>> getVarSpaceVolumeRatio() { return m_subspace_volume_ratio; }
		std::vector<std::vector<size_t>> getPopResources() { return m_pop_resources; }
		size_t getFrontLastGens(){ return m_front_last_gens; }

		MO_HLC& getMO_HLC() const { return *m_mo_hlc; }
		Subspace& getExploitSubspace() const { return *m_exploit_box; }
		MultiPopulation<SPMOEA_pop>& getPop() { return *m_multi_pop; }

		std::vector <std::shared_ptr<Solution<>>>& getHisSols() { return m_historical_sols; }
		std::vector <std::shared_ptr<Solution<>>>& getHisFrontSols() { return m_history_front_sols; }

		std::vector<std::vector <std::shared_ptr<Solution<>>>>& getInteractiveSols() { return m_interactive_sol_pair; }
		std::vector<std::vector <std::shared_ptr<Solution<>>>>& getHisEvolveLocus() { return m_his_evolve_locus; }
		std::vector<std::vector<std::pair<std::shared_ptr<Solution<>>, std::shared_ptr<Solution<>>>>>& getGenEvolveLocus() { return m_gen_evolve_locus; }
		std::vector <std::shared_ptr<Solution<>>>& getArchiveSols() { return m_archive; }

		size_t getPopsize() { return m_pop_size; }
		void setPopSize(size_t v) { m_pop_size = v; }
		size_t archiveNum() { return m_archive_num; }
		void setArchiveNum(size_t v) { m_archive_num = v; }
		void updateEE(size_t explore_num, size_t exploit_num);
		void updateSpaceRatio(const std::vector<Real> &space_ratios){ m_subspace_volume_ratio.emplace_back(space_ratios); }
		void updatePopResource(const std::vector<size_t>& pop_resource) { m_pop_resources.emplace_back(pop_resource); }
		void setFrontSpaceRatio(Real space_ratio) { m_front_space_ratio= space_ratio; }
		Real getFrontSpaceRatio() { return m_front_space_ratio; }
		size_t getAddNumSpace() { return m_add_num_spaces; }
		std::vector<std::vector<Real>>& getOperatorRatio() { return m_operators_ratio; }
		std::vector<std::vector<size_t>>& getIndSpaces() { return m_ind_spaces; }
		std::vector<size_t>& getParentSpaces() { return m_parent_spaces; }
		std::vector<size_t>& getOffspringSpaces() { return m_offspring_spaces; }

		Real getCr() { return m_cr; }
		Real getMr() { return m_mr; }
		Real getCeta() { return m_ceta; }
		Real getMeta() { return m_meta; }
		std::vector<Real> getIGD() { return m_IGD; }
		size_t getSplitGranularity() { return m_split_granularity; }
		std::vector<size_t> getFrontSpace() { return m_front_space; }

		std::vector<Real> getMaxRefPoint() {return m_max_ref_point; }
		bool ifPlotEvolveLocus() { return m_add_evolve_locus; }

		void setCr(Real v) { m_cr = v; }
		void setMr(Real v) { m_mr = v; }
		void setCeta(Real v) { m_ceta = v; }
		void setMeta(Real v) { m_meta = v; }
		void setIGD(std::vector<Real> &v) { m_IGD = v; }
		void setFrontLastGens(size_t v) { m_front_last_gens=v; }
		void setAddNumSpace(size_t v) { m_add_num_spaces = v; }

		void setMaxRefPoint(std::vector<Real>& p) { m_max_ref_point = p; }

#ifdef OFEC_DEMO
		void updateBuffer();
#endif
	protected:
		void run_() override;

	private:
		//std::unique_ptr<MOSP_pop> m_pop;//总种群
		std::unique_ptr<Subspace> m_exploit_box;
		std::unique_ptr<MultiPopulation<SPMOEA_pop>> m_multi_pop;//多种群
		std::unique_ptr<MO_HLC> m_mo_hlc;   //空间分割树
		size_t m_pop_size;
		
		Population<Solution<>> m_new_pops;
		size_t m_num_region_var = 0;
		size_t m_num_region_obj = 0;
		//size_t m_num_series = 10;//用于环境变化时间序列预测的长度
		std::vector<std::pair<Real, Real>> m_pop_obj_range;//当前种群目标值范围
		std::vector<std::pair<Real, Real>> m_pop_var_range;//当前种群搜索范围
		std::vector<std::pair<Real, Real>> m_front_obj_range;//前沿种群目标范围
		std::vector<std::pair<Real, Real>> m_his_obj_range;//历史解目标范围
		
		std::vector <std::shared_ptr<Solution<>>> m_history_front_sols;//store all non-diminated Solutions
		std::vector <std::shared_ptr<Solution<>>> m_cur_front_sols;//store current front sols
		std::vector <std::shared_ptr<Solution<>>> m_archive;//store a number of front sols
		std::vector <std::shared_ptr<Population<Solution<>>>> m_gen_front_sols;//store front sols in each generation
		std::vector <std::shared_ptr<Solution<>>> m_historical_sols;//store historical sols

		std::vector<Solution<>> m_subobj_opt_sol;//子目标最优的边界个体
		
		size_t m_converge_age = 10;//设置收敛与否的阈值
		Real m_R1, m_R2, m_R3;//三个衡量当前种群在目标空间的比值
		std::vector<Real> m_IGD;
		std::vector<size_t> m_over_space_index;//记录选择个体数不足的子空间
		std::vector<std::vector<size_t>> m_front_space_inx;//前沿子空间索引

		size_t m_archive_num;   //存档的规模
		bool m_normalize = true;//目标空间归一化标记

		Real m_cr, m_mr, m_ceta, m_meta;

		Real m_var_space_volume;
		std::vector<std::vector<Real>> m_explore_exploit_ratio;
		std::vector<std::vector<Real>> m_subspace_volume_ratio;
		Real m_front_space_ratio;
		std::vector<std::vector<size_t>> m_pop_resources;

		size_t m_max_split_times = 3;

		size_t m_subtree_num = 5;

		size_t m_init_pop_size = 100;
		size_t m_sub_pop_size = 20;
		std::vector<size_t> m_visit_count;
		//子空间单调性打标

		std::vector<std::vector<std::vector<size_t>>> m_region_link_space;
		std::vector<std::vector<size_t>> m_front_region_link_space;

		size_t m_split_granularity;
		std::vector<size_t> m_front_space;
		size_t m_add_num_spaces;
		size_t m_front_last_gens;
		Real m_niche_radius=0.;

		std::vector<std::vector<Real>> m_operators_ratio;

		//记录交互解
		std::vector<std::vector<std::shared_ptr<Solution<>>>> m_interactive_sol_pair;//2个或3个
		//记录个体历史演化轨迹
		std::vector<std::vector<std::shared_ptr<Solution<>>>> m_his_evolve_locus;
		//记录代际演化轨迹
		std::vector<std::vector<std::pair<std::shared_ptr<Solution<>>, std::shared_ptr<Solution<>>>>> m_gen_evolve_locus;
		//记录父代子代所在子空间信息
		std::vector<std::vector<size_t>> m_ind_spaces;
		std::vector<size_t> m_parent_spaces;
		std::vector<size_t> m_offspring_spaces;

		std::vector<Real> m_max_ref_point;
		bool m_add_evolve_locus = true;

		size_t m_stage_last_time = 0;
	};

}

#endif  // !OFEC_SPMOEA_H
