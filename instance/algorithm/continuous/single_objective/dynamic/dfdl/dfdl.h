/********* Begin Register Information **********
{
	"name": "DFDL",
	"identifier": "DFDL",
	"problem tags": [ "ConOP", "SOP", "GOP", "MMOP", "DOP" ]
}
*********** End Register Information **********/

//
// Created by Mai on 2022/11/16.
//
#ifndef OFEC_DFDL_H
#define OFEC_DFDL_H
#include<cmath>
#include "../metrics_dynamic.h"
#include "../../../../template/framework/hghe/continuous/hghec.h"
#include "../../../../template/classic/de/population.h"
#include "../../../../../../core/algorithm/multi_population.h"
#include "../../../../../../utility/clustering/hslh.h"
#include "../../../../../../utility/clustering/nbc.h"

#include "../../../../template/classic/pso/particle.h"
#include "../../../../template/classic/pso/swarm.h"
#include "../../multi_modal/lips/lips_pop.h"

#include "dfdl_particle.h"
#include "dfdl_swarm.h"


namespace ofec {
    class DFDL : public HGHEC, public MetricsDynamicConOEA {
        enum class PopStart { Niching_Lips, MultiPopulation };
        enum class ClusteringType { HSLH, NBC };

        using TypePop = DFDLSwarm;
        using Solution_type = DFDLParticle;
    private:
        void calculateUncertaintySSPbyLinear();
        void calculateUncertaintySSPbyExponential();
    protected:
        int m_iteration = 0;
        int m_init_num_ssp;
        int m_min_init_pop_size = 5;
        size_t m_outdated_evaluation = 2500;
        PopStart m_pop_start_type;
        ClusteringType m_clustering_type;
        int m_best_ssp_index = -1, m_worst_ssp_index = -1;
        Real m_alpha, m_beta, m_gama;
        int m_obsoleteFactor;
        int m_uncertainty_type;
        enum UncertaintyType{exponential = 0, linear};
        size_t m_size_pop;
        size_t m_num_increase;
        MultiPopulation<TypePop> m_sub_pop;                  // Multi-Population for DFDL
        std::unique_ptr<nanoflann::PartitioningKDTree<Real>> ssp_tree_for_coverage;    //Calculate Coverage
        std::vector<Real> m_probability_SSP;
        std::vector<std::unique_ptr<Solution_type>> m_new_indis;
        std::vector<std::unique_ptr<Solution<>>> m_archive_best_indis;
        std::vector<std::unique_ptr<Solution_type>> m_rest_indis;

        std::vector<std::unique_ptr<Solution<>>> m_archive_best;
        const Real converge_radius = 10e-3;
        const Real converge_vel = 10e-3;
        const Real seed_radius = 10e-1;
        const int min_subpopSize = 1;                        // minimum cluster size
        const int m_min_num_inds = 3;                        // minimum subpop size
        Real m_min_active_hill_volume;

        // for divided pop
        int m_new_pop_evolve_count = 0;
        // for Pop State
        size_t m_origin_cover_rate;   // explore -> exploit
        Real m_divide_factor = 0.5;   // explore -> exploit
        Real m_R;   // explore -> exploit
        Real m_r;   // exploit -> tracker

        // History Memory
        bool m_first_update = true;
        std::vector<std::vector<double>> m_pop_indis;
        size_t m_pre_hill_size = 0;
        size_t m_cur_hill_size = 0;

        size_t m_next_indi_size = 0;
        size_t m_pre_indi_size = 0;
        size_t m_total_indi_size = 0;

        const size_t mc_off_peak = 3;
        const size_t mc_step_indis = 5;

        const size_t mc_max_idxs = 1000;
        const size_t mc_min_idxs = 10;
    protected:
        //Unique Hill
        std::vector<std::vector<size_t>> m_unique_ssp_hills;
        std::vector<std::vector<size_t>> m_hills_pops;

        void updateMemory();
        void updateNumSolutions();
        std::vector<std::vector<size_t>> getUniqueSSPinHills(const std::vector<std::vector<size_t>> & arrays);
        void updateSwarmHillLocated();
        void updateUniqueHillsSwarm();
        bool isAllow2SupplySwarm(int hill);
        bool isAnySolutionInHill(int hill);
        int SolutionLocatedInUniqueHills(const Solution_type & Solution);

            // for Niching method
    protected:
        std::unique_ptr<SwarmLIP> m_swarm_Lips;
        bool detectPopDivision(std::unique_ptr<SwarmLIP> & niche, int type);
    protected:
        void initialize_() override;
        void initSolutions(int Solution_num, int type);
        int createNewSwarms(int Solution_num, int type);
        int createNichingSwarm(int Solution_num, int type);
        int popEvolve();
        int checkOverlapping();
        void checkOverCrowdAndConverge();
        int checkStagnantion();
        void avoidRepeatSearch();
        void hillShielding();

        void run_() override;
        void handleChange(int rf);
        void measureMultiPops(bool flag);
        void updateInfoSSP(const SolutionType *sol, TaskEval task) override;
        void removeOutdatedHisSolutions(size_t outdated_evaluation);

    public:
        enum InitSolType {
            Random = 0, AtHighUncertainty
        };
        int getNumberAboveMeanUncertainty();
        void calculateProbabilitySSP();
        void distinctHisSols();
        void getBestAndWorstSSPIndex();
        bool checkDiversity(Problem *pro);
        void awakeSwarms();
        void updateHills() override;
        int classifyPopulation() override { return 0; };
        int exploitHills() override { return 0; };
        void updatePotentialExploit() override {};
        void groupSubspace() override;
        int exploreHills() override { return 0; };

    public:
        void initSubSwarm();
        int createSupplementarySwarm(int pop_size);
        int createExploreSwarm();
        int createExploitSwarm(int pop_size);
        int createSubSwarm(int num);
        void record() override {};
#ifdef OFEC_DEMO
        // drawSolutionMode = CurSolution | HisSolution
        std::string drawSolutionMode = "CurSolution";
    public:
        virtual void updateBuffer();
#endif

    private:
//        int m_last_ner, m_last_nei;
//        void updateConvergedPoints();
        void archiveEvolveSolution();
        void updateUncertaintySSP();
    };
}


#endif //OFEC_DFDL_H
