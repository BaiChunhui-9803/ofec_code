#ifndef OFEC_NL_SHADE_LBC_OPTIMIZAR_H
#define OFEC_NL_SHADE_LBC_OPTIMIZAR_H


#include <vector>
#include "../../../../../../core/environment/environment.h"
#include "../../../../../../core/random/newran.h"

namespace ofec {

    namespace nl_shade_lbc {

        class Optimizer
        {
        public:
            bool FitNotCalculated;
            int Int_ArchiveSizeParam;
            int MemorySize;
            int MemoryIter;
            int SuccessFilled;
            int MemoryCurrentIndex;
            int NVars;
            int NInds;
            int NIndsMax;
            int NIndsMin;
            int besti;
            int Generation;
            int ArchiveSize;
            int CurrentArchiveSize;
            double MWLp1;
            double MWLp2;
            double MWLm;
            double LBC_fin;
            double F;
            double Cr;
            double bestfit;
            double ArchiveSizeParam;
            double Right;
            double Left;
            std::vector<int> Rands;
            std::vector<int> Indexes;
            std::vector<int> BackIndexes;
            std::vector<double> Weights;
            std::vector<double> Donor;
            std::vector<double> Trial;
            std::vector<double> FitMass;
            std::vector<double> FitMassTemp;
            std::vector<double> FitMassCopy;
            std::vector<double> BestInd;
            std::vector<double> tempSuccessCr;
            std::vector<double> tempSuccessF;
            std::vector<double> FGenerated;
            std::vector<double> CrGenerated;
            std::vector<double> MemoryCr;
            std::vector<double> MemoryF;
            std::vector<double> FitDelta;
            std::vector<double> FitMassArch;
            std::vector<std::vector<double>> Popul;
            std::vector<std::vector<double>> PopulTemp;
            std::vector<std::vector<double>> Archive;

            std::function<double(std::vector<double>& x, ofec::Environment* env)> m_eval_fun;

            void Initialize(int newNInds, int newNVars, int NewMemSize, double NewArchSizeParam, ofec::Random* rnd);
            void Clean();
            void MainCycle(ofec::Environment* env, ofec::Random* rnd);
            void FindNSaveBest(bool init, int ChosenOne);
            inline double GetValue(const int index, const int NInds, const int j);
            void CopyToArchive(double* RefusedParent, double RefusedFitness, ofec::Random* rnd);
            void SaveSuccessCrF(double Cr, double F, double FitD);
            void UpdateMemoryCrF(ofec::Environment* env);
            double MeanWL_general(double* Vector, double* TempWeights, int Size, double g_p, double g_m);
            void RemoveWorst(int NInds, int NewNInds);
        };
    }

}
#endif