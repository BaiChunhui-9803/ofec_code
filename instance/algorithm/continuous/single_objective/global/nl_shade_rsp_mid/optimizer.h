// This is the code of nl_shade_rsp_mid.
// Most of the code is taken from nl_shade_rsp, which was authored by:
// Vladimir Stanovov, Shakhnaz  Akhmedova, Eugene Semenkin 
// and downloaded from: https://github.com/P-N-Suganthan
// Changes were conceptually described in the article: "A Version of NL-SHADE-RSP Algorithm with Midpoint for CEC 2022 Single Objective Bound Constrained Problems"
// In the code, they can be found in conditional compilation blocks enabled by:
// #define K_MEANS_AS_NEAREST
// #define RESAMPLING
// #define COUNT_LIMITS
// Author of the changes: Rafal Biedrzycki

// To compile: (tested under Ubuntu 22.04)
// Use system package manager to install libraries: mlpack and boost (libmlpack-dev, libboost-all-dev)
// compile & link by command: g++ -std=c++17 nl_shade_rsp_mid.cpp -lmlpack -fopenmp

// Before running the program, run command (or-else program will use many threads wchich slowes down in my case):
// export OMP_NUM_THREADS=1


#ifndef OFEC_NL_SHADE_RSP_MID_OPTIMIZAR_H
#define OFEC_NL_SHADE_RSP_MID_OPTIMIZAR_H



#include <vector>
#include "../../../../../../core/environment/environment.h"
#include "../../../../../../core/random/newran.h"



namespace ofec {

    
#define K_MEANS_AS_NEAREST
#define RESAMPLING
#define COUNT_LIMITS
    // only for minization 
    class Optimizer
    {

    protected:


        // ≈∑ œæ‡¿Îº∆À„∫Ø ˝
        struct EuclideanDistance {
            double operator()(const std::vector<double>& a, const std::vector<double>& b) const {
                double sum = 0.0;
                for (size_t i = 0; i < a.size(); ++i) {
                    double diff = a[i] - b[i];
                    sum += diff * diff;
                }
                return std::sqrt(sum);
            }
        };

    public:
        bool FitNotCalculated;
        int Int_ArchiveSizeParam;
        int MemorySize;
        int MemoryIter;
        int SuccessFilled;
        int MemoryCurrentIndex;
        int NVars;
        int popSize;
        int NIndsMax;
        int NIndsMin;
        int besti;


        int Generation;
        int ArchiveSize;
        int CurrentArchiveSize;
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
        std::vector<double> Fitmass;
        std::vector<double> popFitTmp;
        std::vector<double> FitmassCopy;
        std::vector<double> BestInd;
        std::vector<double> tempSuccessCr;
        std::vector<double> tempSuccessF;
        std::vector<double> FGenerated;
        std::vector<double> CrGenerated;
        std::vector<double> MemoryCr;
        std::vector<double> MemoryF;
        std::vector<double> FitDelta;
        std::vector<double> ArchUsages;
        std::vector<std::vector<double>> Popul;
        std::vector<std::vector<double>> populTemp;
        std::vector<std::vector<double>> Archive;
        //  ofstream infoLog;


        bool isRestart = false;
        int evalsAtStart = 0;

        std::function<double(std::vector<double>& x, ofec::Environment* env)> m_eval_fun;


        void Initialize(int newNInds, int newNVars, int NewMemSize, double NewArchSizeParam, ofec::Random* rnd);
        void restart(int newNInds, int newNVars, int NewMemSize, double NewArchSizeParam, ofec::Environment* env, ofec::Random* rnd);
        void Clean();
        void MainCycle(double a, ofec::Environment* env, ofec::Random* rnd);
        void FindNSaveBest(bool init, int ChosenOne);
        inline double GetValue(const int index, const int popSize, const int j);
        void CopyToArchive(double* RefusedParent, ofec::Random* rnd);
        void SaveSuccessCrF(double Cr, double F, double FitD);
        void UpdateMemoryCrF();
        double MeanWL_general(double* Vector, double* TempWeights, int Size, double g_p, double g_m);
        void RemoveWorst(int popSize, int NewNInds);
        void RemoveTooNear(int popSize, int NewNInds, ofec::Environment* env);
    };
}

#endif

