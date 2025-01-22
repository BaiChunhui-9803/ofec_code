#include "optimizer.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <random>
//#include <boost/math/statistics/linear_regression.hpp>
#include <cstring>

#include "../../../../../../core/random/newran.h"
#include "../../../../../../core/problem/continuous/continuous.h"
#include "../../../../../../utility/clustering/kmeans_armilo.h"
#include "../../../../../../utility/metric/silhouette_score.h"


#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS

#include "../../../../../../datum/datum_inclusion.h"



namespace ofec {
	
    int intRandom(int target, ofec::Random* rnd)
    {
        if (target == 0)
            return 0;
        return rnd->uniform.nextNonStd<int>(0, target);
        //return uni_int(generator) % target;
    }
    double RandomD(double minimal, double maximal, ofec::Random* rnd) {
        return rnd->uniform.nextNonStd<double>(minimal, maximal);
        //return uni_real(generator) * (maximal - minimal) + minimal; 
    }
    double NormRand(double mu, double sigma, ofec::Random* rnd) {
        return rnd->normal.next() * sigma + mu;
        //return norm_dist(generator) * sigma + mu; 
    }
    double CachyRand(double mu, double sigma, ofec::Random* rnd) {
        return rnd->cauchy.next() * sigma + mu;
        //return cachy_dist(generator) * sigma + mu; 
    }

    void seqDouble(std::vector<double>& seq, int size) {
        for (int i = 1; i <= size; ++i) {
            seq[i - 1] = i;
        }
    }

    void getMean(const std::vector<std::vector<double>>& popul, int popul_size, int dimensionality, std::vector<double>& mean_indiv) {
        for (int dim = 0; dim < dimensionality; ++dim) {
            double sum = 0;
            for (int indIndx = 0; indIndx < popul_size; ++indIndx) {
                sum += popul[indIndx][dim];
            }
            mean_indiv[dim] = sum / popul_size;
        }
    }



    void getMean(const std::vector<std::vector<double>>& popul, int popul_size, int dimensionality, std::vector<double>& mean_indiv, const std::vector<unsigned>& ids) {
        int idsSize = ids.size();
        for (int dim = 0; dim < dimensionality; ++dim) {
            double sum = 0;
            for (int idsIndx = 0; idsIndx < idsSize; ++idsIndx) {
                int indIndx = ids[idsIndx];
                sum += popul[indIndx][dim];
            }
            mean_indiv[dim] = sum / idsSize;
        }
    }

    //void getMean(double** popul, int popul_size, int dimensionality, double* mean_indiv, arma::uvec ids) {
    //    int idsSize = ids.size();
    //    for (int dim = 0; dim < dimensionality; ++dim) {
    //        double sum = 0;
    //        for (int idsIndx = 0; idsIndx < idsSize; ++idsIndx) {
    //            int indIndx = ids[idsIndx];
    //            sum += popul[indIndx][dim];
    //        }
    //        mean_indiv[dim] = sum / idsSize;
    //    }
    //}

    int getIndxOfNearest(const std::vector<std::vector<double>>& popul, int popul_size, int dimensionality, const std::vector<double>& mean_indiv) {
        int nearestIndx = 0;
        double minDist = 1e20;
        for (int indIndx = 0; indIndx < popul_size; ++indIndx) {
            double sqDist = 0;
            for (int dim = 0; dim < dimensionality; ++dim) {
                sqDist += pow(mean_indiv[dim] - popul[indIndx][dim], 2);
            }
            if (sqDist < minDist) {
                minDist = sqDist;
                nearestIndx = indIndx;
            }
        }
        return nearestIndx;
    }

    double getDist(double* ind1, double* ind2, int dimensionality) {
        double dist = 0.0;
        for (int dim = 0; dim < dimensionality; ++dim) {
            dist += pow(ind1[dim] - ind2[dim], 2);
        }
        return sqrt(dist);
    }

    void qSort1(double* mass, int low, int high)
    {
        int i = low;
        int j = high;
        double x = mass[(low + high) >> 1];
        do
        {
            while (mass[i] < x)    ++i;
            while (mass[j] > x)    --j;
            if (i <= j)
            {
                double temp = mass[i];
                mass[i] = mass[j];
                mass[j] = temp;
                i++;    j--;
            }
        } while (i <= j);
        if (low < j)   qSort1(mass, low, j);
        if (i < high)  qSort1(mass, i, high);
    }

    void qSort2int(double* mass, int* mass2, int low, int high)
    {
        int i = low;
        int j = high;
        double x = mass[(low + high) >> 1];
        do
        {
            while (mass[i] < x)    ++i;
            while (mass[j] > x)    --j;
            if (i <= j)
            {
                double temp = mass[i];
                mass[i] = mass[j];
                mass[j] = temp;
                int temp2 = mass2[i];
                mass2[i] = mass2[j];
                mass2[j] = temp2;
                i++;    j--;
            }
        } while (i <= j);
        if (low < j)
            qSort2int(mass, mass2, low, j);
        if (i < high)
            qSort2int(mass, mass2, i, high);
    }



    void GenerateNextRandUnif(const int num, const int Range, int* Rands, const int Prohib, ofec::Random* rnd)
    {
        for (int j = 0; j != 25; j++)
        {
            bool generateagain = false;
            Rands[num] = intRandom(Range, rnd);
            for (int i = 0; i != num; i++)
                if (Rands[i] == Rands[num])
                    generateagain = true;
            if (!generateagain)
                break;
        }
    }

    void GenerateNextRandUnifOnlyArch(const int num, const int Range, const int Range2, int* Rands, const int Prohib, ofec::Random* rnd)
    {
        for (int j = 0; j != 25; j++)
        {
            bool generateagain = false;
            Rands[num] = intRandom(Range2, rnd) + Range;
            for (int i = 0; i != num; i++)
                if (Rands[i] == Rands[num])
                    generateagain = true;
            if (!generateagain)
                break;
        }
    }

    bool CheckGenerated(const int num, int* Rands, const int Prohib)
    {
        if (Rands[num] == Prohib)
            return false;
        for (int j = 0; j != num; j++)
            if (Rands[j] == Rands[num])
                return false;
        return true;
    }

    //void SaveBestValues(int funcN, int runNum)
    //{
    //    for (int stepFEcount = LastFEcount; stepFEcount < NUM_OF_DUMPS; stepFEcount++)
    //    {
    //        if (NFEval == int(stepsFEval[stepFEcount] * maxFES))
    //        {
    //            double temp = globalBestFit - fopt;
    //            if (temp <= MIN_ERROR)
    //                temp = 0;
    //            ResultsArray[funcN - 1][runNum][stepFEcount] = temp;
    //            LastFEcount = stepFEcount;
    //        }
    //    }
    //
    //}

    void FindLimits(double* Ind, double* Parent, int CurNVars, double CurLeft, double CurRight, ofec::Random* rnd) {
        for (int j = 0; j < CurNVars; ++j) {
            if (Ind[j] < CurLeft)
                Ind[j] = RandomD(CurLeft, CurRight, rnd);
            if (Ind[j] > CurRight)
                Ind[j] = RandomD(CurLeft, CurRight, rnd);
        }
    }

    void countLimits(double* ind, int* indLimCounter, int dimensionality, double lower, double upper) {
        bool onBound = false;
        for (int dim = 0; dim < dimensionality; ++dim) {
            if (ind[dim] <= lower) {
                onBound = true;
                break;
            }
            if (ind[dim] >= upper) {
                onBound = true;
                break;
            }
        }
        if (onBound) {
            ++(*indLimCounter);
        }
        else {
            *indLimCounter = 0;
        }
    }

    //int getNumOfIndOnBounds(double** populTemp, int DIM, int popSize, double lower, double upper) {
    //    int numOfIndOnB = 0;
    //    for (int indIndx = 0; indIndx < popSize; ++indIndx) {
    //        for (int dim = 0; dim < dimensionality; ++dim) {
    //            if (populTemp[indIndx][dim] <= lower || populTemp[indIndx][dim] >= upper) {
    //                ++numOfIndOnB;
    //                break;
    //            }
    //        }
    //    }
    //    return numOfIndOnB;
    //}

    bool IsInfeasible(double* Ind, int CurNVars, double CurLeft, double CurRight) {
        for (int j = 0; j < CurNVars; ++j) {
            if (Ind[j] < CurLeft)
                return true;
            if (Ind[j] > CurRight)
                return true;
        }
        return false;
    }



    void Optimizer::Initialize(int newNInds, int newNVars,
        int NewMemSize, double NewArchSizeParam, ofec::Random* rnd)
    {

        m_eval_fun = [](std::vector<double>& x, ofec::Environment* env) {
            ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
            ofec::Continuous::SolutionType sol(con_pro->numberObjectives(), con_pro->numberConstraints(), con_pro->numberVariables());
            sol.variable() = x;
            sol.evaluate(env);
            ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
            sol.setFitness(pos * sol.objective(0));
            return -sol.fitness();
            };

        isRestart = false;
        evalsAtStart = 0;
        FitNotCalculated = true;
        popSize = newNInds;
        NIndsMax = popSize;
        NIndsMin = 4;
        NVars = newNVars;
        Left = -100;//+10000;
        Right = 100;//+10000;
        Cr = 0.2;
        F = 0.2;
        besti = 0;
        Generation = 0;
        CurrentArchiveSize = 0;
        ArchiveSizeParam = NewArchSizeParam;
        Int_ArchiveSizeParam = ceil(ArchiveSizeParam);
        ArchiveSize = NIndsMax * ArchiveSizeParam;

        //for (int steps_k = 0; steps_k != NUM_OF_DUMPS - 1; steps_k++)
        //    stepsFEval[steps_k] = pow(double(dimensionality), double(steps_k) / 5.0 - 3.0);
        //stepsFEval[NUM_OF_DUMPS - 1] = 1.0;

        Popul.resize(NIndsMax);
        for (auto& it : Popul) {
            it.resize(NVars);
        }
        populTemp = Popul;

        Archive.resize(NIndsMax * Int_ArchiveSizeParam);
        for (auto& it : Archive) {
            it.resize(NVars);
        }
        Fitmass.resize(NIndsMax);
        popFitTmp.resize(NIndsMax);

        FitmassCopy.resize(NIndsMax);
        Indexes.resize(NIndsMax);
        BackIndexes.resize(NIndsMax);

        BestInd.resize(NVars);



        for (int i = 0; i < NIndsMax; i++)
            for (int j = 0; j < NVars; j++)
                Popul[i][j] = RandomD(Left, Right, rnd);
        Donor.resize(NVars);
        Trial.resize(NVars);
        Rands.resize(NIndsMax);
        tempSuccessCr.resize(NIndsMax);
        tempSuccessF.resize(NIndsMax);
        FitDelta.resize(NIndsMax);
        FGenerated.resize(NIndsMax);
        CrGenerated.resize(NIndsMax);
        for (int i = 0; i != NIndsMax; i++)
        {
            tempSuccessCr[i] = 0;
            tempSuccessF[i] = 0;
        }
        MemorySize = NewMemSize;
        MemoryIter = 0;
        SuccessFilled = 0;
        ArchUsages.resize(NIndsMax);
        Weights.resize(NIndsMax);
        MemoryCr.resize(MemorySize);
        MemoryF.resize(MemorySize);
        for (int i = 0; i != MemorySize; i++)
        {
            MemoryCr[i] = 0.2;
            MemoryF[i] = 0.2;
        }
        // foundAt = maxFES;

    }

    void Optimizer::restart(int newNInds, int newNVars,
        int NewMemSize, double NewArchSizeParam, ofec::Environment* env, ofec::Random* rnd)
    {
        auto alg = env->algorithm();
        isRestart = true;
        evalsAtStart = alg->evaluations();
        FitNotCalculated = true;
        popSize = newNInds;
        NIndsMax = popSize;
        NIndsMin = 4;
        NVars = newNVars;
        Cr = 0.2;
        F = 0.2;
        besti = 0;
        Generation = 0;
        CurrentArchiveSize = 0;
        ArchiveSizeParam = NewArchSizeParam;
        Int_ArchiveSizeParam = ceil(ArchiveSizeParam);
        ArchiveSize = NIndsMax * ArchiveSizeParam;

        Popul.resize(NIndsMax);
        for (auto& it : Popul) {
            it.resize(NVars);
        }
        populTemp = Popul;

        Archive.resize(NIndsMax * Int_ArchiveSizeParam);
        for (auto& it : Archive) {
            it.resize(NVars);
        }
        Fitmass.resize(NIndsMax);
        popFitTmp.resize(NIndsMax);
        FitmassCopy.resize(NIndsMax);
        Indexes.resize(NIndsMax);
        BackIndexes.resize(NIndsMax);
        BestInd.resize(NVars);


        for (int i = 0; i < NIndsMax; i++)
            for (int j = 0; j < NVars; j++)
                Popul[i][j] = RandomD(Left, Right, rnd);



        Rands.resize(NIndsMax);
        tempSuccessCr.resize(NIndsMax);
        tempSuccessF.resize(NIndsMax);
        FitDelta.resize(NIndsMax);
        FGenerated.resize(NIndsMax);
        CrGenerated.resize(NIndsMax);
        for (int i = 0; i != NIndsMax; i++)
        {
            tempSuccessCr[i] = 0;
            tempSuccessF[i] = 0;
        }
        MemorySize = NewMemSize;
        MemoryIter = 0;
        SuccessFilled = 0;
        ArchUsages.resize(NIndsMax);
        Weights.resize(NIndsMax);
        MemoryCr.resize(MemorySize);
        MemoryF.resize(MemorySize);
        for (int i = 0; i != MemorySize; i++)
        {
            MemoryCr[i] = 0.2;
            MemoryF[i] = 0.2;
        }


    }

    void Optimizer::SaveSuccessCrF(double Cr, double F, double FitD)
    {
        tempSuccessCr[SuccessFilled] = Cr;
        tempSuccessF[SuccessFilled] = F;
        FitDelta[SuccessFilled] = FitD;
        SuccessFilled++;
    }
    void Optimizer::UpdateMemoryCrF()
    {
        if (SuccessFilled != 0)
        {
            MemoryCr[MemoryIter] = MeanWL_general(tempSuccessCr.data(), FitDelta.data(), SuccessFilled, 2, 1);
            MemoryF[MemoryIter] = MeanWL_general(tempSuccessF.data(), FitDelta.data(), SuccessFilled, 2, 1);
            MemoryIter++;
            if (MemoryIter >= MemorySize)
                MemoryIter = 0;
        }
        else
        {
            MemoryF[MemoryIter] = 0.5;
            MemoryCr[MemoryIter] = 0.5;
        }
    }
    double Optimizer::MeanWL_general(double* Vector, double* TempWeights, int Size, double g_p, double g_m)
    {
        double SumWeight = 0;
        double SumSquare = 0;
        double Sum = 0;
        for (int i = 0; i != SuccessFilled; i++)
            SumWeight += TempWeights[i];
        for (int i = 0; i != SuccessFilled; i++)
            Weights[i] = TempWeights[i] / SumWeight;
        for (int i = 0; i != SuccessFilled; i++)
            SumSquare += Weights[i] * pow(Vector[i], g_p);
        for (int i = 0; i != SuccessFilled; i++)
            Sum += Weights[i] * pow(Vector[i], g_p - g_m);
        if (fabs(Sum) > 0.000001)
            return SumSquare / Sum;
        else
            return 0.5;
    }

    void Optimizer::CopyToArchive(double* RefusedParent, ofec::Random* rnd)
    {
        if (CurrentArchiveSize < ArchiveSize)
        {
            for (int i = 0; i != NVars; i++)
                Archive[CurrentArchiveSize][i] = RefusedParent[i];
            CurrentArchiveSize++;
        }
        else if (ArchiveSize > 0)
        {
            int RandomNum = intRandom(ArchiveSize, rnd);
            for (int i = 0; i != NVars; i++)
                Archive[RandomNum][i] = RefusedParent[i];
        }
    }

    void Optimizer::FindNSaveBest(bool init, int ChosenOne)
    {
        if (Fitmass[ChosenOne] <= bestfit || init)
        {
            bestfit = Fitmass[ChosenOne];
            besti = ChosenOne;
            for (int j = 0; j != NVars; j++)
                BestInd[j] = Popul[besti][j];
        }
        //if (bestfit < globalBestFit) {
        //    globalBestFit = bestfit;
        //    memcpy(bestSol, Popul[ChosenOne], dimensionality * sizeof(double));
        //}
    }

    void Optimizer::RemoveTooNear(int popSize, int NewNInds, ofec::Environment* env) {
        ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
        int dimensionality = con_pro->numberVariables();

        int PointsToRemove = popSize - NewNInds;
        int curPopSize = popSize;

        double TOO_SMALL_DIST = 1e-7;


        for (int L = 0; L != PointsToRemove; L++) {
            double smallestDist = 1e20;
            int best_i = 0;
            int best_j = 0;
            for (int i = 0; i < curPopSize - 1; ++i) {
                for (int j = i + 1; j < curPopSize; ++j) {
                    double dist = getDist(Popul[i].data(), Popul[j].data(), dimensionality);
                    if (dist < smallestDist) {
                        smallestDist = dist;
                        best_i = i;
                        best_j = j;
                    }

                }
            }
            if (smallestDist <= TOO_SMALL_DIST) {
                for (int j = 0; j != NVars; ++j)
                    Popul[best_i][j] = Popul[curPopSize - 1][j];
                Fitmass[best_i] = Fitmass[curPopSize - 1];
                --curPopSize;
            }

        }
        int pointsLeftToRemove = curPopSize - NewNInds;
        if (pointsLeftToRemove > 0) {
            RemoveWorst(curPopSize, NewNInds);
        }
    }

    void Optimizer::RemoveWorst(int popSize, int NewNInds)
    {
        int PointsToRemove = popSize - NewNInds;
        for (int L = 0; L != PointsToRemove; L++)
        {
            double WorstFit = Fitmass[0];
            int WorstNum = 0;
            for (int i = 1; i != popSize; i++)
            {
                if (Fitmass[i] > WorstFit)
                {
                    WorstFit = Fitmass[i];
                    WorstNum = i;
                }
            }
            for (int i = WorstNum; i != popSize - 1; i++)
            {
                for (int j = 0; j != NVars; j++)
                    Popul[i][j] = Popul[i + 1][j];
                Fitmass[i] = Fitmass[i + 1];
            }
        }
    }

    inline double Optimizer::GetValue(const int index, const int popSize, const int j)
    {
        if (index < popSize)
            return Popul[index][j];
        return Archive[index - popSize][j];
    }


    void Optimizer::MainCycle(double optFit, ofec::Environment* env, ofec::Random* rnd)
    {
        using namespace std;

        auto alg = env->algorithm();
        ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
        int dimensionality = con_pro->numberVariables();

        double ArchSuccess;
        double NoArchSuccess;
        double NArchUsages;
        double ArchProbs = 0.5;

        std::vector<std::vector<double>> centroid_matrix (popSize);
        for (int i = 0; i < popSize; ++i) {
            centroid_matrix[i].resize(dimensionality);
        }

        vector<double> lam_seq(popSize);
        seqDouble(lam_seq, popSize);
        vector<tuple<double, double, double>> abTab(dimensionality);
        vector<double> centroids4dim(popSize);
        std::vector<double> est_m(dimensionality, 0);
        std::vector<double> mean_indiv(dimensionality, 0);
        std::vector<double> mean_indivCl0(dimensionality, 0);
        std::vector<double> mean_indivCl1(dimensionality, 0);
        std::vector<double>mean_indiv_old(dimensionality, 0);
        std::vector<int> populLimCount(popSize, 0);
        for (int indIndx = 0; indIndx < popSize; ++indIndx) {
            populLimCount[indIndx] = 0;
        }
        int numOfStagIt = 0;

        const int ITS_MODULO = 17;

        double fit_mean_old = 1e20;


        for (int curIndx = 0; curIndx != popSize; curIndx++)
        {


            Fitmass[curIndx] = m_eval_fun(Popul[curIndx], env);

            FindNSaveBest(curIndx == 0, curIndx);
            //if (!globalBestFitinit || bestfit < globalBestFit)
            //{
            //    globalBestFit = bestfit;
            //    globalBestFitinit = true;
            //    memcpy(bestSol, Popul[curIndx], dimensionality * sizeof(double));
            //}
            //SaveBestValues(func_num, runNum);
        }

#ifdef OFEC_DATUM_MULTI_POP_H
        {

            std::vector<ofec::Continuous::SolutionType> pop(popSize,
                { con_pro->numberObjectives(), con_pro->numberConstraints(),con_pro->numberVariables() });
            for (int idx(0); idx < pop.size(); ++idx) {
                pop[idx].variable().vector() = Popul[idx];
            }

            g_multi_pop.bindPopulation(pop);
            env->algorithm()->datumUpdated(env, g_multi_pop);
        }
#endif
        
        std::vector<double> FitTemp3;
        std::mt19937 generator(1e9*rnd->uniform.next());
        do
        {
            double minfit = Fitmass[0];
            double maxfit = Fitmass[0];
            for (int i = 0; i != popSize; i++)
            {
                FitmassCopy[i] = Fitmass[i];
                Indexes[i] = i;
                if (Fitmass[i] >= maxfit)
                    maxfit = Fitmass[i];
                if (Fitmass[i] <= minfit)
                    minfit = Fitmass[i];
            }
            if (minfit != maxfit)
                qSort2int(FitmassCopy.data(), Indexes.data(), 0, popSize - 1);
            for (int i = 0; i != popSize; i++)
                for (int j = 0; j != popSize; j++)
                    if (i == Indexes[j])
                    {
                        BackIndexes[i] = j;
                        break;
                    }
            FitTemp3.resize(popSize);
            for (int i = 0; i != popSize; i++)
                FitTemp3[i] = exp(-double(i) / (double)popSize);
            std::discrete_distribution<int> ComponentSelector3(FitTemp3.begin(), FitTemp3.end());

            int psizeval = max(2.0, popSize * (0.2 / (double)(alg->maximumEvaluations()) * (double)(alg->evaluations()) + 0.2));

            int CrossExponential = 0;
            if (RandomD(0, 1, rnd) < 0.5)
                CrossExponential = 1;
            for (int curIndx = 0; curIndx != popSize; curIndx++)
            {
                MemoryCurrentIndex = intRandom(MemorySize, rnd);
                Cr = min(1.0, max(0.0, NormRand(MemoryCr[MemoryCurrentIndex], 0.1, rnd)));
                do
                {
                    F = CachyRand(MemoryF[MemoryCurrentIndex], 0.1, rnd);
                } while (F <= 0);
                FGenerated[curIndx] = min(F, 1.0);
                CrGenerated[curIndx] = Cr;
            }
            qSort1(CrGenerated.data(), 0, popSize - 1);
            double iterBestFit = 1E18;
            int iterBestIndx = 0;

            for (int curIndx = 0; curIndx != popSize; curIndx++)
            {
                Rands[0] = Indexes[intRandom(psizeval, rnd)];
                for (int i = 0; i != 25 && !CheckGenerated(0, Rands.data(), curIndx); i++)
                    Rands[0] = Indexes[intRandom(psizeval, rnd)];
                GenerateNextRandUnif(1, popSize, Rands.data(), curIndx, rnd);
                if (RandomD(0, 1, rnd) > ArchProbs || CurrentArchiveSize == 0)
                {
                    Rands[2] = Indexes[ComponentSelector3(generator)];
                    for (int i = 0; i != 25 && !CheckGenerated(2, Rands.data(), curIndx); i++)
                        Rands[2] = Indexes[ComponentSelector3(generator)];
                    ArchUsages[curIndx] = 0;
                }
                else
                {
                    GenerateNextRandUnifOnlyArch(2, popSize, CurrentArchiveSize, Rands.data(), curIndx, rnd);
                    ArchUsages[curIndx] = 1;
                }
                for (int j = 0; j != NVars; j++)
                    Donor[j] = Popul[curIndx][j] +
                    FGenerated[curIndx] * (GetValue(Rands[0], popSize, j) - Popul[curIndx][j]) +
                    FGenerated[curIndx] * (GetValue(Rands[1], popSize, j) - GetValue(Rands[2], popSize, j));

                F = FGenerated[curIndx];//Corrected F
                int WillCrossover = intRandom(NVars, rnd);
                Cr = CrGenerated[BackIndexes[curIndx]];
                double CrToUse = 0;

                if (alg->evaluations() > 0.5 * (alg->maximumEvaluations()))
                    CrToUse = (double(alg->evaluations()) / double(alg->maximumEvaluations()) - 0.5) * 2;

                if (CrossExponential == 0)
                {
                    for (int j = 0; j != NVars; j++)
                    {
                        if (RandomD(0, 1, rnd) < CrToUse || WillCrossover == j)
                            populTemp[curIndx][j] = Donor[j];
                        else
                            populTemp[curIndx][j] = Popul[curIndx][j];
                    }
                }
                else
                {
                    int StartLoc = intRandom(NVars, rnd);
                    int L = StartLoc + 1;
                    while (RandomD(0, 1, rnd) < Cr && L < NVars)
                        L++;
                    for (int j = 0; j != NVars; j++)
                        populTemp[curIndx][j] = Popul[curIndx][j];
                    for (int j = StartLoc; j != L; j++)
                        populTemp[curIndx][j] = Donor[j];
                }
#ifdef RESAMPLING            
                const int MAX_NUM_OF_TRIALS = 100;
                int num_of_trials = 1;
                bool usedRepair = false;
                if (IsInfeasible(populTemp[curIndx].data(), NVars, Left, Right)) {
                    usedRepair = true;
                    const int NUM_OF_TR_WITHOUT_F_HANGE = 10;
                    do {
                        if (num_of_trials > NUM_OF_TR_WITHOUT_F_HANGE) {
                            CrossExponential = 0;
                            if (RandomD(0, 1, rnd) < 0.5)
                                CrossExponential = 1;


                            MemoryCurrentIndex = intRandom(MemorySize, rnd);
                            Cr = min(1.0, max(0.0, NormRand(MemoryCr[MemoryCurrentIndex], 0.1, rnd)));
                            do {
                                F = CachyRand(MemoryF[MemoryCurrentIndex], 0.1, rnd);
                            } while (F <= 0);
                            FGenerated[curIndx] = min(F, 1.0);
                            CrGenerated[curIndx] = Cr;
                        }


                        Rands[0] = Indexes[intRandom(psizeval, rnd)];
                        for (int i = 0; i != 25 && !CheckGenerated(0, Rands.data(), curIndx); i++)
                            Rands[0] = Indexes[intRandom(psizeval, rnd)];
                        GenerateNextRandUnif(1, popSize, Rands.data(), curIndx, rnd);
                        if (RandomD(0, 1, rnd) > ArchProbs || CurrentArchiveSize == 0)
                        {
                            Rands[2] = Indexes[ComponentSelector3(generator)];
                            for (int i = 0; i != 25 && !CheckGenerated(2, Rands.data(), curIndx); i++)
                                Rands[2] = Indexes[ComponentSelector3(generator)];
                            ArchUsages[curIndx] = 0;
                        }
                        else
                        {
                            GenerateNextRandUnifOnlyArch(2, popSize, CurrentArchiveSize, Rands.data(), curIndx, rnd);
                            ArchUsages[curIndx] = 1;
                        }
                        for (int j = 0; j != NVars; j++)
                            Donor[j] = Popul[curIndx][j] +
                            FGenerated[curIndx] * (GetValue(Rands[0], popSize, j) - Popul[curIndx][j]) +
                            FGenerated[curIndx] * (GetValue(Rands[1], popSize, j) - GetValue(Rands[2], popSize, j));
                        F = FGenerated[curIndx]; //corrected F
                        int WillCrossover = intRandom(NVars, rnd);
                        if (num_of_trials <= NUM_OF_TR_WITHOUT_F_HANGE) {
                            Cr = CrGenerated[BackIndexes[curIndx]];
                        }
                        double CrToUse = 0;
                        if (alg->evaluations() > 0.5 * (alg->maximumEvaluations()))
                            CrToUse = (double(alg->maximumEvaluations()) / double(alg->maximumEvaluations()) - 0.5) * 2;
                        if (CrossExponential == 0)
                        {
                            for (int j = 0; j != NVars; j++)
                            {
                                if (RandomD(0, 1, rnd) < CrToUse || WillCrossover == j)
                                    populTemp[curIndx][j] = Donor[j];
                                else
                                    populTemp[curIndx][j] = Popul[curIndx][j];
                            }
                        }
                        else
                        {
                            int StartLoc = intRandom(NVars, rnd);
                            int L = StartLoc + 1;
                            while (RandomD(0, 1, rnd) < Cr && L < NVars)
                                L++;
                            for (int j = 0; j != NVars; j++)
                                populTemp[curIndx][j] = Popul[curIndx][j];
                            for (int j = StartLoc; j != L; j++)
                                populTemp[curIndx][j] = Donor[j];
                        }
                        ++num_of_trials;
                    } while (IsInfeasible(populTemp[curIndx].data(), NVars, Left, Right) && num_of_trials <= MAX_NUM_OF_TRIALS);
                    if (IsInfeasible(populTemp[curIndx].data(), NVars, Left, Right)) {
                        usedRepair = false;
                        //cout << "still infeasible!" << endl;
                        FindLimits(populTemp[curIndx].data(), Popul[curIndx].data(), NVars, Left, Right, rnd);

                    }
                }


#ifdef COUNT_LIMITS
                const int MIN_ITS_ON_BOUND = 9;//best

                //const int MIN_ITS_ON_BOUND=10;
                //const int MIN_ITS_ON_BOUND=11;
                //const int MIN_ITS_ON_BOUND=8;
                //const int MIN_ITS_ON_BOUND=7;


                countLimits(populTemp[curIndx].data(), &(populLimCount[curIndx]), NVars, Left, Right);
                if (populLimCount[curIndx] > MIN_ITS_ON_BOUND) {
                    ////cout<<"Bounds restart, curBestSol:"<<globalBestFit-optFit<<" FES:"<<NFEval<<"-----------------------------------------------"<<endl;
                    return;
                }

#endif  

#endif            

                popFitTmp[curIndx] = m_eval_fun(populTemp[curIndx], env);
                if (popFitTmp[curIndx] < iterBestFit) {
                    iterBestFit = popFitTmp[curIndx];
                    iterBestIndx = curIndx;
                }
                //if (popFitTmp[curIndx] <= globalBestFit) {
                //    globalBestFit = popFitTmp[curIndx];
                //    memcpy(bestSol, populTemp[curIndx], dimensionality * sizeof(double));
                //    if (globalBestFit - optFit <= MIN_ERROR) {
                //        foundAt = NFEval;
                //        ResultsArray[func_num - 1][runNum][NUM_OF_DUMPS] = foundAt;
                //        SaveBestValues(func_num, runNum);
                //        //cout << "Found Opt----------------------------" << endl;
                //        return;
                //    }
                //}

                if (popFitTmp[curIndx] < Fitmass[curIndx])
                    SaveSuccessCrF(Cr, F, fabs(Fitmass[curIndx] - popFitTmp[curIndx]));

                FindNSaveBest(false, curIndx);
                // SaveBestValues(func_num, runNum);

                //if (alg->evaluations() > alg->maximumEvaluations()) {
                //    //cout << "Max Evals----------------------------" << endl;

                //    return;
                //}
            }
            


            int CHOSEN_INDX = 0;


#ifdef K_MEANS_AS_NEAREST

            std::vector<std::vector<double>> data(popSize, std::vector<double>(dimensionality));

            //arma::mat data(dimensionality, popSize); // n_rows, n_cols
            // The assignments will be stored in this vector.
            for (int i = 0; i < popSize; ++i) {
                for (int j = 0; j < dimensionality; ++j) {
                    data[i][j] = populTemp[i][j];//(i,j) at the i-th row and j-th column
                }
            }
            std::vector<size_t> assignments;
            std::vector<size_t> bestAssignments;
            std::vector<std::vector<double>>  centroids;
            std::vector<std::vector<double>>  bestCentroids;

            //const int MAX_K=5;
            const int MAX_K = 2;

            //const int MIN_POP_SIZE_4_SPLIT = 4;
            //const int MIN_POP_SIZE_4_SPLIT = 10;
            const int MIN_POP_SIZE_4_SPLIT = 20;
            double bestSilhouette = -1;
            int best_k = 2;

            EuclideanDistance dis_metric;
            for (int cand_k = 2; cand_k <= MAX_K; ++cand_k) {
                ofec::armilo::KMeans<EuclideanDistance> k(rnd);
                k.Cluster(data, cand_k, centroids, assignments);

                double silhouetteScore = ofec::metrc::Overall(data, assignments, dis_metric);
                if (silhouetteScore > bestSilhouette) {
                    bestSilhouette = silhouetteScore;
                    best_k = cand_k;
                    bestAssignments = assignments;
                    bestCentroids = centroids;
                }
            }
            double bestCandFit = 1e20;

            double MIN_silhouette = 1 / (4 * sqrt(dimensionality));
            //double MIN_silhouette = 1/(5*sqrt(dimensionality));
            //double MIN_silhouette = 1/(3*sqrt(dimensionality));
            if (bestSilhouette > MIN_silhouette && popSize >= MIN_POP_SIZE_4_SPLIT) {
                for (int cur_k = 0; cur_k < best_k; ++cur_k) {
                    for (int i = 0; i < dimensionality; ++i) {
                        mean_indiv[i] = bestCentroids[cur_k][i];
                        FindLimits(mean_indiv.data(), mean_indiv.data(), NVars, Left, Right, rnd);
                    }
                    double fit_mean = m_eval_fun(mean_indiv, env);

                    if (fit_mean < bestCandFit) {
                        bestCandFit = fit_mean;
                        for (int i = 0; i < dimensionality; ++i) {
                            mean_indivCl0[i] = mean_indiv[i];
                        }
                    }

                    //if (fit_mean < globalBestFit) {
                    //    globalBestFit = fit_mean;
                    //    memcpy(bestSol, mean_indiv, dimensionality * sizeof(double));
                    //}
                    //SaveBestValues(func_num, runNum);
                    //if (globalBestFit - optFit <= MIN_ERROR) {
                    //    foundAt = NFEval;
                    //    ResultsArray[func_num - 1][runNum][NUM_OF_DUMPS] = foundAt;
                    //    //cout << "Found Opt Mean----------------------------" << endl;
                    //    return;
                    //}

                    //Sale All in pop:
                    CHOSEN_INDX = getIndxOfNearest(populTemp, popSize, dimensionality, mean_indiv);
                    if (fit_mean < popFitTmp[CHOSEN_INDX]) {
                        popFitTmp[CHOSEN_INDX] = fit_mean;

                        populTemp[CHOSEN_INDX] = mean_indiv;
                        //                    memcpy(populTemp[CHOSEN_INDX], mean_indiv, dimensionality * sizeof(double));
                    }


                }//for( int cur_k=0; cur_k<best_k; ++cur_k){


            }
            else {//kmaeans failed- use just mean
                getMean(populTemp, popSize, dimensionality, mean_indiv);
                FindLimits(mean_indiv.data(), mean_indiv.data(), NVars, Left, Right, rnd);
                double fit_mean = m_eval_fun(mean_indiv, env);
                //  if (fit_mean < globalBestFit) {
                //      globalBestFit = fit_mean;
                //      memcpy(bestSol, mean_indiv, dimensionality * sizeof(double));
                //  }
                ////  SaveBestValues(func_num, runNum);
                //  if (globalBestFit - optFit <= MIN_ERROR) {
                //      foundAt = NFEval;
                //      ResultsArray[func_num - 1][runNum][NUM_OF_DUMPS] = foundAt;
                //      //cout << "Found Opt Mean----------------------------" << endl;
                //      return;
                //  }
                CHOSEN_INDX = getIndxOfNearest(populTemp, popSize, dimensionality, mean_indiv);
                if (fit_mean < popFitTmp[CHOSEN_INDX]) {
                    popFitTmp[CHOSEN_INDX] = fit_mean;
                    populTemp[CHOSEN_INDX] = mean_indiv;
                    //memcpy(populTemp[CHOSEN_INDX], mean_indiv, dimensionality * sizeof(double));
                }



            }//else{//kmaeans failed- use just mean

            const double MIN_DIST = 1e-9;
            //const int MIN_numOfStagIt=7;//best

            const int MIN_numOfStagIt = 8;
            //const int MIN_numOfStagIt=9;
            //const int MIN_numOfStagIt=6;
            //const int MIN_numOfStagIt=5;
            getMean(populTemp, popSize, dimensionality, mean_indiv);
            if (getDist(mean_indiv_old.data(), mean_indiv.data(), dimensionality) < MIN_DIST) {
                numOfStagIt += 1;
                if (numOfStagIt > MIN_numOfStagIt) {
                    //cout << "Restart dist stagIt, curBestSol:" << globalBestFit - optFit << " FES:" << NFEval << "-----------------------------------------------" << endl;
                    return;
                }
            }
            else {
                numOfStagIt = 0;
            }

            mean_indiv_old = mean_indiv;
            //memcpy(mean_indiv_old, mean_indiv, dimensionality * sizeof(double));
            ////cout<<endl;

#endif     


            ArchSuccess = 0;
            NoArchSuccess = 0;
            NArchUsages = 0;
            for (int curIndx = 0; curIndx != popSize; curIndx++)
            {
                if (popFitTmp[curIndx] <= Fitmass[curIndx])
                {
                    if (ArchUsages[curIndx] == 1)
                    {
                        ArchSuccess += (Fitmass[curIndx] - popFitTmp[curIndx]) / Fitmass[curIndx];
                        NArchUsages += 1;
                    }
                    else
                        NoArchSuccess += (Fitmass[curIndx] - popFitTmp[curIndx]) / Fitmass[curIndx];
                    CopyToArchive(Popul[curIndx].data(), rnd);
                }
                if (popFitTmp[curIndx] <= Fitmass[curIndx])
                {
                    for (int j = 0; j != NVars; j++)
                        Popul[curIndx][j] = populTemp[curIndx][j];
                    Fitmass[curIndx] = popFitTmp[curIndx];
                }
            }
            if (NArchUsages != 0)
            {
                ArchSuccess = ArchSuccess / NArchUsages;
                NoArchSuccess = NoArchSuccess / (popSize - NArchUsages);
                ArchProbs = ArchSuccess / (ArchSuccess + NoArchSuccess);
                ArchProbs = max(0.1, min(0.9, ArchProbs));
                if (ArchSuccess == 0)
                    ArchProbs = 0.5;
            }
            else
                ArchProbs = 0.5;
            int newNInds;
            if (isRestart) {//HomoAtReset
                double shapeConst = 0.1;//jak większe to szybszy spadek, mozna tez pokręcić minimalnym rozmiarem pop
                //double shapeConst = 0.2;//jak większe to szybszy spadek, mozna tez pokręcić minimalnym rozmiarem pop
                //double shapeConst = 0.05;

                //const int MIN_POP=4;

                const int MIN_POP = 20;//best

                //const int MIN_POP=30;
                //const int MIN_POP=10;

                double divider = (alg->maximumEvaluations() - evalsAtStart) / shapeConst;

                double delta = pow((MIN_POP * (alg->maximumEvaluations() - evalsAtStart) / divider - (alg->maximumEvaluations() - evalsAtStart) / divider * popSize), 2) - 4 * (alg->maximumEvaluations() - evalsAtStart) / divider * (MIN_POP - popSize);
                if (delta <= 0) {
                    newNInds = MIN_POP;
                }
                else {
                    double b1 = (-(MIN_POP * (alg->maximumEvaluations() - evalsAtStart) / divider - (alg->maximumEvaluations() - evalsAtStart) / divider * popSize) - sqrt(delta)) / (2 * (MIN_POP - popSize));
                    double a1 = popSize - 1 / b1;
                    newNInds = round(a1 + 1 / ((alg->evaluations() - evalsAtStart) / divider + b1));
                }
            }
            else {
                newNInds = round((NIndsMin - NIndsMax) * pow((double(alg->evaluations()) / double(alg->maximumEvaluations())), (1.0 - double(alg->evaluations()) / double(alg->maximumEvaluations()))) + NIndsMax);
            }

#ifdef OFEC_DATUM_MULTI_POP_H
            {

                std::vector<ofec::Continuous::SolutionType> pop(popSize,
                    { con_pro->numberObjectives(), con_pro->numberConstraints(),con_pro->numberVariables() });
                for (int idx(0); idx < pop.size(); ++idx) {
                    pop[idx].variable().vector() = Popul[idx];
                }

                g_multi_pop.bindPopulation(pop);
                env->algorithm()->datumUpdated(env, g_multi_pop);
            }
#endif

            if (newNInds < NIndsMin)
                newNInds = NIndsMin;
            if (newNInds > NIndsMax)
                newNInds = NIndsMax;
            int newArchSize = round((NIndsMin - NIndsMax) * pow((double(alg->evaluations()) / double(alg->maximumEvaluations())), (1.0 - double(alg->evaluations()) / double(alg->maximumEvaluations()))) + NIndsMax) * ArchiveSizeParam;
            if (newArchSize < NIndsMin)
                newArchSize = NIndsMin;
            ArchiveSize = newArchSize;
            if (CurrentArchiveSize >= ArchiveSize)
                CurrentArchiveSize = ArchiveSize;
            RemoveWorst(popSize, newNInds);
            popSize = newNInds;
            UpdateMemoryCrF();
            SuccessFilled = 0;
            Generation++;



        } while (alg->evaluations() < alg->maximumEvaluations());
        //delete[] est_m;
        //delete[] mean_indiv;
    }

    void Optimizer::Clean()
    {
        Donor.clear();
        Trial.clear();
        Rands.clear();
        for (int i = 0; i != NIndsMax; i++)
        {
            Popul[i].clear();
            populTemp[i].clear();
        }
        for (int i = 0; i != NIndsMax * Int_ArchiveSizeParam; i++)
            Archive[i].clear();
        ArchUsages.clear();
        Archive.clear();
        Popul.clear();
        populTemp.clear();
        Fitmass.clear();
        popFitTmp.clear();
        FitmassCopy.clear();
        BestInd.clear();
        Indexes.clear();
        BackIndexes.clear();
        tempSuccessCr.clear();
        tempSuccessF.clear();
        FGenerated.clear();
        CrGenerated.clear();
        FitDelta.clear();
        MemoryCr.clear();
        MemoryF.clear();
        Weights.clear();
    }


}