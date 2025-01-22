#include "optimizer.h"
#include <math.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <random>

#include "../../../../../../core/environment/environment.h"
#include "../../../../../../core/problem/continuous/continuous.h"

#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS
#include "../../../../../../datum/datum_inclusion.h"

namespace ofec {

    namespace nl_shade_lbc {

        int IntRandom(int target, ofec::Random* rnd)
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


        void qSort1(double* Mass, int low, int high)
        {
            int i = low;
            int j = high;
            double x = Mass[(low + high) >> 1];
            do
            {
                while (Mass[i] < x)    ++i;
                while (Mass[j] > x)    --j;
                if (i <= j)
                {
                    double temp = Mass[i];
                    Mass[i] = Mass[j];
                    Mass[j] = temp;
                    i++;    j--;
                }
            } while (i <= j);
            if (low < j)   qSort1(Mass, low, j);
            if (i < high)  qSort1(Mass, i, high);
        }
        void qSort2int(double* Mass, int* Mass2, int low, int high)
        {
            int i = low;
            int j = high;
            double x = Mass[(low + high) >> 1];
            do
            {
                while (Mass[i] < x)    ++i;
                while (Mass[j] > x)    --j;
                if (i <= j)
                {
                    double temp = Mass[i];
                    Mass[i] = Mass[j];
                    Mass[j] = temp;
                    int temp2 = Mass2[i];
                    Mass2[i] = Mass2[j];
                    Mass2[j] = temp2;
                    i++;    j--;
                }
            } while (i <= j);
            if (low < j)   qSort2int(Mass, Mass2, low, j);
            if (i < high)  qSort2int(Mass, Mass2, i, high);
        }



        void GenerateNextRandUnif(const int num, const int Range, int* Rands, const int Prohib, ofec::Random* rnd)
        {
            for (int j = 0; j != 25; j++)
            {
                bool generateagain = false;
                Rands[num] = IntRandom(Range, rnd);
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
                Rands[num] = IntRandom(Range2, rnd) + Range;
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
        //void SaveBestValues(int funcN, int RunN, double newbestfit)
        //{
        //    double temp = globalbest - fopt;
        //    if (temp <= 10E-8 && ResultsArray[funcN - 1][RunN][16] == MaxFEval)
        //    {
        //        ResultsArray[funcN - 1][RunN][16] = NFEval;
        //    }
        //    for (int stepFEcount = LastFEcount; stepFEcount < 16; stepFEcount++)
        //    {
        //        if (NFEval == int(stepsFEval[stepFEcount] * MaxFEval))
        //        {
        //            if (temp <= 10E-8)
        //                temp = 0;
        //            ResultsArray[funcN - 1][RunN][stepFEcount] = temp;
        //            LastFEcount = stepFEcount;
        //        }
        //    }
        //}

        void FindLimits(double* Ind, double* Parent, int CurNVars, double CurLeft, double CurRight)
        {
            for (int j = 0; j < CurNVars; j++)
            {
                for (int j = 0; j < CurNVars; j++)
                {
                    if (Ind[j] < CurLeft)
                        Ind[j] = (CurLeft + Parent[j]) / 2.0;
                    if (Ind[j] > CurRight)
                        Ind[j] = (CurRight + Parent[j]) / 2.0;
                }
            }
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

            FitNotCalculated = true;
            NInds = newNInds;
            NIndsMax = NInds;
            NIndsMin = 4;
            NVars = newNVars;

            Left = -100;
            Right = 100;
            Cr = 0.9;
            F = 0.5;
            besti = 0;
            Generation = 0;
            CurrentArchiveSize = 0;
            MWLp1 = 3.5;
            MWLp2 = 1.0;
            MWLm = 1.5;
            LBC_fin = 1.5;
            ArchiveSizeParam = NewArchSizeParam;
            Int_ArchiveSizeParam = ceil(ArchiveSizeParam);
            ArchiveSize = NIndsMax * ArchiveSizeParam;

            //for (int steps_k = 0; steps_k != 15; steps_k++)
            //    stepsFEval[steps_k] = pow(double(GNVars), double(steps_k) / 5.0 - 3.0);
            //stepsFEval[15] = 1.0;

            Popul.resize(NIndsMax);
            for (auto& it : Popul) {
                it.resize(NVars);
            }
            PopulTemp = Popul;

            Archive.resize(NIndsMax * Int_ArchiveSizeParam);
            for (auto& it : Archive) {
                it.resize(NVars);
            }
            //Popul = new double* [NIndsMax];
            //for (int i = 0; i != NIndsMax; i++)
            //    Popul[i] = new double[NVars];
            //PopulTemp = new double* [NIndsMax];
            //for (int i = 0; i != NIndsMax; i++)
            //    PopulTemp[i] = new double[NVars];
            //Archive = new double* [NIndsMax * Int_ArchiveSizeParam];
            //for (int i = 0; i != NIndsMax * Int_ArchiveSizeParam; i++)
            //    Archive[i] = new double[NVars];
            FitMass.resize(NIndsMax);
            FitMassTemp.resize(NIndsMax);
            FitMassCopy.resize(NIndsMax);
            FitMassArch.resize(NIndsMax * Int_ArchiveSizeParam);
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
            Weights.resize(NIndsMax);
            MemoryCr.resize(MemorySize);
            MemoryF.resize(MemorySize);
            for (int i = 0; i != MemorySize; i++)
            {
                MemoryCr[i] = 0.9;
                MemoryF[i] = 0.5;
            }
        }
        void Optimizer::SaveSuccessCrF(double Cr, double F, double FitD)
        {
            tempSuccessCr[SuccessFilled] = Cr;
            tempSuccessF[SuccessFilled] = F;
            FitDelta[SuccessFilled] = FitD;
            SuccessFilled++;
        }
        void Optimizer::UpdateMemoryCrF(ofec::Environment* env)
        {
            auto alg = env->algorithm();
            auto MaxFEval = alg->maximumEvaluations();
            if (SuccessFilled != 0)
            {
                double FMWL = LBC_fin + (MWLp1 - LBC_fin) * double(MaxFEval - alg->evaluations()) / (double)MaxFEval;
                double CrMWL = LBC_fin + (MWLp2 - LBC_fin) * double(MaxFEval - alg->evaluations()) / (double)MaxFEval;
                MemoryF[MemoryIter] = (MemoryF[MemoryIter] + MeanWL_general(tempSuccessF.data(), FitDelta.data(), SuccessFilled, FMWL, MWLm)) * 0.5;
                MemoryCr[MemoryIter] = (MemoryCr[MemoryIter] + MeanWL_general(tempSuccessCr.data(), FitDelta.data(), SuccessFilled, CrMWL, MWLm)) * 0.5;
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
        void Optimizer::CopyToArchive(double* RefusedParent, double RefusedFitness, ofec::Random* rnd)
        {
            if (CurrentArchiveSize < ArchiveSize)
            {
                for (int i = 0; i != NVars; i++)
                    Archive[CurrentArchiveSize][i] = RefusedParent[i];
                FitMassArch[CurrentArchiveSize] = RefusedFitness;
                CurrentArchiveSize++;
            }
            else if (ArchiveSize > 0)
            {
                int RandomNum = IntRandom(ArchiveSize, rnd);
                int counter = 0;
                while (FitMassArch[RandomNum] < RefusedFitness)
                {
                    RandomNum = IntRandom(ArchiveSize, rnd);
                    counter++;
                    if (counter == ArchiveSize)
                        break;
                }
                for (int i = 0; i != NVars; i++)
                    Archive[RandomNum][i] = RefusedParent[i];
                FitMassArch[RandomNum] = RefusedFitness;
            }
        }
        void Optimizer::FindNSaveBest(bool init, int ChosenOne)
        {
            if (FitMass[ChosenOne] <= bestfit || init)
            {
                bestfit = FitMass[ChosenOne];
                besti = ChosenOne;
                for (int j = 0; j != NVars; j++)
                    BestInd[j] = Popul[besti][j];
            }
            //if (bestfit < globalbest)
            //    globalbest = bestfit;
        }
        void Optimizer::RemoveWorst(int NInds, int NewNInds)
        {
            int PointsToRemove = NInds - NewNInds;
            for (int L = 0; L != PointsToRemove; L++)
            {
                double WorstFit = FitMass[0];
                int WorstNum = 0;
                for (int i = 1; i != NInds; i++)
                {
                    if (FitMass[i] > WorstFit)
                    {
                        WorstFit = FitMass[i];
                        WorstNum = i;
                    }
                }
                for (int i = WorstNum; i != NInds - 1; i++)
                {
                    for (int j = 0; j != NVars; j++)
                        Popul[i][j] = Popul[i + 1][j];
                    FitMass[i] = FitMass[i + 1];
                }
            }
        }
        inline double Optimizer::GetValue(const int index, const int NInds, const int j)
        {
            if (index < NInds)
                return Popul[index][j];
            return Archive[index - NInds][j];
        }
        void Optimizer::MainCycle(ofec::Environment* env, ofec::Random* rnd)
        {
            ofec::Continuous* con_pro = dynamic_cast<ofec::Continuous*>(env->problem());
            auto alg = env->algorithm();
            auto MaxFEval = alg->maximumEvaluations();
            auto GNVars = con_pro->numberVariables();

            double ArchProbs = 0.5;
            for (int TheChosenOne = 0; TheChosenOne != NInds; TheChosenOne++)
            {
                FitMass[TheChosenOne] = m_eval_fun(Popul[TheChosenOne], env);
                FindNSaveBest(TheChosenOne == 0, TheChosenOne);
                //if (!globalbestinit || bestfit < globalbest)
                //{
                //    globalbest = bestfit;
                //    globalbestinit = true;
                //}
              //  SaveBestValues(func_num, RunN, bestfit);
            }

#ifdef OFEC_DATUM_MULTI_POP_H
            {

                std::vector<ofec::Continuous::SolutionType> pop(NInds,
                    { con_pro->numberObjectives(), con_pro->numberConstraints(),con_pro->numberVariables() });
                for (int idx(0); idx < pop.size(); ++idx) {
                    pop[idx].variable().vector() = Popul[idx];
                }

                g_multi_pop.bindPopulation(pop);
                env->algorithm()->datumUpdated(env, g_multi_pop);
            }
#endif


            std::mt19937 generator_uni_i_2(1e9 * rnd->uniform.next());
            do
            {
                double minfit = FitMass[0];
                double maxfit = FitMass[0];
                for (int i = 0; i != NInds; i++)
                {
                    FitMassCopy[i] = FitMass[i];
                    Indexes[i] = i;
                    if (FitMass[i] >= maxfit)
                        maxfit = FitMass[i];
                    if (FitMass[i] <= minfit)
                        minfit = FitMass[i];
                }
                if (minfit != maxfit)
                    qSort2int(FitMassCopy.data(), Indexes.data(), 0, NInds - 1);
                for (int i = 0; i != NInds; i++)
                    for (int j = 0; j != NInds; j++)
                        if (i == Indexes[j])
                        {
                            BackIndexes[i] = j;
                            break;
                        }
                std::vector<double> FitTemp3(NInds, 0);
                //FitTemp3.resize(NInds);
                for (int i = 0; i != NInds; i++)
                    FitTemp3[i] = exp(-double(i) / (double)NInds);
                std::discrete_distribution<int> ComponentSelector3(FitTemp3.begin(), FitTemp3.end());
                int psizeval = std::max(2.0, NInds * (0.1 / (double)MaxFEval * (double)alg->evaluations() + 0.2));
                for (int TheChosenOne = 0; TheChosenOne != NInds; TheChosenOne++)
                {
                    MemoryCurrentIndex = IntRandom(MemorySize, rnd);
                    Cr = std::min(1.0, std::max(0.0, NormRand(MemoryCr[MemoryCurrentIndex], 0.1, rnd)));
                    do
                    {
                        F = CachyRand(MemoryF[MemoryCurrentIndex], 0.1, rnd);
                    } while (F <= 0);
                    FGenerated[TheChosenOne] = std::min(F, 1.0);
                    CrGenerated[TheChosenOne] = Cr;
                }
                qSort1(CrGenerated.data(), 0, NInds - 1);
                for (int TheChosenOne = 0; TheChosenOne != NInds; TheChosenOne++)
                {
                    for (int Repeat = 0; Repeat != 100; Repeat++)
                    {
                        if (Repeat > 0)
                        {
                            do
                            {
                                F = CachyRand(MemoryF[MemoryCurrentIndex], 0.1, rnd);
                            } while (F <= 0);
                            FGenerated[TheChosenOne] = std::min(F, 1.0);
                        }
                        Rands[0] = Indexes[IntRandom(psizeval, rnd)];
                        for (int i = 0; i != 25 && !CheckGenerated(0, Rands.data(), TheChosenOne); i++)
                            Rands[0] = Indexes[IntRandom(psizeval, rnd)];
                        GenerateNextRandUnif(1, NInds, Rands.data(), TheChosenOne, rnd);
                        if (RandomD(0, 1, rnd) > ArchProbs || CurrentArchiveSize == 0)
                        {
                            Rands[2] = Indexes[ComponentSelector3(generator_uni_i_2)];
                            for (int i = 0; i != 25 && !CheckGenerated(2, Rands.data(), TheChosenOne); i++)
                                Rands[2] = Indexes[ComponentSelector3(generator_uni_i_2)];
                        }
                        else
                            GenerateNextRandUnifOnlyArch(2, NInds, CurrentArchiveSize, Rands.data(), TheChosenOne, rnd);
                        F = FGenerated[TheChosenOne];
                        for (int j = 0; j != NVars; j++)
                            Donor[j] = Popul[TheChosenOne][j] +
                            FGenerated[TheChosenOne] * (GetValue(Rands[0], NInds, j) - Popul[TheChosenOne][j]) +
                            FGenerated[TheChosenOne] * (GetValue(Rands[1], NInds, j) - GetValue(Rands[2], NInds, j));

                        int WillCrossover = IntRandom(NVars, rnd);
                        Cr = CrGenerated[BackIndexes[TheChosenOne]];
                        for (int j = 0; j != NVars; j++)
                        {
                            if (RandomD(0, 1, rnd) < Cr || WillCrossover == j)
                                PopulTemp[TheChosenOne][j] = Donor[j];
                            else
                                PopulTemp[TheChosenOne][j] = Popul[TheChosenOne][j];
                        }
                        bool stopRep = true;
                        for (int j = 0; j != GNVars; j++)
                        {
                            if (PopulTemp[TheChosenOne][j] > Right)
                                stopRep = false;
                            if (PopulTemp[TheChosenOne][j] < Left)
                                stopRep = false;
                        }
                        if (stopRep)
                            break;
                    }
                    FindLimits(PopulTemp[TheChosenOne].data(), Popul[TheChosenOne].data(), NVars, Left, Right);
                    FitMassTemp[TheChosenOne] = m_eval_fun(PopulTemp[TheChosenOne], env);
                    //if (FitMassTemp[TheChosenOne] <= globalbest)
                    //    globalbest = FitMassTemp[TheChosenOne];

                    if (FitMassTemp[TheChosenOne] < FitMass[TheChosenOne])
                        SaveSuccessCrF(Cr, F, fabs(FitMass[TheChosenOne] - FitMassTemp[TheChosenOne]));
                    FindNSaveBest(false, TheChosenOne);
                    //SaveBestValues(func_num, RunN, bestfit);
                }
                for (int TheChosenOne = 0; TheChosenOne != NInds; TheChosenOne++)
                {
                    if (FitMassTemp[TheChosenOne] <= FitMass[TheChosenOne])
                    {
                        CopyToArchive(Popul[TheChosenOne].data(), FitMass[TheChosenOne], rnd);
                        for (int j = 0; j != NVars; j++)
                            Popul[TheChosenOne][j] = PopulTemp[TheChosenOne][j];
                        FitMass[TheChosenOne] = FitMassTemp[TheChosenOne];
                    }
                }



#ifdef OFEC_DATUM_MULTI_POP_H
                {

                    std::vector<ofec::Continuous::SolutionType> pop(NInds,
                        { con_pro->numberObjectives(), con_pro->numberConstraints(),con_pro->numberVariables() });
                    for (int idx(0); idx < pop.size(); ++idx) {
                        pop[idx].variable().vector() = Popul[idx];
                    }

                    g_multi_pop.bindPopulation(pop);
                    env->algorithm()->datumUpdated(env, g_multi_pop);
                }
#endif

                int newNInds = round((NIndsMin - NIndsMax) * pow((double(alg->evaluations()) / double(MaxFEval)), (1.0 - double(alg->evaluations()) / double(MaxFEval))) + NIndsMax);
                if (newNInds < NIndsMin)
                    newNInds = NIndsMin;
                if (newNInds > NIndsMax)
                    newNInds = NIndsMax;
                int newArchSize = round((NIndsMin - NIndsMax) * pow((double(alg->evaluations()) / double(MaxFEval)), (1.0 - double(alg->evaluations()) / double(MaxFEval))) + NIndsMax) * ArchiveSizeParam;
                if (newArchSize < NIndsMin)
                    newArchSize = NIndsMin;
                ArchiveSize = newArchSize;
                if (CurrentArchiveSize >= ArchiveSize)
                    CurrentArchiveSize = ArchiveSize;
                RemoveWorst(NInds, newNInds);
                NInds = newNInds;
                UpdateMemoryCrF(env);
                SuccessFilled = 0;
                Generation++;
            } while (alg->evaluations() < MaxFEval);
        }
        void Optimizer::Clean()
        {
            Donor.clear();
            Trial.clear();
            Rands.clear();
            for (int i = 0; i != NIndsMax; i++)
            {
                Popul[i].clear();
                PopulTemp[i].clear();
            }
            for (int i = 0; i != NIndsMax * Int_ArchiveSizeParam; i++)
                Archive[i].clear();
            Archive.clear();
            Popul.clear();
            PopulTemp.clear();
            FitMass.clear();
            FitMassTemp.clear();
            FitMassCopy.clear();
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


}