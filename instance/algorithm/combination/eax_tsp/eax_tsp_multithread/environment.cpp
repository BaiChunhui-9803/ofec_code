#include "environment.h"


#include <cmath> 
#include <iostream>
#include "cross.h"
#include "../../../../../core/problem/problem.h"

using namespace ofec::eax_tsp_mt;

//TEnvironment::TEnvironment()
//{
//   // fEvaluator.reset(new TEvaluator());
//    fFlagC.resize(10);
//}

//void TEnvironment::setPopulation()
//{
//    int N = fEvaluator->Ncity;
//
//    fIndexForMating.resize(Npop + 1);
//   // tCurPop.resize(Npop);
// //   for (int i = 0; i < Npop; ++i) tCurPop[i].define(N);
//
//    tCross->Npop = Npop;
//}
//

void TEnvironment::define(ofec::Environment *env, const std::shared_ptr<Random> &rnd)
{
    m_random = rnd;
    fEvaluator.reset(new TEvaluator());
    fEvaluator->setInstance(env);
    optimum = -1;
    //if (pro->optimaBase()->numberObjectives() > 0) {
    //    optimum = pro->optimaBase()->objective(0)[0];
    //}
    m_maxGeneration = -1;

    int N = fEvaluator->Ncity;

   // Npop = 300;
    //fIndexForMating = new int [ Npop + 1 ];
  //  fIndexForMating.resize(Npop + 1);
  //  tCurPop.resize(Npop);
   // for ( int i = 0; i < Npop; ++i ) tCurPop[i].define( N );

    gBestValue = -1;
    gBest.define(N);

    tBest.define( N );
    tCross.reset( new TCross( N));
    tCross->eval = fEvaluator;                 
    tCross->Npop = Npop;

    tKopt.reset(new TKopt( N ));
    tKopt->eval = fEvaluator;
    tKopt->tRand = tCross->tRand;
    tKopt->setInvNearList();


    fEdgeFreq.resize(N);
    for (auto& it : fEdgeFreq) it.resize(N);
}




//void TEnvironment::doItOrigin()
//{
//
//    int m_iter(0);
//
//    this->initPop();
//    this->init();
//    this->getEdgeFreq();
//
//
//    while (true)
//    {
//        this->setAverageBest();
//
//        if (gBestValue == -1 || fBestValue < gBestValue)
//        {
//            gBestValue = fBestValue;
//            gBest = tBest;
//            // printf("find better solution %d\n", gBestValue);
//            if (gBestValue <= this->optimum)
//            {
//                printf("Find optimal solution %lf, exit\n", gBestValue);
//                this->terminate = true;
//                break;
//            }
//        }
//        if (fCurNumOfGen % 50 == 0)
//        {
//            printf("%d:\t%lf\t%lf\n", fCurNumOfGen, fBestValue, fAverageValue);
//            // record time every 50 gens
//            //this->fTimeEnd = clock();
//            //duration = (int)((double)(this->fTimeEnd - this->fTimeStart) / (double)CLOCKS_PER_SEC);
//            //if (duration >= tmax)
//            //    break;
//        }
//
//        if (this->terminationCondition())
//            break;
//
//        this->selectForMating();
//        ++fCurNumOfGen;
//
//        ++m_iter;
//    }
//}
//


void TEnvironment::init()
{
    //m_generation = 0;
    fAccumurateNumCh = 0;
    fCurNumOfGen = 0;
    fStagBest = 0;
    fMaxStagBest = 0;
    fStage = 1;             /* Stage I */
    fFlagC[0] = 4;          /* Diversity preservation: 1:Greedy, 2:--- , 3:Distance, 4:Entropy (see Section 4) */
    fFlagC[1] = 1;          /* Eset Type: 1:Single-AB, 2:Block2 (see Section 3) */

    m_stagnation_iteration = 0;
    m_curBestCost = std::numeric_limits<int>::max();
} 



bool TEnvironment::terminationCondition(const std::vector<TIndi*>& curPop){

    double curBestVal = std::numeric_limits<double>::max();

    for (auto& it : curPop) {
        if (curBestVal > it->fEvaluationValue) {
            curBestVal = it->fEvaluationValue;
        }
    }

    if (m_curBestCost > curBestVal) {
        m_stagnation_iteration = 0;
        m_curBestCost = curBestVal;
    }
    else {
        ++m_stagnation_iteration;
    }
    if (m_maxStagnationIter != -1) {
        if (m_stagnation_iteration >= m_maxStagnationIter) return true;
    }
    if (optimum != -1) {
        if (curBestVal <= optimum) return true;
    }
    if (m_maxGeneration != -1) {
        if (fCurNumOfGen >= m_maxGeneration) return true;
    }

    if (fAverageValue - fBestValue < 0.001)  return true;
    if (fStage == 1) /* Stage I */
    {
        /* 1500/N_ch (See Section 2.2) */
        if (fStagBest == int(1500 / Nch) && fMaxStagBest == 0)
        {
            /* fMaxStagBest = G/10 (See Section 2.2) */
            fMaxStagBest = int(fCurNumOfGen / 10);
        }
        /* Terminate Stage I (proceed to Stage II) */
        else if (fMaxStagBest != 0 && fMaxStagBest <= fStagBest) {
            fStagBest = 0;
            fMaxStagBest = 0;
            fCurNumOfGen1 = fCurNumOfGen;
            fFlagC[1] = 2;
            fStage = 2;
        }
        return false;
    }
    if (fStage == 2) /* Stage II */
    {
        /* 1500/N_ch */
        if (fStagBest == int(1500 / Nch) && fMaxStagBest == 0)
        {
            /* fMaxStagBest = G/10 (See Section 2.2) */
            fMaxStagBest = int((fCurNumOfGen - fCurNumOfGen1) / 10);
        }
        /* Terminate Stage II and GA */
        else if (fMaxStagBest != 0 && fMaxStagBest <= fStagBest)
        {
            return true;
        }
        return false;
    }

    return true;
}


void TEnvironment::setAverageBest(const std::vector<TIndi*>& curPop)
{
    double stockBest = tBest.fEvaluationValue;
    fAverageValue = 0.0;
    fBestIndex = 0;
    fBestValue = curPop[0]->fEvaluationValue;
    for(int i = 0; i < Npop; ++i )
    {
        fAverageValue += curPop[i]->fEvaluationValue;
        if(curPop[i]->fEvaluationValue < fBestValue )
        {
            fBestIndex = i;
            fBestValue = curPop[i]->fEvaluationValue;
        }
    }
    tBest = *curPop[ fBestIndex ];
    fAverageValue /= (double)Npop;
    if( tBest.fEvaluationValue < stockBest )
    {
        fStagBest = 0;
        fBestNumOfGen = fCurNumOfGen;
        fBestAccumeratedNumCh = fAccumurateNumCh;
    }
    else ++fStagBest;
}

//void TEnvironment::initPop(){
//    for ( int i = 0; i < Npop; ++i )
//    {
//        tKopt->makeRandSol(tCurPop[i], m_random.get()); /* Make a random tour */
//        tKopt->doIt(tCurPop[i], m_random.get());        /* Apply the local search with the 2-opt neighborhood */
//    }
//}

//void TEnvironment::selectForMating()
//{
//    /* fIndexForMating[] <-- a random permutation of 0, ..., fNumOfPop-1 */
//    tCross->tRand->permutation( fIndexForMating, Npop, Npop, m_random.get() );
//    fIndexForMating[ Npop ] = fIndexForMating[ 0 ];
//}

//void TEnvironment::generateKids( int s )
//{
//    /* Note: tCurPop[fIndexForMating[s]] is replaced with a best offspring solutions in tCorss->DoIt(). 
//     fEegeFreq[][] is also updated there. */
//    tCross->setParents( tCurPop[fIndexForMating[s]], tCurPop[fIndexForMating[s+1]], fFlagC, Nch, m_random.get() );  
//    tCross->doIt( tCurPop[fIndexForMating[s]], tCurPop[fIndexForMating[s+1]], Nch, 1, fFlagC, fEdgeFreq, m_random.get());
//    fAccumurateNumCh += tCross->fNumOfGeneratedCh;
//}

void TEnvironment::clearEdgeFreq() {
    int  k0, k1, N = fEvaluator->Ncity;
    for (int j1 = 0; j1 < N; ++j1)
        for (int j2 = 0; j2 < N; ++j2)
            fEdgeFreq[j1][j2] = 0;
}
void TEnvironment::calEdgeFreq(std::vector< std::vector<std::vector<int>>*>& links) {
    int  k0, k1, N = fEvaluator->Ncity;
    for (int j1 = 0; j1 < N; ++j1)
        for (int j2 = 0; j2 < N; ++j2)
            fEdgeFreq[j1][j2] = 0;

    for (int i = 0; i < links.size(); ++i)
        for (int j = 0; j < N; ++j) {
            k0 = (*links[i])[j][0];
            k1 = (*links[i])[j][1];
            ++fEdgeFreq[j][k0];
            ++fEdgeFreq[j][k1];
        }
}

//void TEnvironment::getEdgeFreq()
//{
//    int  k0, k1, N = fEvaluator->Ncity;
//    for( int j1 = 0; j1 < N; ++j1 )
//        for( int j2 = 0; j2 < N; ++j2 )
//            fEdgeFreq[ j1 ][ j2 ] = 0;
//
//    for( int i = 0; i < Npop; ++i )
//        for(int j = 0; j < N; ++j ){
//            k0 = tCurPop[ i ].fLink[ j ][ 0 ];
//            k1 = tCurPop[ i ].fLink[ j ][ 1 ];
//            ++fEdgeFreq[ j ][ k0 ];
//            ++fEdgeFreq[ j ][ k1 ];
//        }
//}
