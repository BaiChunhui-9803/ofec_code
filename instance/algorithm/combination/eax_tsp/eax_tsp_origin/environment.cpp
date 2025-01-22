#include "environment.h"


#include <cmath> 
#include <iostream>
#include "cross.h"
#include "../../../../../core/problem/problem.h"

using namespace ofec::eax_tsp;

//void MakeRandSol(TEvaluator* eval, TIndi& indi);
//void Make2optSol(TEvaluator* eval, TIndi& indi);

TEnvironment::TEnvironment()
{
    fEvaluator.reset(new TEvaluator());
    fFlagC.resize(10);
}

TEnvironment::~TEnvironment()
{

   // int N = fEvaluator->Ncity;

    //delete [] fIndexForMating;
    //delete [] tCurPop;
    //delete fEvaluator;
    //delete tCross;


    //for( int i = 0; i < N; ++i ) delete [] fEdgeFreq[ i ];
    //delete [] fEdgeFreq;
}

void ofec::eax_tsp::TEnvironment::setPopulation()
{
    int N = fEvaluator->Ncity;

    fIndexForMating.resize(Npop + 1);
    tCurPop.resize(Npop);
    for (int i = 0; i < Npop; ++i) tCurPop[i].define(N);

    tCross->Npop = Npop;
}


//void TEnvironment::updateLEDDFile()
//{
//        
//    std::string curfilename;
//    curfilename = "tspPop_"+ std::to_string(m_RunId) + "_solutions.txt";
//    std::cout <<m_filepath<<"\t" << curfilename << std::endl;
//
//    m_output_popSol.open(m_filepath + "/" + curfilename);
//
//    curfilename = "tspPop_" + std::to_string(m_RunId) + "_attributes.txt";
//    m_output_solAttributes.open(m_filepath + "/" + curfilename);
//    m_output_solAttributes << "solution_id\trun_id\tgen_num\tobjective" << std::endl;
//   
//}
//
//void TEnvironment::outputLEDDInfo()
//{
//    // ++m_generation;
//    for (int i = 0; i < Npop; ++i) {
//        fEvaluator->transferSol(tCurPop[i], m_outputSol);
//        ++m_solId;
//        m_output_popSol << m_solId << "\t";
//        for (int idx(0); idx < m_outputSol.size(); ++idx) {
//            if (idx) {
//                m_output_popSol << ", ";
//            }
//            m_output_popSol << m_outputSol[idx] - 1;
//        }
//        m_output_popSol << std::endl;
//        m_output_solAttributes << m_solId << "\t" << m_RunId << "\t" << m_generation << "\t" << tCurPop[i].fEvaluationValue << std::endl;
//    }
//}

//void TEnvironment::setPopulation(int numPop) {
//
//}

void TEnvironment::define(ofec::Environment *env, const std::shared_ptr<Random> &rnd)
{

    fEvaluator.reset(new TEvaluator());
    fEvaluator->setInstance(env);
    optimum = -1;
    //if (pro->optimaBase()->numberObjectives() > 0) {
    //    optimum = pro->optimaBase()->objective(0)[0];
    //}
    m_maxGeneration = -1;

  //  m_N= 
    int N = fEvaluator->Ncity;


   // Npop = 300;
    //fIndexForMating = new int [ Npop + 1 ];
    fIndexForMating.resize(Npop + 1);
    tCurPop.resize(Npop);
    for ( int i = 0; i < Npop; ++i ) tCurPop[i].define( N );

    gBestValue = -1;
    gBest.define(N);

    tBest.define( N );
    tCross.reset( new TCross( N ,rnd));
    tCross->eval = fEvaluator;                 
    tCross->Npop = Npop;

    tKopt.reset(new TKopt( N ));
    tKopt->eval = fEvaluator;
    tKopt->tRand = tCross->tRand;
    tKopt->setInvNearList();


    fEdgeFreq.resize(N);
    for (auto& it : fEdgeFreq) it.resize(N);
//    fEdgeFreq = new int* [ N ]; 
//    for( int i = 0; i < N; ++i ) fEdgeFreq[ i ] = new int [ N ];
 //   this->fTimeStart = clock();
}




void TEnvironment::doItOrigin()
{

    int m_iter(0);
    //m_pop_filename = std::to_string(m_RunId) + ".txt";
    //m_output_pop.open(m_filepath+"/"+m_pop_filename);

    //m_popBest_filename = std::to_string(m_RunId) + "_best.txt";
    //m_output_popBest.open(m_filepath + "/" + m_popBest_filename);

    this->initPop();
    this->init();
    this->getEdgeFreq();
   // this->fTimeEnd = clock();
   // int duration = (int)((double)(this->fTimeEnd - this->fTimeStart) / (double)CLOCKS_PER_SEC);

    while (true)
    {
        this->setAverageBest();
      //  m_output_pop << "iteration\t" << m_iter << "\tbest\t" << tBest.fEvaluationValue << std::endl;
  //      for (int s = 0; s < Npop; ++s) {
//            fEvaluator->writeTo(m_output_pop, tCurPop[s]);
  //      }
   //     m_output_popBest << "iteration\t" << m_iter << "\tbest\t" << tBest.fEvaluationValue << std::endl;
  //      fEvaluator->writeTo(m_output_popBest, tBest);


        if (gBestValue == -1 || fBestValue < gBestValue)
        {
            gBestValue = fBestValue;
            gBest = tBest;
            // printf("find better solution %d\n", gBestValue);
            if (gBestValue <= this->optimum)
            {
                printf("Find optimal solution %lf, exit\n", gBestValue);
                this->terminate = true;
                break;
            }
        }
        if (fCurNumOfGen % 50 == 0)
        {
            printf("%d:\t%lf\t%lf\n", fCurNumOfGen, fBestValue, fAverageValue);
            // record time every 50 gens
            //this->fTimeEnd = clock();
            //duration = (int)((double)(this->fTimeEnd - this->fTimeStart) / (double)CLOCKS_PER_SEC);
            //if (duration >= tmax)
            //    break;
        }

        if (this->terminationCondition())
            break;

        this->selectForMating();
        //for (int s = 0; s < Npop; ++s) {
        //    this->generateKids(s);
        //    m_popChildren_filename = std::to_string(m_RunId) + "_" + std::to_string(m_iter) + "_" + std::to_string(s) + ".txt";
        //    m_output_popChildren.open(m_filepath + "/" + m_popChildren_filename);
        //    m_output_popChildren << tCross->Localsearch_indis().size() << std::endl;
        //    for (auto& it : tCross->Localsearch_indis()) {
        //  //      fEvaluator->writeTo(m_output_popChildren, it);
        //    }
        //    m_output_popChildren.close();
        //}


        ++fCurNumOfGen;

        ++m_iter;
    }

    //if (duration >= tmax)
    //    this->terminate = true;

    //++m_RunId;

    //m_output_pop.close();
    //m_output_popBest.close();
}



void TEnvironment::doIt() {


    //{
    //    std::string curfilename;
    //    curfilename = "tspPop_" + std::to_string(m_RunId) + "_solutions.txt";
    //    std::cout << m_filepath << "\t" << curfilename << std::endl;

    //    m_output_popSol.open(m_filepath + "/" + curfilename);

    //    curfilename = "tspPop_" + std::to_string(m_RunId) + "_attributes.txt";
    //    m_output_solAttributes.open(m_filepath + "/" + curfilename);
    //    m_output_solAttributes << "solution_id\trun_id\tgen_num\tobjective" << std::endl;

    //}



    this->initPop();
    this->init();
    this->getEdgeFreq();

    //this->writeBest();

    while (true)
    {
        this->setAverageBest();
        if (this->terminationCondition()) break;

        this->selectForMating();

        for (int s = 0; s < Npop; ++s)
        {
            this->generateKids(s);
            //  this->selectForSurvival(s);
        }

       // this->writeBest();

        ++fCurNumOfGen;
    }

    // this->fTimeEnd = clock();
}

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



bool TEnvironment::terminationCondition()
{
    double curBestVal = std::numeric_limits<double>::max();
    
    for (auto& it : tCurPop) {
        if (curBestVal > it.fEvaluationValue) {
            curBestVal = it.fEvaluationValue;
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

    if ( fAverageValue - fBestValue < 0.001 )  return true;
    if (fStage == 1) /* Stage I */
    {
        /* 1500/N_ch (See Section 2.2) */
        if( fStagBest == int(1500/Nch) && fMaxStagBest == 0)
        {
            /* fMaxStagBest = G/10 (See Section 2.2) */
            fMaxStagBest = int(fCurNumOfGen / 10);
        }
        /* Terminate Stage I (proceed to Stage II) */
        else if( fMaxStagBest != 0 && fMaxStagBest <= fStagBest ){ 
            fStagBest = 0;
            fMaxStagBest = 0;
            fCurNumOfGen1 = fCurNumOfGen;
            fFlagC[ 1 ] = 2;
            fStage = 2;
        }
        return false;
    }
    if (fStage == 2) /* Stage II */
    {
        /* 1500/N_ch */
        if( fStagBest == int(1500/Nch) && fMaxStagBest == 0 )
        {
            /* fMaxStagBest = G/10 (See Section 2.2) */
            fMaxStagBest = int( (fCurNumOfGen - fCurNumOfGen1) / 10 );
        }
        /* Terminate Stage II and GA */
        else if( fMaxStagBest != 0 && fMaxStagBest <= fStagBest )
        {
            return true;
        }
        return false;
    }

    return true;
}

void TEnvironment::setAverageBest()
{
    double stockBest = tBest.fEvaluationValue;
    fAverageValue = 0.0;
    fBestIndex = 0;
    fBestValue = tCurPop[0].fEvaluationValue;
    for(int i = 0; i < Npop; ++i )
    {
        fAverageValue += tCurPop[i].fEvaluationValue;
        if( tCurPop[i].fEvaluationValue < fBestValue )
        {
            fBestIndex = i;
            fBestValue = tCurPop[i].fEvaluationValue;
        }
    }
    tBest = tCurPop[ fBestIndex ];
    fAverageValue /= (double)Npop;
    if( tBest.fEvaluationValue < stockBest )
    {
        fStagBest = 0;
        fBestNumOfGen = fCurNumOfGen;
        fBestAccumeratedNumCh = fAccumurateNumCh;
    }
    else ++fStagBest;
}

void TEnvironment::initPop(){
    for ( int i = 0; i < Npop; ++i )
    {
        tKopt->makeRandSol(tCurPop[i]); /* Make a random tour */
        tKopt->doIt(tCurPop[i]);        /* Apply the local search with the 2-opt neighborhood */
    }
}

void TEnvironment::selectForMating()
{
    /* fIndexForMating[] <-- a random permutation of 0, ..., fNumOfPop-1 */
    tCross->tRand->permutation( fIndexForMating, Npop, Npop );
    fIndexForMating[ Npop ] = fIndexForMating[ 0 ];
}

void TEnvironment::generateKids( int s )
{
    /* Note: tCurPop[fIndexForMating[s]] is replaced with a best offspring solutions in tCorss->DoIt(). 
     fEegeFreq[][] is also updated there. */
    tCross->setParents( tCurPop[fIndexForMating[s]], tCurPop[fIndexForMating[s+1]], fFlagC, Nch );  
    tCross->doIt( tCurPop[fIndexForMating[s]], tCurPop[fIndexForMating[s+1]], Nch, 1, fFlagC, fEdgeFreq );
    fAccumurateNumCh += tCross->fNumOfGeneratedCh;
}

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

void TEnvironment::getEdgeFreq()
{
    int  k0, k1, N = fEvaluator->Ncity;
    for( int j1 = 0; j1 < N; ++j1 )
        for( int j2 = 0; j2 < N; ++j2 )
            fEdgeFreq[ j1 ][ j2 ] = 0;

    for( int i = 0; i < Npop; ++i )
        for(int j = 0; j < N; ++j ){
            k0 = tCurPop[ i ].fLink[ j ][ 0 ];
            k1 = tCurPop[ i ].fLink[ j ][ 1 ];
            ++fEdgeFreq[ j ][ k0 ];
            ++fEdgeFreq[ j ][ k1 ];
        }
}

void TEnvironment::printOn()
{
    //printf("Total time: %d\n", duration);
    //printf("bestval = %d, optimum = %d \n", gBestValue, this->optimum);
    //fEvaluator->writeToStdout(gBest);
    //if (gBestValue != -1 && gBestValue <= this->optimum)
    //    printf("Successful\n");
    //else
    //    printf("Unsuccessful\n");
    //fflush(stdout);
}
//
//void TEnvironment::writeBest()
//{
//    ++m_generation;
//    for (int i = 0; i < Npop; ++i) {
//
//        {
//            auto& sol(m_outputSol);
//            auto& indi(tCurPop[i]);
//            sol.resize(fEvaluator->Ncity);
//            //std::vector<int> Array(Ncity);
//            auto& Array(sol);
//            int curr = 0, st = 0, count = 0, pre = -1, next;
//            while (1)
//            {
//                Array[count++] = curr + 1;
//                if (count > fEvaluator->Ncity)
//                {
//                    //  printf("Invalid\n");
//                    return;
//                }
//                if (indi.fLink[curr][0] == pre)
//                    next = indi.fLink[curr][1];
//                else
//                    next = indi.fLink[curr][0];
//
//                pre = curr;
//                curr = next;
//                if (curr == st)
//                    break;
//            }
//        }
//
//      //  fEvaluator->transferSol(tCurPop[i], m_outputSol);
//        ++m_solId;
//        m_output_popSol << m_solId << "\t";
//        for (int idx(0); idx < m_outputSol.size(); ++idx) {
//            if (idx) {
//                m_output_popSol << ", ";
//            }
//            m_output_popSol << m_outputSol[idx] - 1;
//        }
//        m_output_popSol << std::endl;
//        m_output_solAttributes << m_solId << "\t" << m_RunId << "\t" << m_generation << "\t" << tCurPop[i].fEvaluationValue << std::endl;
//    }
//
//
//    //FILE *fp;
//    //char filename[ 80 ];
//
//    //sprintf( filename, "bestSolution.txt" );
//    //fp = fopen( filename, "a");
//    //fEvaluator->writeTo( fp, gBest );
//    //fclose( fp );
//}
