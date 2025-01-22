#include "randomize.h"
#include <cmath>
#include "environment.h"
using namespace ofec::eax_tsp2;


void TRandom::setEnv(TEnvironment* env) {
    m_env = env;
}
TRandom::TRandom(){}
TRandom::~TRandom(){}

int TRandom::Integer( int minNumber, int maxNumber, Random* rnd){
    using namespace ofec;
    m_env->testRandom(rnd, "TRandom::Integer before");
    int value = rnd->uniform.nextNonStd<int>(minNumber, maxNumber + 1);
    m_env->testRandom(rnd, "TRandom::Integer after");
    return value;
    //return m_random.uniform.nextNonStd<int>(minNumber, maxNumber+1);
 //   return minNumber + (rand() % (maxNumber - minNumber + 1));
}

double TRandom::Double( double minNumber, double maxNumber, Random* rnd){
    using namespace ofec;
    m_env->testRandom(rnd, "TRandom::Double before");
    double value = rnd->uniform.nextNonStd<double>(minNumber, maxNumber);
    m_env->testRandom(rnd, "TRandom::Double after");
    return value;
   // return minNumber + rand() % (int)(maxNumber - minNumber);
}

void TRandom::permutation( std::vector<int>& arr, int numOfElement, int numOfSample, Random* rnd){

    if( numOfElement <= 0 ) return;
    using namespace ofec;
    int i, j, k, r;
//    int *b = new int[numOfElement];
    std::vector<int> b(numOfElement);
    for(j=0;j<numOfElement;j++) b[j]=0;
    for(i=0;i<numOfSample;i++){
        //r=rand()%(numOfElement-i);

        m_env->testRandom(rnd, "TRandom::permutation before");
        r = rnd->uniform.nextNonStd<int>(0, numOfElement - i);
        m_env->testRandom(rnd, "TRandom::permutation after");
        k=0;
        for(j=0;j<=r;j++){
            while(b[k]==1) ++k;
            k++;
        }
        arr[i]=k-1;
        b[k-1]=1;
    }
  //  delete [] b;
}

double TRandom::normalDistribution( double mu, double sigma, Random* rnd){
    double U1,U2,X;
    double PI = 3.1415926;
    while( 1 ){
        U1 = this->Double( 0.0, 1.0 , rnd);
        if( U1 != 0.0 ) break;
    }
    U2 = this->Double( 0.0, 1.0, rnd);
    X = sqrt(-2.0*log(U1)) * cos(2*PI*U2);
    return( mu + sigma*X );
}

void TRandom::shuffle(std::vector<int>& arr , int numOfElement, Random* rnd){
   // int *a = new int[numOfElement];
    std::vector<int> a(numOfElement);
 //   int *b = new int[numOfElement];
    std::vector<int> b(numOfElement);
    this->permutation( b, numOfElement, numOfElement, rnd);
    for( int i = 0; i < numOfElement; ++i ) a[ i ] = arr[ i ];
    for( int i = 0; i < numOfElement; ++i ) arr[ i ] = a[ b[ i ] ];
 //   delete [] a;
  //  delete [] b;
}
