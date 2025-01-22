#include "randomize.h"
#include <cmath>
using namespace ofec::eax_tsp;

TRandom::TRandom(const std::shared_ptr<Random> &rnd) : m_random(rnd){}
TRandom::~TRandom(){}

int TRandom::Integer( int minNumber, int maxNumber ){
    using namespace ofec;
    return m_random->uniform.nextNonStd<int>(minNumber, maxNumber + 1);
    //return m_random.uniform.nextNonStd<int>(minNumber, maxNumber+1);
 //   return minNumber + (rand() % (maxNumber - minNumber + 1));
}

double TRandom::Double( double minNumber, double maxNumber ){
    using namespace ofec;
    return m_random->uniform.nextNonStd<double>(minNumber, maxNumber);
   // return minNumber + rand() % (int)(maxNumber - minNumber);
}

void TRandom::permutation( std::vector<int>& arr, int numOfElement, int numOfSample ){

    if( numOfElement <= 0 ) return;
    using namespace ofec;
    int i, j, k, r;
//    int *b = new int[numOfElement];
    std::vector<int> b(numOfElement);
    for(j=0;j<numOfElement;j++) b[j]=0;
    for(i=0;i<numOfSample;i++){
        //r=rand()%(numOfElement-i);
        r = m_random->uniform.nextNonStd<int>(0, numOfElement - i);
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

double TRandom::normalDistribution( double mu, double sigma ){
    double U1,U2,X;
    double PI = 3.1415926;
    while( 1 ){
        U1 = this->Double( 0.0, 1.0 );
        if( U1 != 0.0 ) break;
    }
    U2 = this->Double( 0.0, 1.0 );
    X = sqrt(-2.0*log(U1)) * cos(2*PI*U2);
    return( mu + sigma*X );
}

void TRandom::shuffle(std::vector<int>& arr , int numOfElement ){
   // int *a = new int[numOfElement];
    std::vector<int> a(numOfElement);
 //   int *b = new int[numOfElement];
    std::vector<int> b(numOfElement);
    this->permutation( b, numOfElement, numOfElement );
    for( int i = 0; i < numOfElement; ++i ) a[ i ] = arr[ i ];
    for( int i = 0; i < numOfElement; ++i ) arr[ i ] = a[ b[ i ] ];
 //   delete [] a;
  //  delete [] b;
}
