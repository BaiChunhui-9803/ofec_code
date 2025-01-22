#include "sort.h"
#include <algorithm>
#include "environment.h"

//using namespace eax_tsp2;

void ofec::eax_tsp2::TSort::setEnv(TEnvironment* env) {
    m_env = env;
}


void ofec::eax_tsp2::TSort::selectionSort(std::vector<int>& Arg, int l, int r){
    int id;
    for(int i=l;i<r;++i){
        id=i;
        for(int j=i+1;j<=r;++j)
            if(Arg[j]<Arg[id]) id=j;
        std::swap(Arg[i], Arg[id]);
        //eax_tsp2::swap(Arg[i], Arg[id]);
    }
}

int ofec::eax_tsp2::TSort::partition(std::vector<int>& Arg, int l, int r, Random *rnd){
    using namespace ofec;
//    int id=l+rand()%(r-l+1);
    m_env->testRandom(rnd, "TSort::partition before");
    int id = rnd->uniform.nextNonStd<int>(l, r + 1);
    m_env->testRandom(rnd, "TSort::partition after");
    std::swap(Arg[l], Arg[id]);
    id=l;
    for(int i=l+1;i<=r;++i)
        if(Arg[i]<Arg[l]) std::swap(Arg[++id], Arg[i]);
    std::swap(Arg[l], Arg[id]);
    return id;
}

void ofec::eax_tsp2::TSort::quickSort(std::vector<int>& Arg, int l, int r, Random *rnd){
    if(l<r){
        if(r-l<20)
        {
            selectionSort(Arg, l, r);
            return ;
        }
        int mid= partition(Arg, l, r, rnd);
        quickSort(Arg, l, mid-1, rnd);
        quickSort(Arg, mid+1, r, rnd);
    }
}

ofec::eax_tsp2::TSort::TSort(){}
ofec::eax_tsp2::TSort::~TSort(){}

void ofec::eax_tsp2::TSort::index( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd ){
    int indexBest = 0;
    double valueBest;
    //int *checked = new int [ numOfArg ];
    std::vector<int> checked(numOfArg);
    for( int i = 0 ; i < numOfArg ; ++i ) checked[ i ] = 0;
    for( int i = 0; i < numOfOrd; ++i ){
        valueBest = 99999999999.9;
        for( int j = 0; j < numOfArg; ++j ){
            if( ( Arg[j] < valueBest ) && checked[j]==0){
                valueBest = Arg[j];
                indexBest = j;
            }
        }
        indexOrderd[ i ]=indexBest;
        checked[ indexBest ]=1;
    }
  //  delete [] checked;
}

void ofec::eax_tsp2::TSort::indexB(const std::vector<double>& Arg, int numOfArg, std::vector<int>& indexOrderd, int numOfOrd )
{
    int indexBest = 0;
    double valueBest;
    //int *checked = new int [ numOfArg ];
    std::vector<int> checked(numOfArg);
    for( int i = 0 ; i < numOfArg ; ++i ) checked[ i ] = 0;
    for( int i = 0; i < numOfOrd; ++i ){
        valueBest = -99999999999.9;
        for( int j = 0; j < numOfArg; ++j ){
            if( ( Arg[j] > valueBest ) && checked[j]==0){
                valueBest = Arg[j];
                indexBest = j;
            }
        }
        indexOrderd[ i ]=indexBest;
        checked[ indexBest ]=1;
    }
   // delete [] checked;
}

void ofec::eax_tsp2::TSort::index( int* Arg, int numOfArg, int* indexOrderd, int numOfOrd ){
    int indexBest = 0;
    int valueBest;
    //int *checked = new int [ numOfArg ];
    std::vector<int> checked(numOfArg);
    for( int i = 0 ; i < numOfArg ; ++i ) checked[ i ] = 0;
    for( int i = 0; i < numOfOrd; ++i ){
        valueBest = 99999999;
        for( int j = 0; j < numOfArg; ++j ){
            if( ( Arg[j] < valueBest ) && checked[j]==0){
                valueBest = Arg[j];
                indexBest = j;
            }
        }
        indexOrderd[ i ]=indexBest;
        checked[ indexBest ]=1;
    }
   // delete [] checked;
}

void ofec::eax_tsp2::TSort::indexB(const std::vector<int>& Arg, int numOfArg, std::vector<int>& indexOrderd, int numOfOrd ){
    int indexBest = 0;
    int valueBest;
    //int *checked = new int [ numOfArg ];
    std::vector<int> checked(numOfArg);
    for( int i = 0 ; i < numOfArg ; ++i ) checked[ i ] = 0;
    for( int i = 0; i < numOfOrd; ++i ){
        valueBest = -999999999;
        for( int j = 0; j < numOfArg; ++j ){
            if( ( Arg[j] > valueBest ) && checked[j]==0){
                valueBest = Arg[j];
                indexBest = j;
            }
        }
        indexOrderd[ i ]=indexBest;
        checked[ indexBest ]=1;
    }
  //  delete [] checked;
}

void ofec::eax_tsp2::TSort::sort(std::vector<int>& Arg, int numOfArg, Random* rnd){
    
   // std::sort(Arg, Arg + numOfArg);
    //selectionSort(Arg, 0, numOfArg-1);
    quickSort(Arg, 0, numOfArg-1, rnd);
}

