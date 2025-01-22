#include "sort.h"
#include <algorithm>


//using namespace eax_tsp;

//TSort* tSort = NULL;

void ofec::eax_tsp::TSort::InitSort(){
  //  tSort = new TSort();
}

//void eax_tsp::TSort::swap(int &x, int &y){
//    int s=x;
//    x=y;
//    y=s;
//}

void ofec::eax_tsp::TSort::selectionSort(std::vector<int>& Arg, int l, int r){
    int id;
    for(int i=l;i<r;++i){
        id=i;
        for(int j=i+1;j<=r;++j)
            if(Arg[j]<Arg[id]) id=j;
        std::swap(Arg[i], Arg[id]);
        //eax_tsp::swap(Arg[i], Arg[id]);
    }
}

int ofec::eax_tsp::TSort::partition(std::vector<int>& Arg, int l, int r, Random *rnd){
    using namespace ofec;
//    int id=l+rand()%(r-l+1);
    int id = rnd->uniform.nextNonStd<int>(l, r + 1);
    std::swap(Arg[l], Arg[id]);
    id=l;
    for(int i=l+1;i<=r;++i)
        if(Arg[i]<Arg[l]) std::swap(Arg[++id], Arg[i]);
    std::swap(Arg[l], Arg[id]);
    return id;
}

void ofec::eax_tsp::TSort::quickSort(std::vector<int>& Arg, int l, int r, Random *rnd){
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

ofec::eax_tsp::TSort::TSort(const std::shared_ptr<Random> &rnd) : m_random(rnd){}
ofec::eax_tsp::TSort::~TSort(){}

void ofec::eax_tsp::TSort::index( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd ){
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

void ofec::eax_tsp::TSort::indexB(const std::vector<double>& Arg, int numOfArg, std::vector<int>& indexOrderd, int numOfOrd )
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

void ofec::eax_tsp::TSort::index( int* Arg, int numOfArg, int* indexOrderd, int numOfOrd ){
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

void ofec::eax_tsp::TSort::indexB(const std::vector<int>& Arg, int numOfArg, std::vector<int>& indexOrderd, int numOfOrd ){
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

void ofec::eax_tsp::TSort::sort(std::vector<int>& Arg, int numOfArg ){
    
   // std::sort(Arg, Arg + numOfArg);
    //selectionSort(Arg, 0, numOfArg-1);
    quickSort(Arg, 0, numOfArg-1, m_random.get());
}

