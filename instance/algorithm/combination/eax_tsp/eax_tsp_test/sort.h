#ifndef OFEC_EAX_TSP2_TSORT_H
#define OFEC_EAX_TSP2_TSORT_H

#include <algorithm>
#include <vector>
#include <memory>
#include "../../../../../core/random/newran.h"

namespace ofec::eax_tsp2 {
    
    class TEnvironment;
    

    class TSort {
    protected:
        TEnvironment* m_env;
        
    public:
        void setEnv(TEnvironment* env);
        
        TSort();
        ~TSort();
        void index(double* Arg, int numOfArg, int* indexOrderd, int numOfOrd);
        void index(int* Arg, int numOfArg, int* indexOrderd, int numOfOrd);
        void indexB(const std::vector<double>& Arg, int numOfArg, std::vector<int>& indexOrderd, int numOfOrd);
        void indexB(const std::vector<int>& Arg, int numOfArg, std::vector<int>& indexOrderd, int numOfOrd);
        void sort(std::vector<int>& Arg, int numOfArg, Random *rnd);

        //void swap(int& x, int& y);
        void selectionSort(std::vector<int>& Arg, int l, int r);
        int partition(std::vector<int>& Arg, int l, int r, Random *rnd);
        void quickSort(std::vector<int>& Arg, int l, int r, Random *rnd);

    };

 //   extern TSort* tSort;

}

#endif // !OFEC_EAX_TSP_TSORT_H
