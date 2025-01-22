#ifndef EAX_TSP_MULTITHREAD_RAND_H
#define EAX_TSP_MULTITHREAD_RAND_H

#include <memory>
#include "../../../../../core/random/newran.h"


namespace ofec::eax_tsp_mt {

    class TRandom {
    public:


        int Integer(int minNumber, int maxNumber, Random* rnd);
        double Double(double minNumber, double maxNumber, Random* rnd);
        double normalDistribution(double mu, double sigma, Random* rnd);
        void permutation(std::vector<int>& arr, int numOfelement, int numOfSample, Random* rnd);
        void shuffle(std::vector<int>& arr, int numOfElement, Random* rnd);
    };
}
#endif


