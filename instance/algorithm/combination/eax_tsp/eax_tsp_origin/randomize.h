#ifndef EAX_TSP_RAND_H
#define EAX_TSP_RAND_H

#include <memory>
#include "../../../../../core/random/newran.h"


namespace ofec::eax_tsp {

    class TRandom {
    public:
        std::shared_ptr<Random> m_random;

        TRandom(const std::shared_ptr<Random> &rnd);
        ~TRandom();
        int Integer(int minNumber, int maxNumber);
        double Double(double minNumber, double maxNumber);
        double normalDistribution(double mu, double sigma);
        void permutation(std::vector<int>& arr, int numOfelement, int numOfSample);
        void shuffle(std::vector<int>& arr, int numOfElement);
    };
}
#endif


