
#ifndef EAX_TSP_TEVAL_H
#define EAX_TSP_TEVAL_H


#include "indi.h"
#include <fstream>
#include <string>
#include "../../../../../core/environment/environment.h"


namespace ofec {
    class Problem;
    namespace eax_tsp {


        class TEvaluator
        {
        public:
            TEvaluator();
            ~TEvaluator();

            void setInstance(Environment *env);
            //  void setInstance(const std::string& filename);/* Set the instance */
            void doIt(TIndi& indi);           /* Set the value of indi.fEvaluationValue */
          //  void writeTo(FILE* fp, TIndi& indi);     /* Write an tour to a file*/
         //   void writeToStdout(TIndi& indi); /* Write a tour to stdout */
         //   void writeTo(std::ofstream& out, const TIndi& indi);
            void transferSol(const TIndi& indi, std::vector<int>& sol);
            bool checkValid(int* array, double value);  /* Check an tour */

            int fNearNumMax; /* Maximum number of k (see below) */
            std::vector<std::vector<int>> fNearCity; /* NearCity[i][k]: k-th nearest city from */
            std::vector<std::vector<double>> fEdgeDis;  /* EdgeDis[i][j]: distance between i and j */
            int Ncity;       /* Number of cities */
        //    double* x;       /* x[i]: x-coordinate of */
        //    double* y;       /* y[i]: x-coordinate of */
         //   int* Array;
        };

    }
}

#endif
