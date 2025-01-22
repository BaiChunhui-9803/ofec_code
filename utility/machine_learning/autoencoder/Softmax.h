

/*
*
*  from: https://github.com/OCEChain/denosing_autoencoder
*
*/

#ifndef OFEC_SOFTMAX_H
#define OFEC_SOFTMAX_H

#include "utils.h"
#include <ctime>
#include <fstream>
#include <iostream>

namespace ofec {
    class Softmax
    {
    private:
        Eigen::MatrixXd theta;
        int inputSize;
        int numClasses;
        Eigen::MatrixXd randomInitialize(int lIn, int lOut);
        double computeCost(double lambda, Eigen::MatrixXd& data, Eigen::MatrixXi& labels, Eigen::MatrixXd& thetaGrad);
        void miniBatchSGD(Eigen::MatrixXd& trainData, Eigen::MatrixXi& labels, double lambda, double alpha, int maxIter, int batchSize);
    public:
        Softmax(int inputSize, int numClasses);
        Eigen::MatrixXi predict(Eigen::MatrixXd& data);
        void train(Eigen::MatrixXd& data, Eigen::MatrixXi& labels, double lambda, double alpha, int maxIter, int miniBatchSize);
        Eigen::MatrixXd getTheta();
    };
}



#endif // OFEC_SOFTMAX_H
