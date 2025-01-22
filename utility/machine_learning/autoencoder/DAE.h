/********* Begin Register Information **********
{
    "dependency on libraries": [ "Eigen" ]
}
*********** End Register Information **********/
/*
*
*  from: https://github.com/OCEChain/denosing_autoencoder
*
*/

#ifndef OFEC_DAE_H
#define OFEC_DAE_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <ctime>
#include <iostream>
#include <fstream>
#include "utils.h"

namespace ofec {
    class DAE {
    public:
        Eigen::MatrixXd theta1;
        Eigen::MatrixXd theta2;
        Eigen::MatrixXd b1;
        Eigen::MatrixXd b2;
        int inputSize;
        int hiddenSize;
        DAE(int inputSize, int hiddenSize);
        void train(Eigen::MatrixXd& trainData, double noiseRatio, double alpha, int maxIter, int miniBatchSize);
        Eigen::MatrixXd getTheta();
        Eigen::MatrixXd getBias();
    private:
        Eigen::MatrixXd noiseInput(Eigen::MatrixXd& z, double noiseRatio);
        Eigen::MatrixXd randomInitialize(int lIn, int lOut);
        void updateParameters(Eigen::MatrixXd& theta1Grad1, Eigen::MatrixXd& theta2Grad2, Eigen::MatrixXd& b1Grad, Eigen::MatrixXd& b2Grad, double alpha);
        void miniBatchSGD(Eigen::MatrixXd& trainData, Eigen::MatrixXd& noiseData, double alpha, int maxIter, int batchSize);
        double computeCost(Eigen::MatrixXd& data, Eigen::MatrixXd& noiseData, Eigen::MatrixXd& theta1Grad, Eigen::MatrixXd& theta2Grad, Eigen::MatrixXd& b1Grad, Eigen::MatrixXd& b2Grad);
    };
}


#endif // OFEC_DAE_H
