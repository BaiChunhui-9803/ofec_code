

/*
*
*  from: https://github.com/OCEChain/denosing_autoencoder
*
*/

#ifndef OFEC_STACKEDAE_H
#define OFEC_STACKEDAE_H

#include "DAE.h"
#include "utils.h"
#include "softMax.h"

namespace ofec {
    class StackedAE {
    private:
        Eigen::MatrixXd aeTheta1;
        Eigen::MatrixXd aeTheta2;
        Eigen::MatrixXd aeB1;
        Eigen::MatrixXd aeB2;
        Eigen::MatrixXd softMaxTheta;
        //momentum;
        Eigen::MatrixXd V_theta1;
        Eigen::MatrixXd V_theta2;
        Eigen::MatrixXd V_B1;
        Eigen::MatrixXd V_B2;
        Eigen::MatrixXd V_softMaxtheta;
        //Historic Gradients
        Eigen::MatrixXd G_theta1;
        Eigen::MatrixXd G_theta2;
        Eigen::MatrixXd G_B1;
        Eigen::MatrixXd G_B2;
        Eigen::MatrixXd G_softMaxtheta;
        //parameters
        int numClasses;
        int ae1HiddenSize;
        int ae2HiddenSize;
        int inputSize;
        Eigen::MatrixXd softmaxGradient(Eigen::MatrixXd& x);
        Eigen::MatrixXd feedForward(Eigen::MatrixXd& theta, Eigen::MatrixXd& b, Eigen::MatrixXd data);
        void SGD_updateParameters(Eigen::MatrixXd& theta1Grad, Eigen::MatrixXd& theta2Grad, Eigen::MatrixXd& b1Grad, Eigen::MatrixXd& b2Grad, Eigen::MatrixXd& softmaxTheta, double alpha);
        double computeCost(Eigen::MatrixXd& theta1Grad, Eigen::MatrixXd& b1Grad, Eigen::MatrixXd& theta2Grad, Eigen::MatrixXd& b2Grad, Eigen::MatrixXd& softmaxThetaGrad, Eigen::MatrixXd& data, Eigen::MatrixXi& labels, double lambda);
    public:
        StackedAE(int ae1HiddenSize, int ae2HiddenSize, int numClasses);
        Eigen::MatrixXi predict(Eigen::MatrixXd& data);
        void fineTune(Eigen::MatrixXd& data, Eigen::MatrixXi& labels, double lambda, double alpha, int maxIter, int batchSize);
        void preTrain(Eigen::MatrixXd& data, Eigen::MatrixXi& labels, double lambda[], double alpha[], int miniBatchSize[], int maxIter[], double noiseRatio[] = NULL, double beta[] = NULL);
        Eigen::MatrixXd getAe1Theta();
        Eigen::MatrixXd getAe2Theta();
    };
}



#endif // OFEC_STACKEDAE_H
