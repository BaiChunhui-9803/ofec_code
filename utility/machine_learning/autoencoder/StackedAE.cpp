#include "stackedAE.h"

namespace ofec {
    StackedAE::StackedAE(int ae1HiddenSize, int ae2HiddenSize, int numClasses)
    {
        this->ae1HiddenSize = ae1HiddenSize;
        this->ae2HiddenSize = ae2HiddenSize;
        this->numClasses = numClasses;
    }

    Eigen::MatrixXd StackedAE::getAe1Theta()
    {
        return aeTheta1;
    }

    Eigen::MatrixXd StackedAE::getAe2Theta()
    {
        return aeTheta2;
    }

    //forward calculation
    Eigen::MatrixXd StackedAE::feedForward(Eigen::MatrixXd& theta, Eigen::MatrixXd& b,
        Eigen::MatrixXd data)
    {
        int m = data.cols();
        Eigen::MatrixXd z2 = theta * data + b.replicate(1, m);
        Eigen::MatrixXd a2 = sigmoid(z2);
        return a2;
    }

    //predict
    Eigen::MatrixXi StackedAE::predict(
        Eigen::MatrixXd& data)
    {
        //forward calculation
        char str[200] = { 0 };
        /*sprintf_s(str,"aetheta1 %d %d; %d %d",aeTheta1.rows(),aeTheta1.cols(),data.rows(),data.cols());
        MessageBoxA(NULL,str,"",MB_OK);*/
        Eigen::MatrixXd term1 = aeTheta1 * data;
        Eigen::MatrixXd z2 = bsxfunPlus(term1, aeB1);
        Eigen::MatrixXd a2 = sigmoid(z2);
        Eigen::MatrixXd term2 = aeTheta2 * a2;
        /*sprintf_s(str,"aetheta2 %d %d; %d %d",aeTheta2.rows(),aeTheta2.cols(),a2.rows(),a2.cols());
        MessageBoxA(NULL,str,"",MB_OK);*/
        Eigen::MatrixXd z3 = bsxfunPlus(term2, aeB2);
        Eigen::MatrixXd a3 = sigmoid(z3);
        Eigen::MatrixXd z4 = softMaxTheta * a3;
        /*sprintf_s(str,"softmaxTheta %d %d; %d %d",softMaxTheta.rows(),softMaxTheta.cols(),a3.rows(),a3.cols());
        MessageBoxA(NULL,str,"",MB_OK);*/
        //char str[200];
        /*sprintf_s(str,"%d %d",z4.rows(),z4.cols());
        MessageBoxA(NULL,str,"",MB_OK);*/
        Eigen::MatrixXi pred(z4.cols(), 1);
        for (int i = 0; i < z4.cols(); i++)
        {
            double max = INT_MIN;
            int idx = 0;
            for (int j = 0; j < z4.rows(); j++)
            {
                if (z4(j, i) > max)
                {
                    idx = j;
                    max = z4(j, i);
                }
            }
            pred(i, 0) = idx;
        }
        return pred;
    }

    //component wise softmax gradient
    Eigen::MatrixXd StackedAE::softmaxGradient(Eigen::MatrixXd& x)
    {
        Eigen::MatrixXd negX = x * (-1);
        Eigen::MatrixXd expX = expMat(negX);
        Eigen::MatrixXd term1 = (Eigen::MatrixXd::Ones(expX.rows(), expX.cols())
            + expX).array().square();
        Eigen::MatrixXd grad = expX.cwiseQuotient(term1);
        return grad;
    }

    //update all parameters
    void StackedAE::SGD_updateParameters(Eigen::MatrixXd& theta1Grad, Eigen::MatrixXd& theta2Grad,
        Eigen::MatrixXd& b1Grad, Eigen::MatrixXd& b2Grad,
        Eigen::MatrixXd& softmaxThetaGrad, double alpha)
    {
        V_theta1 = 0.9 * V_theta1 - alpha * theta1Grad;
        V_theta2 = 0.9 * V_theta2 - alpha * theta2Grad;
        V_B1 = 0.9 * V_B1 - alpha * b1Grad;
        V_B2 = 0.9 * V_B2 - alpha * b2Grad;
        V_softMaxtheta = 0.9 * V_softMaxtheta - alpha * softmaxThetaGrad;
        aeTheta1 += V_theta1;
        aeTheta2 += V_theta2;
        aeB1 += V_B1;
        aeB2 += V_B2;
        softMaxTheta += V_softMaxtheta;
    }

    //cost function
    double StackedAE::computeCost(Eigen::MatrixXd& theta1Grad,
        Eigen::MatrixXd& b1Grad, Eigen::MatrixXd& theta2Grad,
        Eigen::MatrixXd& b2Grad, Eigen::MatrixXd& softmaxThetaGrad,
        Eigen::MatrixXd& data, Eigen::MatrixXi& labels, double lambda)
    {
        Eigen::MatrixXd groundTruth = binaryCols(labels, numClasses);
        int M = labels.rows();
        //forward calculate
        Eigen::MatrixXd term1 = aeTheta1 * data;
        Eigen::MatrixXd z2 = bsxfunPlus(term1, aeB1);
        Eigen::MatrixXd a2 = sigmoid(z2);
        Eigen::MatrixXd term2 = aeTheta2 * a2;
        Eigen::MatrixXd z3 = bsxfunPlus(term2, aeB2);
        Eigen::MatrixXd a3 = sigmoid(z3);
        Eigen::MatrixXd z4 = softMaxTheta * a3;
        Eigen::MatrixXd a4 = expMat(z4);
        Eigen::MatrixXd a4ColSum = a4.colwise().sum();
        a4 = bsxfunRDivide(a4, a4ColSum);
        //calculate delta
        Eigen::MatrixXd delta4 = a4 - groundTruth;
        Eigen::MatrixXd delta3 = (softMaxTheta.transpose() * delta4).cwiseProduct(sigmoidGradient(z3));
        Eigen::MatrixXd delta2 = (aeTheta2.transpose() * delta3).cwiseProduct(sigmoidGradient(z2));

        //calculate delta
        softmaxThetaGrad = (groundTruth - a4) * a3.transpose() * (-1.0 / M) + softMaxTheta * lambda;

        theta2Grad = delta3 * a2.transpose() * (1.0 / M) + aeTheta2 * lambda;
        b2Grad = delta3.rowwise().sum() * (1.0 / M);
        theta1Grad = delta2 * data.transpose() * (1.0 / M) + aeTheta1 * lambda;
        b1Grad = delta2.rowwise().sum() * (1.0 / M);

        //compute cost
        double cost = (-1.0 / M) * (groundTruth.cwiseProduct(logMat(a4))).array().sum()
            + lambda / 2.0 * softMaxTheta.array().square().sum()
            + lambda / 2.0 * aeTheta1.array().square().sum()
            + lambda / 2.0 * aeTheta2.array().square().sum();

        return cost;
    }

    //fine tune the model
    void StackedAE::fineTune(Eigen::MatrixXd& data, Eigen::MatrixXi& labels,
        double lambda, double alpha, int maxIter, int batchSize)
    {
        Eigen::MatrixXd theta1Grad(aeTheta1.rows(), aeTheta1.cols());
        Eigen::MatrixXd theta2Grad(aeTheta2.rows(), aeTheta2.cols());
        Eigen::MatrixXd b1Grad(aeB1.rows(), aeB1.cols());
        Eigen::MatrixXd b2Grad(aeB2.rows(), aeB2.cols());
        Eigen::MatrixXd softmaxThetaGrad(softMaxTheta.rows(), softMaxTheta.cols());
        Eigen::MatrixXd miniTrainData(data.rows(), batchSize);
        Eigen::MatrixXi miniLabels(batchSize, 1);
        int iter = 1;
        int numBatches = data.cols() / batchSize;

        //mini batch stochastic gradient decent
        for (int i = 0; i < maxIter; i++)
        {
            double J = 0;
            // compute the cost
            for (int j = 0; j < numBatches; j++)
            {
                miniTrainData = data.middleCols(j * batchSize, batchSize);
                miniLabels = labels.middleRows(j * batchSize, batchSize);
                J += computeCost(theta1Grad, b1Grad, theta2Grad,
                    b2Grad, softmaxThetaGrad, miniTrainData, miniLabels, lambda);
                if (miniTrainData.cols() < 1 || miniTrainData.rows() < 1)
                {
                    std::cout << "Too few training examples!" << std::endl;
                }


                if (fabs(J) < 0.001)
                {

                }
                if (j == 0) {
                    V_theta1.setZero(theta1Grad.rows(), theta1Grad.cols());
                    V_theta2.setZero(theta2Grad.rows(), theta2Grad.cols());
                    V_B1.setZero(b1Grad.rows(), b1Grad.cols());
                    V_B2.setZero(b2Grad.rows(), b2Grad.cols());
                    V_softMaxtheta.setZero(softmaxThetaGrad.rows(), softmaxThetaGrad.cols());
                }
                SGD_updateParameters(theta1Grad, theta2Grad, b1Grad, b2Grad, softmaxThetaGrad, alpha);
            }
            J = J / numBatches;
            std::cout << "iter: " << iter++ << "  cost: " << J << std::endl;
        }
    }

    //pretrain the model
    void StackedAE::preTrain(Eigen::MatrixXd& data, Eigen::MatrixXi& labels,
        double lambda[], double alpha[], int miniBatchSize[],
        int maxIter[], double noiseRatio[],
        double beta[])
    {
        int numOfExamples = data.cols();
        int ndim = data.rows();
        inputSize = ndim;
        //stacked denoising autoencoders
        std::cout << "PreTraining with denoising autoencoder ..." << std::endl;
        //train the first denoising autoencoder
        DAE ae1(ndim, ae1HiddenSize);
        std::cout << "PreTraining ae1 ..." << std::endl;
        ae1.train(data, noiseRatio[0], alpha[0], maxIter[0], miniBatchSize[0]);

        Eigen::MatrixXd theta1 = ae1.getTheta();
        aeTheta1.resize(theta1.rows(), theta1.cols());
        aeTheta1 = theta1;
        Eigen::MatrixXd b1 = ae1.getBias();
        aeB1.resize(b1.rows(), b1.cols());
        aeB1 = b1;

        //train the second denoising autoencoder
        Eigen::MatrixXd ae1Features = feedForward(aeTheta1, aeB1, data);
        DAE ae2(ae1HiddenSize, ae2HiddenSize);
        std::cout << "PreTraining ae2 ..." << std::endl;
        ae2.train(ae1Features, noiseRatio[1], alpha[1], maxIter[1], miniBatchSize[1]);

        Eigen::MatrixXd theta2 = ae2.getTheta();
        aeTheta2.resize(theta2.rows(), theta2.cols());
        aeTheta2 = theta2;
        Eigen::MatrixXd b2 = ae2.getBias();
        aeB2.resize(b2.rows(), b2.cols());
        aeB2 = b2;
        //train the softmax regression
        Eigen::MatrixXd ae2Features = feedForward(aeTheta2, aeB2, ae1Features);
        std::cout << "PreTraining softmax ..." << std::endl;
        Softmax softmax(ae2HiddenSize, numClasses);
        softmax.train(ae2Features, labels, lambda[2], alpha[2], maxIter[2], miniBatchSize[2]);
        Eigen::MatrixXd smTheta = softmax.getTheta();
        softMaxTheta.resize(smTheta.rows(), smTheta.cols());
        softMaxTheta = smTheta;
    }
}





