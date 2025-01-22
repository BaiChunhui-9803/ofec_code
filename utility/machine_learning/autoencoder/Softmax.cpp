#include "softmax.h"

namespace ofec {
    Softmax::Softmax(int inputSize, int numClasses)
    {
        this->inputSize = inputSize;
        this->numClasses = numClasses;
        theta = randomInitialize(numClasses, inputSize);
    }

    Eigen::MatrixXd Softmax::getTheta()
    {
        return theta;
    }

    //random initialize the weight
    Eigen::MatrixXd Softmax::randomInitialize(int lIn, int lOut)
    {
        //Random initialize the weight in a specific range
        int i, j;
        double epsilon = sqrt(6.0 / (lIn + lOut + 1));
        Eigen::MatrixXd result(lIn, lOut);
        srand((unsigned int)time(NULL));
        for (i = 0; i < lOut; i++)
        {
            for (j = 0; j < lIn; j++)
            {
                result(j, i) = ((double)rand() / (double)RAND_MAX) * 2 * epsilon - epsilon;
            }
        }
        return result;
    }

    //Predict
    Eigen::MatrixXi Softmax::predict(Eigen::MatrixXd& data)
    {
        //cout << theta.rows() << " " << theta.cols() << endl;
        //cout << data.rows() << " " << data.cols() << endl;
        Eigen::MatrixXd m = theta * data;
        Eigen::MatrixXi pred(m.cols(), 1);
        for (int i = 0; i < m.cols(); i++)
        {
            double max = 0;
            int idx = 0;
            for (int j = 0; j < m.rows(); j++)
            {
                if (m(j, i) > max)
                {
                    max = m(j, i);
                    idx = j;
                }
            }
            pred(i, 0) = idx;
        }
        return pred;
    }

    //cost function
    double Softmax::computeCost(double lambda, Eigen::MatrixXd& data,
        Eigen::MatrixXi& labels, Eigen::MatrixXd& thetaGrad)
    {
        int numCases = data.cols();
        Eigen::MatrixXd groundTruth = binaryCols(labels, numClasses);
        //
        Eigen::MatrixXd M = theta * data;
        Eigen::MatrixXd maxM = M.colwise().maxCoeff();

        M = bsxfunMinus(M, maxM);
        Eigen::MatrixXd expM = expMat(M);

        Eigen::MatrixXd tmp1 = (expM.colwise().sum()).replicate(numClasses, 1);

        Eigen::MatrixXd p = expM.cwiseQuotient(tmp1);
        //compute cost
        double cost = (groundTruth.cwiseProduct(logMat(p))).sum() * (-1.0 / numCases)
            + (lambda / 2.0) * theta.array().square().sum();
        //compute the gradient of theta
        thetaGrad = (groundTruth - p) * data.transpose() * (-1.0 / numCases)
            + theta * lambda;
        return cost;
    }

    //mini batch stochastic gradient descent
    void Softmax::miniBatchSGD(Eigen::MatrixXd& trainingData, Eigen::MatrixXi& labels, double lambda,
        double alpha, int maxIter, int batchSize)
    {
        //get the binary code of labels
        Eigen::MatrixXd thetaGrad(theta.rows(), theta.cols());
        Eigen::MatrixXd miniTrainingData(trainingData.rows(), batchSize);
        Eigen::MatrixXi miniLabels(batchSize, 1);
        int iter = 1;
        int numBatches = trainingData.cols() / batchSize;

        //test
        std::cout << "numBatches: " << numBatches << std::endl;
        //end test

        //mini batch stochastic gradient decent
        for (int i = 0; i < maxIter; i++)
        {
            double J = 0;
            // compute the cost
            for (int j = 0; j < numBatches; j++)
            {
                //test
                std::cout << i << " " << j << std::endl;
                //end test
                miniTrainingData = trainingData.middleCols(j * batchSize, batchSize);
                miniLabels = labels.middleRows(j * batchSize, batchSize);
                J += computeCost(lambda, miniTrainingData, miniLabels, thetaGrad);
                if (miniTrainingData.cols() < 1 || miniTrainingData.rows() < 1)
                {
                    std::cout << "Too few training examples!" << std::endl;
                }

                if (fabs(J) < 0.001)
                {
                    break;
                }
                //update theta
                theta -= thetaGrad * alpha;
            }
            J = J / numBatches;
            std::cout << "iter: " << iter++ << "  cost: " << J << std::endl;
        }
    }

    //train the model
    void Softmax::train(Eigen::MatrixXd& data, Eigen::MatrixXi& labels,
        double lambda, double alpha,
        int maxIter, int miniBatchSize)
    {
        miniBatchSGD(data, labels, lambda, alpha, maxIter, miniBatchSize);
    }
}




