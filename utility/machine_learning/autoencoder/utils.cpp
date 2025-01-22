#include "utils.h"
#include<functional>
#include<algorithm>

using namespace std;

namespace ofec {
    Eigen::MatrixXd bsxfunMinus(Eigen::MatrixXd& m, Eigen::MatrixXd& x)
    {
        Eigen::MatrixXd r = m;
        if (x.rows() == 1)
        {
            r = x.replicate(m.rows(), 1);
        }
        if (x.cols() == 1)
        {
            r = x.replicate(1, m.cols());
        }
        return m - r;
    }

    Eigen::MatrixXd bsxfunRDivide(Eigen::MatrixXd& m, Eigen::MatrixXd& x)
    {
        Eigen::MatrixXd r = m;
        if (x.rows() == 1)
        {
            r = x.replicate(m.rows(), 1);
        }
        if (x.cols() == 1)
        {
            r = x.replicate(1, m.cols());
        }
        return m.cwiseQuotient(r);
    }



    Eigen::MatrixXd bsxfunPlus(Eigen::MatrixXd& m, Eigen::MatrixXd& x)
    {
        Eigen::MatrixXd r = m;
        if (x.rows() == 1)
        {
            r = x.replicate(m.rows(), 1);
        }
        if (x.cols() == 1)
        {
            r = x.replicate(1, m.cols());
        }
        return m + r;
    }

    Eigen::MatrixXd sigmoidGradient(Eigen::MatrixXd& z)
    {
        //return sigmoid(z) .* (1 - sigmoid(z))
        Eigen::MatrixXd result;
        Eigen::MatrixXd sigm = sigmoid(z);
        Eigen::MatrixXd item = Eigen::MatrixXd::Ones(z.rows(), z.cols()) - sigm;
        result = sigm.cwiseProduct(item);
        return result;
    }

    

    Eigen::MatrixXd binaryCols(Eigen::MatrixXi& labels, int numOfClasses)
    {
        // return binary code of labels
        //eye function
        Eigen::MatrixXd e = Eigen::MatrixXd::Identity(numOfClasses, numOfClasses);
        int numOfExamples = labels.rows();
        int inputSize = e.cols();
        Eigen::MatrixXd result(inputSize, numOfExamples);
        for (int i = 0; i < numOfExamples; i++)
        {
            int idx = labels(i, 0);
            result.col(i) = e.col(idx);
        }
        return result;
    }


    double calcAccurancy(Eigen::MatrixXi& pred, Eigen::MatrixXi& labels)
    {
        int numOfExamples = pred.rows();
        double sum = 0;
        for (int i = 0; i < numOfExamples; i++)
        {
            if (pred(i, 0) == labels(i, 0))
            {
                sum += 1;
            }
        }
        return sum / numOfExamples;
    }

    

    double reciprocalScalar(double x)
    {
        return 1.0 / x;
    }

    //scalar sigmoid function
    double sigmoidScalar(double x)
    {
        return 1.0 / (1 + exp(-x));
    }

    //scalar log function
    double logScalar(double x)
    {
        return log(x);
    }

    //scalar exp function
    double expScalar(double x)
    {
        return exp(x);
    }

    //scalar sqrt function
    double sqrtScalar(double x)
    {
        return sqrt(x);
    }

    //component wise sigmoid function
    Eigen::MatrixXd sigmoid(Eigen::MatrixXd& z)
    {
        return z.unaryExpr(&sigmoidScalar);
    }

    Eigen::MatrixXd sqrtMat(Eigen::MatrixXd& z)
    {
        return z.unaryExpr(&sqrtScalar);
    }

    //return 1.0 ./ z
    Eigen::MatrixXd reciprocal(Eigen::MatrixXd& z)
    {
        return z.unaryExpr(&reciprocalScalar);
    }

    //component wise exp function
    Eigen::MatrixXd expMat(Eigen::MatrixXd& z)
    {
        return z.unaryExpr(&expScalar);
    }

    //component wise log function
    Eigen::MatrixXd logMat(Eigen::MatrixXd& z)
    {
        return z.unaryExpr(&logScalar);
    }
}


