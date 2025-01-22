/********* Begin Register Information **********
{
    "dependency on libraries": [ "Eigen" ]
}
*********** End Register Information **********/

#ifndef OFEC_UTILS_H
#define OFEC_UTILS_H

#include <cmath>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>



namespace ofec {
    double reciprocalScalar(double x);
    double sigmoidScalar(double x);
    double logScalar(double x);
    double expScalar(double x);
    double sqrtScalar(double x);
    Eigen::MatrixXd bsxfunMinus(Eigen::MatrixXd& m, Eigen::MatrixXd& x);
    Eigen::MatrixXd bsxfunRDivide(Eigen::MatrixXd& m, Eigen::MatrixXd& x);
    Eigen::MatrixXd bsxfunPlus(Eigen::MatrixXd& m, Eigen::MatrixXd& x);
    Eigen::MatrixXd sigmoid(Eigen::MatrixXd& z);
    Eigen::MatrixXd sigmoidGradient(Eigen::MatrixXd& z);
    Eigen::MatrixXd binaryCols(Eigen::MatrixXi& labels, int numOfClasses);
    Eigen::MatrixXd expMat(Eigen::MatrixXd& z);
    Eigen::MatrixXd logMat(Eigen::MatrixXd& z);
    Eigen::MatrixXd sqrtMat(Eigen::MatrixXd& z);
    Eigen::MatrixXd powMat(Eigen::MatrixXd& z, int power);
    Eigen::MatrixXd reciprocal(Eigen::MatrixXd& z);
    double calcAccurancy(Eigen::MatrixXi& pred, Eigen::MatrixXi& labels);
}


#endif // OFEC_UTILS_H
