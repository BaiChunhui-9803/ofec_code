#include"peak.h"
#include "../../../../../core/definition.h"

namespace ofec {

    Peak::Peak(int numvar, double height, double shape, double radius, std::vector<double> position, bool rotation, std::string peakShape)
    {
        gen = std::mt19937(rd());
        uni_domain = std::uniform_real_distribution<double>(M_domain_min, M_domain_max);//M_domain_min, M_domain_max);//用于生成位置
        nor_domain = std::normal_distribution<double>(M_domain_min, M_domain_max);//M_domain_min, M_domain_max);
        uni_rand_value = std::uniform_real_distribution<double>(M_scaledDiag_min, M_scaledDiag_max);

        m_numVar = numvar;
        m_height = height;
        m_shape = shape;//由外部给的--与shape height correlation 相关
        m_radius = radius;//由外部给的
        m_rotation = rotation;
        m_peakShape = peakShape;


        if (m_position.empty()) {
            for (int i = 0; i < m_numVar; i++)
                m_position.push_back(uni_domain(gen));
        }

        m_rotationmatrix = Eigen::MatrixXd::Identity(m_numVar, m_numVar);
        if (m_rotation)//generate random rotation matrix 原文中没有纯随机，是根据-某R名人士-提出的方法
        {
            double quarterPi = OFEC_PI / 4.0;
            for (int j = 0; j < m_numVar - 1; j++) {
                for (int k = j + 1; k < m_numVar; k++) {
                    Eigen::MatrixXd r = Eigen::MatrixXd::Zero(m_numVar, m_numVar);
                    double alpha = uni_domain(gen) * (-quarterPi) + quarterPi;
                    r(j, j) = cos(alpha);
                    r(j, k) = sin(alpha);
                    r(k, j) = -sin(alpha);
                    r(k, k) = cos(alpha);
                    for (int i = 0; i < m_numVar; i++) {
                        for (int l = 0; l < m_numVar; l++) {
                            m_rotationmatrix(i, l) += r(i, l);
                        }
                    }
                }
            }
        }

        //a vector of random values-- uniform distribution
        Eigen::VectorXd v = Eigen::VectorXd::Zero(m_numVar, 1);
        for (size_t i = 0; i < v.size(); i++)
        {
            v[i] = uni_rand_value(gen);
        }



        if (m_peakShape == "sphere")
        {
            scaledDiagValues = Eigen::VectorXd::Zero(m_numVar, 1).unaryExpr([v](double dummy) {return v[0]; });//所以要用捕获列表
        }
        else if (m_peakShape == "ellipse")
        {
            scaledDiagValues = v;
        }
        else
        {
            throw "undefined shape";
        }
        Eigen::DiagonalMatrix<double, Eigen::Dynamic> dia(scaledDiagValues);
        m_covariance = m_rotationmatrix.transpose() * dia * m_rotationmatrix;
        m_covariance = m_covariance.inverse();
    }

}