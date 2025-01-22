
#ifndef OFEC_PEAK_H
#define OFEC_PEAK_H

//#define M_PI 3.14
#define M_domain_min 0.0 //ÿһάȡֵ��Χ
#define M_domain_max 1.0 
#define M_height 1.0//ȫ�����ŵķ�
//#define M_numVar 2
#define M_factor 0.8
#define M_radii 0.65
// scaledDiagValues 
#define M_scaledDiag_min 0.0025
#define M_scaledDiag_max 0.0525

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/Core>

namespace ofec {
    class Peak {
    protected:
        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<double> uni_domain;//M_domain_min, M_domain_max);//��������λ��
        std::normal_distribution<double> nor_domain;//M_domain_min, M_domain_max);
        std::uniform_real_distribution<double> uni_rand_value;
    public:
        //һ�������������
        int m_numVar;
        double m_height;
        std::vector<double> m_position;//�÷������  Ӧ��eigen vector ��һЩ
        std::string m_peakShape;//ֻ��Բ����Բ  �����Խ�������ֵ�Ƿ����
        bool m_rotation; //�÷��Ӧ����ת����rotation  �����ǶԽ���Ԫ���Ƿ�Ϊ0
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_rotationmatrix;//����ת���ǵ�λ��
        Eigen::VectorXd scaledDiagValues;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_covariance;

        //���ڼ���g��x��
        double m_radius;
        double m_shape;
        //����
        double m_distance = 0.0;//����make funnel����

    public:

        Peak(int numvar = 2, double height = 1.0, double shape = 0.5, double radius = 0.5, std::vector<double> position = {}, bool rotation = true, std::string peakShape = "ellipse");//Ĭ���βη���β��

    };

}


#endif // !OFEC_PEAK_H
