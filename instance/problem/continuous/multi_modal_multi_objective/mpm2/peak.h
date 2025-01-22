
#ifndef OFEC_PEAK_H
#define OFEC_PEAK_H

//#define M_PI 3.14
#define M_domain_min 0.0 //每一维取值范围
#define M_domain_max 1.0 
#define M_height 1.0//全局最优的峰
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
        std::uniform_real_distribution<double> uni_domain;//M_domain_min, M_domain_max);//用于生成位置
        std::normal_distribution<double> nor_domain;//M_domain_min, M_domain_max);
        std::uniform_real_distribution<double> uni_rand_value;
    public:
        //一个峰的自身属性
        int m_numVar;
        double m_height;
        std::vector<double> m_position;//该峰的中心  应该eigen vector 好一些
        std::string m_peakShape;//只有圆或椭圆  决定对角线特征值是否相等
        bool m_rotation; //该峰对应的旋转矩阵rotation  决定非对角线元素是否为0
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_rotationmatrix;//不旋转就是单位阵
        Eigen::VectorXd scaledDiagValues;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m_covariance;

        //用于计算g（x）
        double m_radius;
        double m_shape;
        //其他
        double m_distance = 0.0;//用于make funnel距离

    public:

        Peak(int numvar = 2, double height = 1.0, double shape = 0.5, double radius = 0.5, std::vector<double> position = {}, bool rotation = true, std::string peakShape = "ellipse");//默认形参放置尾部

    };

}


#endif // !OFEC_PEAK_H
