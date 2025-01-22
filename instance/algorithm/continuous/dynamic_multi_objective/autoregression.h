
#ifndef AUTOREGRESSION_H
#define AUTOREGRESSION_H

#include "../../../utility/linear_algebra/vector.h"
#include "../../../utility/linear_algebra/matrix.h"
#include "../../../utility/functional.h"

namespace ofec {
	class Autoregression {
	public:
		template<typename T1>
		Real k_autocorrelation(std::vector<T1> &v, int k);
		template<typename T2>
		Vector& P_autocorrelation(std::vector<T2> &v, int p);
		template<typename T3>
		std::vector<T3>& auto_regression(std::vector<T3> &v, int N, int P, std::vector<T3> &v2);
		template<typename T4>
		Real manifold_distance(std::vector<T4> &v1, std::vector<T4> &v2);
	private:

	};

	template<typename T1>
	Real Autoregression::k_autocorrelation(std::vector<T1> &v, int k) {//k�������ϵ��
		//�ȼ��������ֵ
		Real ave = 0;
		for (int i = 0; i < v.size(); i++) {
			ave += v[i];
		}
		ave = ave / v.size();
		//�ټ��������ϵ��
		Real den = 0;//��ĸ
		for (int i = 0; i < v.size(); i++) {
			den += std::pow(v[i] - ave, 2);
		}
		Real coff = 0;
		for (int i = 0; i < v.size() - k; i++) {
			coff = coff + (v[i] - ave)*(v[i + k] - ave) / den;
		}
		return coff;
	}

	//p��ƫ�����ϵ��
	template<typename T2>
	Vector& Autoregression::P_autocorrelation(std::vector<T2> &v, int p) {
		Vector coff;
		//�������е�ǰp���ͺ�������ϵ��
		Matrix rho1(p, p);
		Matrix rho2(p, p);
		Vector a(p);
		for (int i = 0; i < p; i++) {
			for (int j = 0; j < p; j++) {
				if (i < j)
					rho1[i][j] = k_autocorrelation(v, j - i);
				else if (i == j)
					rho1[i][j] = 1;
				else
					rho1[i][j] = k_autocorrelation(v, i - j);
			}
			a[i] = k_autocorrelation(v, i + 1);
		}
		//Ȼ����k�����ϵ����ɵľ��󷽳����k��ƫ���ϵ��
		rho2 = rho1.transpose();
		for (int i = 0; i < p; i++) {
			Matrix temp = rho2;
			temp[i] = a;
			temp.transpose();
			coff.pushBack(temp.determinant() / rho1.determinant());
		}
		return coff;
	}

	//P�׵������Իع�ϵ��������Ԥ��ֵ
	template<typename T3>
	std::vector<T3>& Autoregression::auto_regression(std::vector<T3> &v, int N, int P, std::vector<T3> &v2) {//pΪ������NΪʱ�����г���
		//�״θ���Ϊ3����ͨ��ʱ�����е������ϵ����ƫ�����ϵ����ȡ
		//Y=X*a+c,cΪ��ֵΪ0����̬�ֲ�
		//int P = 3;
		size_t num = v.size();
		if (num < N)
			throw "the length of the time series is not enough";
		Matrix Y(N - P, 1);
		Matrix x(N - P, P);
		//��x��Y��ֵ
		for (int i = 0; i < N - P; i++) {
			for (int j = 0; j < P; j++) {
				x[i][j] = v[i + P - j];
			}
			Y[i][0] = v[i + P];
		}
		Matrix temp1 = x.transpose();
		Matrix temp2 = temp1 * x;
		//std::cout << temp2.determinant() << std::endl;
		temp2.inverse();
		temp1 = temp2 * temp1*Y;
		Vector coff = temp1.transpose().data()[0];
		Real predict_value;
		Vector series(P);
		for (int i = 0; i < P; i++) {
			series[i] = v[num - P + i];
		}
		predict_value = series * coff;
		coff.pushBack(predict_value);
		v2 = coff.vect();
		return v2;
	}

	template<typename T4>
	Real Autoregression::manifold_distance(std::vector<T4> &v1, std::vector<T4> &v2) {
		Real distance = 0;
		for (auto &i : v1) {
			Real min_d = std::numeric_limits<Real>::max();
			for (size_t j = 0; j < v2.size(); j++) {
				Real d = euclideanDistance(v2[j].begin(), v2[j].end(), i.begin());
				if (d<min_d)
					min_d = d;
			}
			distance += min_d;
		}
		return distance / v1.size();
	}
}

#endif // !OFEC_AUTOREGRESSION_H