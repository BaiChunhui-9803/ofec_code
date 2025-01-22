#ifndef OFEC_PREDICT_PF_H
#define OFEC_PREDICT_PF_H

#include "../../../../../utility/linear_algebra/matrix.h"
#include "../../../../../utility/functional.h"
#include "../../../template/classic/de/individual.h"
#include <cmath>

namespace ofec {
	//template <typename Solution = OFEC::DE::Solution>
	//class Kalman_Filter {
	//protected:
	//	Matrix m_A;
	//	Matrix m_P;
	//	Matrix m_Q;
	//	Matrix m_R;
	//	Matrix m_H;
	//	Matrix m_predict_P;//
	//	Matrix m_K;//kalman gain
	//	size_t m_index;//indicate where have a change
	//	std::vector<std::vector<Matrix>> m_pre_estimate;//store previous estimate state
	//	std::vector<std::vector<Matrix>> m_post_estimate;//store posterior estimate state
	//public:
	//	Kalman_Filter(){ }
	//	void initialize_para(const std::vector<std::vector<Solution>>& v);
	//	void update_para(const std::vector<std::vector<Solution>>& v);
	//	void predict_pre_estimate(const std::vector<std::vector<Solution>>& v);
	//	void predict(const std::vector<std::vector<Solution>>& v, std::vector<std::vector<Real>>& prediction_result);
	//	void set_A(Matrix& m) { m_A = m; }
	//	void set_P(Matrix& m) { m_P = m; }
	//	void set_Q(Matrix& m) { m_Q = m; }
	//	void set_R(Matrix& m) { m_R = m; }
	//	void set_H(Matrix& m) { m_H = m; }
	//	void set_K(Matrix& m) { m_K = m; }
	//	void set_predict_P(Matrix& m) { m_predict_P = m; }
	//};
	//template <typename Solution>
	//void Kalman_Filter<Solution>::initialize_para(const std::vector<std::vector<Solution>>& v) {
	//	Matrix K_temp(2, 2);
	//	K_temp[0][0] = 1; K_temp[0][1] = 0; K_temp[1][0] = 0; K_temp[1][1] = 1;
	//	set_K(K_temp);
	//	Matrix A_temp(2, 2);
	//	A_temp[0][0] = 1; A_temp[0][1] = 1; A_temp[1][0] = 0; A_temp[1][1] = 1;
	//	set_A(A_temp);
	//	Matrix P_temp(2, 2);
	//	P_temp[0][0] = 1; P_temp[0][1] = 0; P_temp[1][0] = 0; P_temp[1][1] = 1;
	//	set_P(P_temp);
	//	Matrix Q_temp(2, 2);
	//	Q_temp[0][0] = 0.04; Q_temp[0][1] = 0; Q_temp[1][0] = 0; Q_temp[1][1] = 0.04;
	//	set_Q(Q_temp);
	//	Matrix R_temp(2, 2);
	//	R_temp[0][0] = 0.01; R_temp[0][1] = 0; R_temp[1][0] = 0; R_temp[1][1] = 0.01;
	//	set_R(R_temp);
	//	set_H(P_temp);
	//	//initialize m_post_estimate and m_pre_estimate
	//	if (v.size() == 0 || v[0].size() < 3)
	//		std::cout << "prediction data too little" << std::endl;
	//	else {
	//		//allocate memory
	//		m_post_estimate.resize(v.size());
	//		m_pre_estimate.resize(v.size());
	//		for (size_t seq = 0; seq < v.size(); ++seq) {
	//			m_post_estimate[seq].resize(v[0][0].objectiveSize());
	//			m_pre_estimate[seq].resize(v[0][0].objectiveSize());
	//			Matrix temp(2, 1);
	//			for (int j = 0; j < v[0][0].objectiveSize(); ++j) {
	//				m_post_estimate[seq][j] = temp;
	//				m_pre_estimate[seq][j] = temp;
	//			}
	//		}
	//		//initialize m_post_estimate
	//		for (size_t seq = 0; seq < v.size(); ++seq) {
	//			for (size_t j = 0; j < v[0][0].objectiveSize(); ++j) {
	//				m_post_estimate[seq][j][0][0] = v[seq][1].objective(j);
	//				m_post_estimate[seq][j][1][0] = v[seq][1].objective(j) - v[seq][0].objective(j);
	//			}
	//		}
	//		/*for (size_t seq = 0; seq < v.size(); ++seq) {
	//			std::vector<Matrix> temp(v[0][0].objectiveSize());
	//			for (size_t j = 0; j < v[0][0].objectiveSize(); ++j) {
	//				Matrix state(2, 1);
	//				state[0][0] = v[seq][1].objective()[j];
	//				state[1][0] = v[seq][1].objective()[j] - v[seq][0].objective()[j];
	//				temp[j] = state;					
	//			}
	//			m_post_estimate[seq] = temp;
	//		}*/
	//	}		
	//}
	//template <typename Solution>
	//void Kalman_Filter<Solution>::update_para(const std::vector<std::vector<Solution>> &v) {
	//	m_predict_P = m_A * m_P * m_A.transpose() + m_Q;
	//	Matrix S = m_H * m_predict_P * m_H.transpose() + m_R;
	//	S.inverse();		
	//	m_K = m_predict_P * m_H.transpose() * S;
	//	m_P = m_predict_P - m_K * m_H * m_predict_P;
	//	//update posterior estimate state
	//	for (size_t seq = 0; seq < v.size(); ++seq) {
	//		for (size_t j = 0; j < v[0][0].objectiveSize(); j++) {				
	//			Matrix observe_Matrix(2, 1);
	//			observe_Matrix[0][0] = v[seq][2].objective(j);
	//			observe_Matrix[1][0] = v[seq][2].objective(j) - v[seq][1].objective(j);				
	//			m_post_estimate[seq][j] = m_pre_estimate[seq][j] + m_K * (observe_Matrix - m_H * m_pre_estimate[seq][j]);
	//		}
	//	}
	//	/*for (size_t seq = 0; seq < v.size(); seq++) {
	//		std::vector<Matrix> temp_Matrix(v[0][0].objectiveSize());
	//		for (size_t j = 0; j < v[0][0].objectiveSize(); j++) {
	//			Matrix state;
	//			Matrix observe_Matrix(2, 1);
	//			observe_Matrix[0][0] = v[seq][2].objective(j);
	//			observe_Matrix[1][0] = v[seq][2].objective(j) - v[seq][1].objective(j);
	//			state = m_pre_estimate[seq][j] + m_K * (observe_Matrix - m_H * m_pre_estimate[seq][j]);
	//			temp_Matrix[j] = state;
	//		}
	//		m_post_estimate[seq] = temp_Matrix;
	//	}*/
	//	/*for (size_t seq = 0; seq < v.size(); seq++) {
	//		std::vector<Matrix> temp_Matrix;
	//		for (size_t j = 0; j < v[0][0].objectiveSize(); j++) {
	//			Matrix temp;
	//			Matrix observe_Matrix(2, 1);
	//			observe_Matrix[0][0] = v[seq][2].objective(j);
	//			observe_Matrix[1][0] = v[seq][2].objective(j) - v[seq][1].objective(j);
	//			temp = m_pre_estimate[seq][j] + m_K * (observe_Matrix - m_H * m_pre_estimate[seq][j]);
	//			temp_Matrix.emplace_back(temp);
	//		}
	//		m_post_estimate[seq] = temp_Matrix;
	//	}*/
	//}
	//template <typename Solution>
	//void Kalman_Filter<Solution>::predict_pre_estimate(const std::vector<std::vector<Solution>> &v) {
	//	for (size_t seq = 0; seq < v.size(); ++seq) {
	//		for (size_t j = 0; j < v[0][0].objectiveSize(); ++j) {
	//			m_pre_estimate[seq][j] = m_A * m_post_estimate[seq][j];
	//		}
	//	}
	//}
	//template <typename Solution>
	//void Kalman_Filter<Solution>::predict(const std::vector<std::vector<Solution>>& v, std::vector<std::vector<Real>>& prediction_result) {
	//	initialize_para(v);
	//	predict_pre_estimate(v);
	//	update_para(v);
	//	predict_pre_estimate(v);
	//	if (prediction_result.size() == 0) {
	//		prediction_result.resize(v.size());
	//		for (size_t j = 0; j < v.size(); ++j)
	//			prediction_result[j].resize(v[0][0].objectiveSize());
	//	}
	//	for (size_t seq = 0; seq < v.size(); ++seq) {
	//		for (size_t j = 0; j < v[0][0].objectiveSize(); ++j) {
	//			prediction_result[seq][j] = m_pre_estimate[seq][j][0][0];
	//		}			
	//	}
	//}

	template <typename Solution = ofec::IndDE>
	class kalmanFilter {
	protected:
		Matrix m_A;
		Matrix m_P;
		Matrix m_Q;
		Matrix m_R;
		Matrix m_H;
		Matrix m_predict_P;//
		Matrix m_K;//kalman gain
		size_t m_index;//indicate where have a change
		std::vector<Matrix> m_pre_estimate;//store previous estimate state
		std::vector<Matrix> m_post_estimate;//store posterior estimate state
		std::vector<Real> m_predict_result;//store predict result, i.e., m_pre_estimate
		bool m_type = false;
	public:
		//Kalman_Filter() { }
		void initialize(const std::vector<Solution>& v);
		void updatePara(const std::vector<Solution>& v);
		const std::vector<Real>& predict(const std::vector<Solution>& v);

		//const std::vector<Real>& predict_result() const { return m_predict_result; }
		bool type() { return m_type; }

		void set_A(Matrix& m) { m_A = m; }
		void set_P(Matrix& m) { m_P = m; }
		void set_Q(Matrix& m) { m_Q = m; }
		void set_R(Matrix& m) { m_R = m; }
		void set_H(Matrix& m) { m_H = m; }
		void set_K(Matrix& m) { m_K = m; }
		void set_predict_P(Matrix& m) { m_predict_P = m; }
	};

	template <typename Solution>
	void kalmanFilter<Solution>::initialize(const std::vector<Solution>& v) {
		Matrix K_temp(2, 2);
		K_temp[0][0] = 1; K_temp[0][1] = 0; K_temp[1][0] = 0; K_temp[1][1] = 1;
		set_K(K_temp);
		Matrix A_temp(2, 2);
		A_temp[0][0] = 1; A_temp[0][1] = 1; A_temp[1][0] = 0; A_temp[1][1] = 1;
		set_A(A_temp);
		Matrix P_temp(2, 2);
		P_temp[0][0] = 1; P_temp[0][1] = 0; P_temp[1][0] = 0; P_temp[1][1] = 1;
		set_P(P_temp);
		Matrix Q_temp(2, 2);
		Q_temp[0][0] = 0.04; Q_temp[0][1] = 0; Q_temp[1][0] = 0; Q_temp[1][1] = 0.04;
		set_Q(Q_temp);
		Matrix R_temp(2, 2);
		R_temp[0][0] = 0.01; R_temp[0][1] = 0; R_temp[1][0] = 0; R_temp[1][1] = 0.01;
		set_R(R_temp);
		set_H(P_temp);
		//initialize m_post_estimate and m_pre_estimate
		if (v.size() < 2)
			std::cout << "prediction data too little" << std::endl;
		else {
			//allocate memory	
			m_predict_result.resize(v[0].objectiveSize());
			Matrix temp(2, 1);
			m_post_estimate.resize(v[0].objectiveSize(), temp);
			m_pre_estimate.resize(v[0].objectiveSize(), temp);
			//initialize m_post_estimate			
			for (size_t j = 0; j < v[0].objectiveSize(); ++j) {
				m_post_estimate[j][0][0] = v[1].objective(j);
				m_post_estimate[j][1][0] = v[1].objective(j) - v[0].objective(j);
			}
			/*for (size_t seq = 0; seq < v.size(); ++seq) {
				std::vector<Matrix> temp(v[0][0].objectiveSize());
				for (size_t j = 0; j < v[0][0].objectiveSize(); ++j) {
					Matrix state(2, 1);
					state[0][0] = v[seq][1].objective()[j];
					state[1][0] = v[seq][1].objective()[j] - v[seq][0].objective()[j];
					temp[j] = state;
				}
				m_post_estimate[seq] = temp;
			}*/
		}
	}
	template <typename Solution>
	void kalmanFilter<Solution>::updatePara(const std::vector<Solution>& v) {
		m_predict_P = m_A * m_P * m_A.transpose() + m_Q;
		Matrix S = m_H * m_predict_P * m_H.transpose() + m_R;
		S.inverse();
		m_K = m_predict_P * m_H.transpose() * S;
		m_P = m_predict_P - m_K * m_H * m_predict_P;
		//update posterior estimate state		
		for (size_t j = 0; j < v[0].objectiveSize(); j++) {
			Matrix observe_Matrix(2, 1);
			observe_Matrix[0][0] = v[1].objective(j);	//observe value f
			observe_Matrix[1][0] = v[1].objective(j) - v[0].objective(j); //observe value delta(f)
			m_post_estimate[j] = m_pre_estimate[j] + m_K * (observe_Matrix - m_H * m_pre_estimate[j]);
		}
		/*for (size_t seq = 0; seq < v.size(); seq++) {
			std::vector<Matrix> temp_Matrix(v[0][0].objectiveSize());
			for (size_t j = 0; j < v[0][0].objectiveSize(); j++) {
				Matrix state;
				Matrix observe_Matrix(2, 1);
				observe_Matrix[0][0] = v[seq][2].objective(j);
				observe_Matrix[1][0] = v[seq][2].objective(j) - v[seq][1].objective(j);
				state = m_pre_estimate[seq][j] + m_K * (observe_Matrix - m_H * m_pre_estimate[seq][j]);
				temp_Matrix[j] = state;
			}
			m_post_estimate[seq] = temp_Matrix;
		}*/
		/*for (size_t seq = 0; seq < v.size(); seq++) {
			std::vector<Matrix> temp_Matrix;
			for (size_t j = 0; j < v[0][0].objectiveSize(); j++) {
				Matrix temp;
				Matrix observe_Matrix(2, 1);
				observe_Matrix[0][0] = v[seq][2].objective(j);
				observe_Matrix[1][0] = v[seq][2].objective(j) - v[seq][1].objective(j);
				temp = m_pre_estimate[seq][j] + m_K * (observe_Matrix - m_H * m_pre_estimate[seq][j]);
				temp_Matrix.emplace_back(temp);
			}
			m_post_estimate[seq] = temp_Matrix;
		}*/
	}
	template <typename Solution>
	const std::vector<Real>& kalmanFilter<Solution>::predict(const std::vector<Solution>& v) {
		if (m_type == true) {
			updatePara(v);	//update paremeter(m_K,m_P) and posterior estimate state (m_post_estimate)	
		}
		if (m_type == false) {
			initialize(v);
			m_type = true;
		}
		for (size_t j = 0; j < v[0].objectiveSize(); ++j) {
			m_pre_estimate[j] = m_A * m_post_estimate[j];	//update prediction estimate (m_pre_estimate)
			m_predict_result[j] = m_pre_estimate[j][0][0];	//record prediction value
		}
		return m_predict_result;
	}


	template <typename Solution = ofec::IndDE>
	class leastSquare {
	protected:
		std::vector<Real> m_predict_result;//store predict result
	public:
		//LeastSquare() {  }
		const std::vector<Real>& predict_L_theta(const std::vector<Solution>& v);	//predict according to L and theta
		const std::vector<Real>& predict_liner(const std::vector<Solution>& v);//predict according to every dimension objective, use y=ax+b
		const std::vector<Real>& predict_e(const std::vector<Solution>& v);//predict according to every dimension objective, use y=b*e^(ax)
		const std::vector<Real>& predict_liner_e(const std::vector<Solution>& v);//predict according to every dimension objective, use y=ax+b and y=b*e^(ax)
		Real least_square_estimation_liner(const std::vector<Real>& data);	//y=ax+b
		Real least_square_estimation_e(const std::vector<Real>& data);		//y=b*e^(ax), lny=ax+lnb
	};

	template <typename Solution>
	const std::vector<Real>& leastSquare<Solution>::predict_L_theta(const std::vector<Solution>& v) {
		if (v.size() < 3)
			std::cout << "prediction data too little" << std::endl;
		int M = v[0].objectiveSize();
		m_predict_result.resize(M);
		std::vector<std::vector<Real>> predict_result_2D(v.size() - 1, std::vector<Real>(2));	//(f1,f2),(f1',f3')...(f1',fM')
		for (int j = 1; j < M; ++j) {
			std::vector<Real> length(v.size() - 1);	//record norm between Solutions
			std::vector<Real> angle(v.size() - 1);	//record angle between norm and f1
			//std::vector<Real> predict_result(2);
			for (int i = 0; i < v.size() - 1; ++i) {
				length[i] = pow(fabs(v[i + 1].objective(0) - v[i].objective(0)), 2) + pow(fabs(v[i + 1].objective(j) - v[i].objective(j)), 2);	// (f1,fj)
				length[i] = sqrt(length[i]);
				if (fabs(length[i]) < 1e-10)
					angle[i] = 0;
				else {
					angle[i] = fabs(v[i + 1].objective(0) - v[i].objective(0)) / length[i];
					angle[i] = acos(angle[i]);
				}
			}
			//use least square to predict the position (f1,fj) by (length, index) and (angle, index)
			Real l = least_square_estimation_liner(length);
			if (l < 0)
				l = 0;
			Real theta = least_square_estimation_liner(angle);
			if (theta > acos(0))
				theta = acos(0);
			/*predict_result[0] = v.back().objective(0) + l * cos(theta);
			predict_result[1] = v.back().objective(j) + l * sin(theta);*/
			predict_result_2D[j - 1][0] = v.back().objective(0) + l * cos(theta);
			predict_result_2D[j - 1][1] = v.back().objective(j) + l * sin(theta);
			//calculate prediction result
			if (j == 1) {
				m_predict_result[j - 1] = predict_result_2D[j - 1][0];
				m_predict_result[j] = predict_result_2D[j - 1][1];
			}
			else
				m_predict_result[j] = predict_result_2D[j - 1][0] * predict_result_2D[0][1] / predict_result_2D[0][0]; //f1/f2=f1'/fj', fj=f1'*f2/f1
		}

		return m_predict_result;
	}

	template <typename Solution>
	const std::vector<Real>& leastSquare<Solution>::predict_liner(const std::vector<Solution>& v) {
		if (v.size() < 3)
			std::cout << "prediction data too little" << std::endl;
		int M = v[0].objectiveSize();
		m_predict_result.resize(M);
		for (int j = 0; j < M; ++j) {
			std::vector<Real> one_dim_obj(v.size());
			for (int i = 0; i < v.size(); ++i) {
				one_dim_obj[i] = v[i].objective(j);
			}
			m_predict_result[j] = least_square_estimation_liner(one_dim_obj);
		}
		return m_predict_result;
	}

	template <typename Solution>
	const std::vector<Real>& leastSquare<Solution>::predict_e(const std::vector<Solution>& v) {
		if (v.size() < 3)
			std::cout << "prediction data too little" << std::endl;
		int M = v[0].objectiveSize();
		m_predict_result.resize(M);
		for (int j = 0; j < M; ++j) {
			std::vector<Real> one_dim_obj(v.size());
			for (int i = 0; i < v.size(); ++i) {
				one_dim_obj[i] = v[i].objective(j);
			}
			m_predict_result[j] = least_square_estimation_e(one_dim_obj);
		}
		return m_predict_result;
	}

	template <typename Solution>
	const std::vector<Real>& leastSquare<Solution>::predict_liner_e(const std::vector<Solution>& v) {
		if (v.size() < 3)
			std::cout << "prediction data too little" << std::endl;
		int M = v[0].objectiveSize();
		m_predict_result.resize(M);
		for (int j = 0; j < M; ++j) {
			bool flag = true;
			std::vector<Real> one_dim_obj(v.size());
			for (int i = 0; i < v.size(); ++i) {
				one_dim_obj[i] = v[i].objective(j);
			}
			if (one_dim_obj[0] < one_dim_obj[1])	//it is increasing, then use liner, i.e., y=ax+b
				flag = false;
			if (flag)
				m_predict_result[j] = least_square_estimation_e(one_dim_obj);
			else
				m_predict_result[j] = least_square_estimation_liner(one_dim_obj);
		}
		return m_predict_result;
	}

	template <typename Solution>
	Real leastSquare<Solution>::least_square_estimation_liner(const std::vector<Real>& data) {
		Real sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
		for (size_t i = 0; i < data.size(); ++i) {
			sum1 += ((i + 1) * data[i]);
			sum2 += (i + 1);
			sum3 += data[i];
			sum4 = sum4 + (i + 1) * (i + 1);
		}
		Real a = (sum1 - sum2 * sum3 / data.size()) / (sum4 - sum2 * sum2 / data.size());
		Real b = sum3 / data.size() - a * sum2 / data.size();

		return a * (data.size() + 1) + b;
	}

	template <typename Solution>
	Real leastSquare<Solution>::least_square_estimation_e(const std::vector<Real>& data) {
		std::vector<Real> data_(data.size());
		for (size_t i = 0; i < data_.size(); ++i) {
			data_[i] = log(data[i]); //ln(data[i])
		}
		Real sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0;
		for (size_t i = 0; i < data_.size(); ++i) {
			sum1 += ((i + 1) * data_[i]);
			sum2 += (i + 1);
			sum3 += data_[i];
			sum4 = sum4 + (i + 1) * (i + 1);
		}
		Real a = (sum1 - sum2 * sum3 / data_.size()) / (sum4 - sum2 * sum2 / data_.size());
		Real b = sum3 / data_.size() - a * sum2 / data_.size();
		Real result = a * (data_.size() + 1) + b;
		result = exp(result);
		return result;
	}
}

#endif
