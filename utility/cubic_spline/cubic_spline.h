/**
 * @author Jose Davis Nidhin
 * modified by DiaoYiya
 */

#pragma once

#include <math.h>
#include <iostream>
#include<vector>



namespace ofec {
	template <typename T>
	class CubicSpline
	{
	public:
		CubicSpline() {};
		~CubicSpline() {};


		void Initialize(const std::vector<T> &srcX,const std::vector<T>& srcY);
		T Interpolate(T x);

		void set_interval(T interval, T from) {
			m_interval = interval;
			m_from = from;
		}

		T Interpolate_static_inte(T x);

	
	private:
		std::vector<T> m_B,  m_C,  m_D,  m_X,  m_Y;
		int m_size;
		T m_interval = 0.;
		T m_from = 0;

		int GetIndex(T x);
		T Interpolate(T x, int index);
	};

	template<typename T>
	inline T CubicSpline<T>::Interpolate_static_inte(T x)
	{
		int index =  (x-m_from)/ m_interval;
		if (index >= m_size) {
			index = m_size - 1;
		}
		////int index = static_cast<int>((x - m_from + 0.5* m_interval) / m_interval);
		//if (index < 0) {
		//	index = 0;
		//}
		//else if (index >= m_size) {
		//	index = m_size - 1;
		//}
		return Interpolate(x, index);
	}


	template<typename T>
	inline void CubicSpline<T>::Initialize(const std::vector<T>& srcX, const std::vector<T>& srcY)
	{

		m_X = srcX;
		m_Y = srcY;
		m_size = srcX.size();
		if (m_size == 1) {
			m_B.resize(m_size);
			m_C.resize(m_size);
			m_D.resize(m_size);
			m_B[0] = m_C[0] = m_D[0] = 0;
		}
		else if (m_size == 2) {
			m_B.resize(m_size);
			m_C.resize(m_size);
			m_D.resize(m_size);
			for (int i(0); i < 2; ++i) {
				m_C[i] = m_D[i] = 0;
				m_B[i] = (m_Y[1] - m_Y[0]) / (m_X[1] - m_X[0]);
			}

		}
		else {



			m_B.resize(m_size );
			m_B[m_size - 1] = 0;
			m_C.resize(m_size);
			m_D.resize(m_size);
			m_D[m_size - 1] = 0;

			std::vector<float> H(m_size - 1);
			std::vector<float> alpha(m_size - 2);

			std::vector<float> L(m_size - 2);
			std::vector<float> U(m_size - 1);
			std::vector<float> Z(m_size - 1);

			for (int i = 0; i < m_size - 1; i++) {
				H[i] = srcX[i + 1] - srcX[i];
			}

			for (int i = 1; i < m_size - 1; i++) {
				alpha[i - 1] = (3.0f / H[i] * (srcY[i + 1] - srcY[i])) - (3.0f / H[i - 1] * (srcY[i] - srcY[i - 1]));
			}

			U[0] = Z[0] = 0.0f;

			for (int i = 1; i < m_size - 1; i++) {
				L[i - 1] = (2.0f * (srcX[i + 1] - srcX[i - 1])) - (H[i - 1] * U[i - 1]);
				U[i] = H[i] / L[i - 1];
				Z[i] = (alpha[i - 1] - (H[i - 1] * Z[i - 1])) / L[i - 1];
			}

			m_C[m_size - 1] = 0.0f;

			for (int i = m_size - 2; i >= 0; i--) {
				m_C[i] = Z[i] - (U[i] * m_C[i + 1]);
				m_B[i] = ((srcY[i + 1] - srcY[i]) / H[i]) - (H[i] * (m_C[i + 1] + 2 * m_C[i]) / 3);
				m_D[i] = (m_C[i + 1] - m_C[i]) / (3 * H[i]);
			}

		}

	}
	template<typename T>
	inline T CubicSpline<T>::Interpolate(T x)
	{
		int index = GetIndex(x);
		return Interpolate(x, index);
	}

	template<typename T>
	inline int CubicSpline<T>::GetIndex(T x)
	{
		int last = m_size - 1;
		int first = 0;
		int mid;

		while (last - first > 1) {
			mid = (last + first) / 2;

			if (x > m_X[mid])
				first = mid;
			else
				last = mid;
		}

		return first;
	}
	template<typename T>
	inline T CubicSpline<T>::Interpolate(T x, int index)
	{
		T y(0);
		y = m_Y[index] +
			m_B[index] * (x - m_X[index]) +
			m_C[index] * pow((x - m_X[index]), 2) +
			m_D[index] * pow((x - m_X[index]), 3);

		return y;
	}
}


