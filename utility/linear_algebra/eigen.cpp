#include "eigen.h"
#include "../../core/exception.h"
#include <unordered_set>
#include <numeric>
#include <algorithm>

using namespace Eigen;

namespace ofec {
	/*-----------------------------------------------------------------------------------------------
	* https://www.mathworks.com/help/matlab/ref/prctile.html
	* -----------------------------------------------------------------------------------------------
	* For an n-element vector A, prctile returns percentiles by using a sorting-based algorithm:
	*	1. The sorted elements in A are taken as the 100(0.5/n)th, 100(1.5/n)th, ..., 100([n ¨C 0.5]/n)th percentiles. For example:
	*		- For a data vector of five elements such as {6, 3, 2, 10, 1}, the sorted elements {1, 2, 3, 6, 10} respectively correspond to the 10th, 30th, 50th, 70th, and 90th percentiles.
	*		- For a data vector of six elements such as {6, 3, 2, 10, 8, 1}, the sorted elements {1, 2, 3, 6, 8, 10} respectively correspond to the (50/6)th, (150/6)th, (250/6)th, (350/6)th, (450/6)th, and (550/6)th percentiles.
	*	2. prctile uses linear interpolation to compute percentiles for percentages between 100(0.5/n) and 100([n ¨C 0.5]/n).
	*	3. prctile assigns the minimum or maximum values of the elements in A to the percentiles corresponding to the percentages outside that range.
	* prctile treats NaNs as missing values and removes them.
	*------------------------------------------------------------------------------------------------*/
	double prctile(const Ref<const VectorXd> &vec, size_t p) {
		size_t n = vec.size();
		if (n == 0)
			throw Exception("Vector must not be empty@ Vector::prctile");
		else if (n == 1)
			return vec(0);
		else {
			VectorXd s = vec;
			std::sort(s.begin(), s.end());
			std::vector<double> q(n);
			for (size_t i = 0; i < n; ++i)
				q[i] = (i + 0.5) / n;
			double q_p = p / 100.0;
			if (q_p < q.front())
				return s(0);
			else if (q_p > q.back())
				return s(indexing::last);
			else {
				size_t l = floor(q_p * n), r = l + 1;
				return s[l] + (q_p - q[l]) / (q[r] - q[l]) * (s[r] - s[l]);
			}
		}
	}

	double median(const Ref<const VectorXd> &vec) {
		return prctile(vec, 50);
	}

	MatrixXd randMatXd(int rows, int cols, Uniform *unif) {
		MatrixXd new_mat(rows, cols);
		for (size_t i = 0; i < rows; ++i) {
			for (size_t j = 0; j < cols; ++j) {
				new_mat(i, j) = unif->next();
			}
		}
		return new_mat;
	}

	VectorXd randVecXd(int size, Uniform *unif) {
		VectorXd new_vec(size);
		for (size_t i = 0; i < size; ++i)
			new_vec(i) = unif->next();
		return new_vec;
	}

	VectorXd randnVecXd(int size, Normal *norm) {
		VectorXd new_vec(size);
		for (size_t i = 0; i < size; ++i)
			new_vec(i) = norm->next();
		return new_vec;
	}
	
	VectorXi join(const Ref<const VectorXi> &vec, int val) {
		VectorXi new_vec(vec.size() + 1);
		new_vec << vec, val;
		return new_vec;
	}

	VectorXi join(const Ref<const VectorXi> &vec1, const Ref<const VectorXi> &vec2) {
		VectorXi new_vec(vec1.size() + vec2.size());
		new_vec << vec1, vec2;
		return new_vec;
	}

	VectorXd join(const Ref<const VectorXd> &vec, double val) {
		VectorXd new_vec(vec.size() + 1);
		new_vec << vec, val;
		return new_vec;
	}

	VectorXd join(const Ref<const VectorXd> &vec1, const Ref<const VectorXd> &vec2) {
		VectorXd new_vec(vec1.size() + vec2.size());
		new_vec << vec1, vec2;
		return new_vec;
	}

	VectorXi setdiff(const Ref<const VectorXi> &vec1, const Ref<const VectorXi> &vec2) {
		std::unordered_set<int> set;
		for (int val : vec1) set.insert(val);
		for (int val : vec2) set.erase(val);
		VectorXi new_vec(set.size());
		size_t i = 0;
		for (int val : set) new_vec(i++) = val;
		return new_vec;
	}

	VectorXi sort(const Ref<const VectorXd> &vec, bool ascend) {
		VectorXi seq(vec.size());
		std::iota(seq.begin(), seq.end(), 0);
		if (ascend)
			std::sort(seq.begin(), seq.end(), [&vec](int i, int j) {return vec(i) < vec(j); });
		else
			std::sort(seq.begin(), seq.end(), [&vec](int i, int j) {return vec(i) > vec(j); });
		return seq;
	}

	VectorXi find(const Ref<const ArrayXb> &vec) {
		std::vector<int> v;
		for (int i = 0; i < vec.size(); ++i)
			if (vec(i))	v.push_back(i);
		return Map<VectorXi>(v.data(), v.size());
	}

	Eigen::VectorXd normcdf(const Ref<const VectorXd> &vec) {
		VectorXd result(vec.size());
		for (size_t i = 0; i < vec.size(); ++i)
			result[i] = 0.5 * erfc(-vec(i) / sqrt(2));
		return result;
	}

	MatrixXd joinRow(const Ref<const MatrixXd> &mat1, const Ref<const MatrixXd> &mat2) {
		if (mat1.cols() != mat2.cols())
			throw Exception("The matrices to be joined by rows must have the same number of columns@joinRow");
		MatrixXd new_mat(mat1.rows() + mat2.rows(), mat1.cols());
		new_mat << mat1, mat2;
		return new_mat;
	}

	MatrixXd joinCol(const Ref<const MatrixXd> &mat1, const Ref<const MatrixXd> &mat2) {
		if (mat1.rows() != mat2.rows())
			throw Exception("The matrices to be joined by columns must have the same number of rows@joinRow");
		MatrixXd new_mat(mat1.rows(), mat1.cols() + mat2.cols());
		new_mat << mat1, mat2;
		return new_mat;
	}
}