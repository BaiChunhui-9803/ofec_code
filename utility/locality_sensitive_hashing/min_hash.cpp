#include "min_hash.h"
#include "../hash-library/hash_generator.h"
#include "../function/custom_function.h"

namespace ofec {

	void MinHash::initialize(Random *rnd, HashType hash_type) {
		HashGenerator(m_hash_fun, hash_type);
		//m_hash_values.resize(m_num_perm);
		m_permutations.resize(m_num_perm);
		for (int idx(0); idx < m_num_perm; ++idx) {
			m_permutations[idx].first = rnd->uniform.nextNonStd<unsigned long long>(1e3, m_prime);
			m_permutations[idx].second = rnd->uniform.nextNonStd<unsigned long long>(1e3, m_prime);
		}
		m_max_hash = static_cast<unsigned long long>(powl(2, 31)) - 1;

		//	m_max_hash = static_cast<unsigned long long>(1) << 31 - 1;

		//	int stop = -1;

	}
	std::vector<int> MinHash::getMinHash(std::vector<std::string>& vals)const
	{
		std::vector<int> min_hash(m_permutations.size(), m_max_hash);
		std::vector<long long> valInt(vals.size(), 0);
		for (int idx(0); idx < valInt.size(); ++idx) {
			valInt[idx] = UTILITY::HexadecimalStrToInt((*m_hash_fun)(vals[idx]));
		}

		long long phv(0);
		for (int idp(0); idp < m_permutations.size(); ++idp) {
			for (int idv(0); idv < vals.size(); ++idv) {
				phv = valInt[idv];
				phv = ((m_permutations[idp].first * phv + m_permutations[idp].second) % m_prime) & m_max_hash;
				min_hash[idp] = std::min<int>(phv, min_hash[idp]);
			}
		}
		return std::move(min_hash);
	}



	std::vector<int> MinHash::getMinHash(std::vector<unsigned  long long>& vals)const
	{
		std::vector<int> min_hash(m_permutations.size(), m_max_hash);
		const auto & valInt = vals;
		//for (int idx(0); idx < valInt.size(); ++idx) {
		//	valInt[idx] = UTILITY::HexadecimalStrToInt((*m_hash_fun)(std::to_string(vals[idx])));
		//}
		long long phv(0);
		for (int idp(0); idp < m_permutations.size(); ++idp) {
			for (int idv(0); idv < vals.size(); ++idv) {
				phv = valInt[idv];
				phv = ((m_permutations[idp].first * phv + m_permutations[idp].second) % m_prime) & m_max_hash;
				min_hash[idp] = std::min<int>(phv, min_hash[idp]);
			}
		}
		return std::move(min_hash);
	}
}