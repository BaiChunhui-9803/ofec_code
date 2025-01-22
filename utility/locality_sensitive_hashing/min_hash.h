#ifndef  MIN_HASH_H
#define  MIN_HASH_H

#include<vector>
#include<memory>
//#include"../hash-library/sha1.h"
#include"../hash-library/hash.h"
#include"../random/newran.h"

namespace ofec {
	class MinHash {
	public:

	protected:
		/*
   _mersenne_prime = np.uint64((1 << 61) - 1)
_max_hash = np.uint64((1 << 32) - 1)
_hash_range = (1 << 32)
		 */
		int m_num_perm = 128;
		//	std::vector<int> m_hash_values;
		std::vector<std::pair<unsigned long long, unsigned long long>> m_permutations;
		unsigned long long m_prime = 11003;
		unsigned long long m_max_hash = long(1 << 32) - 1;
		//	unsigned long long m_hash_range = long long(1 << 32);
		std::unique_ptr<Hash> m_hash_fun;

	public:
		unsigned long long numberPrime()const {
			return m_prime;
		}
		void setNumberPerm(int num_perm) {
			m_num_perm = num_perm;
		}

		void setNumberPrime(unsigned long long numPrime) {
			m_prime = numPrime;
		}
		void initialize(Random *rnd, HashType hash_type);
		std::vector<int> getMinHash(std::vector<std::string>& vals)const;
		std::vector<int> getMinHash(std::vector<unsigned long long>& vals)const;
	};
}

#endif // ! MIN_HASH_H
