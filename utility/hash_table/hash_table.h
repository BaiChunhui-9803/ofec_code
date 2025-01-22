#ifndef  UTILITY_HASH_TABLE_H
#define  UTILITY_HASH_TABLE_H

#include<vector>
#include "../function/custom_function.h"
#include "../../core/random/newran.h"
#include <algorithm>
#include <iomanip>
#include <map>

namespace utility {



	struct HashValueSequenceMapper {
		std::vector<unsigned> m_hash_val;
		int m_from = 0;
		int m_to;
		void initialize(ofec::Random *rnd, int from, int to) {
			m_from = from;
			m_to = to;
			m_hash_val.resize(to-from+1);
			unsigned maxVal = std::numeric_limits<unsigned>::max();
			maxVal = sqrt(maxVal);
			for (auto& it : m_hash_val) {
				it = rnd->uniform.nextNonStd<unsigned>(0, maxVal);
			}

		}
		unsigned getHash(int idx) {
			return m_hash_val[idx - m_from];
		}
		unsigned calTSPHash(const std::vector<int>& sol)const {
			unsigned Hash(0);
			for (int idx(1); idx < sol.size(); ++idx) {
				Hash ^= m_hash_val[sol[idx - 1]-m_from] * m_hash_val[sol[idx]-m_from];
			}
			Hash ^= m_hash_val[sol.back()- m_from] * m_hash_val[sol.front()- m_from];
			return Hash;
		}

	};

	template<typename T>
	struct HashMap {

		//unsigned m_solId = 0;
		std::vector<std::pair<unsigned, std::vector<int>>> m_hash_ids;
		std::vector<T> m_totalSol;
		std::map<unsigned, unsigned> m_hashToIds;

		void initialize() {
			m_hash_ids.clear();
			m_totalSol.clear();
			m_hashToIds.clear();
		}
		// return flag_exist
		bool insertValue(const T& val, unsigned hashVal, int& id) {
			if (m_hashToIds.find(hashVal) == m_hashToIds.end()) {
				int mapId = m_hashToIds[hashVal] = m_hash_ids.size();
				id = m_totalSol.size();
				m_totalSol.push_back(val);
				std::vector<int> curIds;
				std::pair<unsigned, std::vector<int>> curInfo;
				curInfo.first = hashVal;
				curInfo.second = curIds;
				curIds.push_back(id);
				m_hash_ids.push_back(curInfo);
				return false;
			}
			else {
				int mapId = m_hashToIds[hashVal];
				id = -1;
				for (auto& posId : m_hash_ids[mapId].second) {
					if (m_totalSol[posId] == val) {
						id = posId;
						break;
					}
				}
				if (id == -1) {
					id = m_totalSol.size();
					m_totalSol.push_back(val);
					m_hash_ids[mapId].second.push_back(id);
					return false;

				}
				else {
					return true;
				}
			}
		}
	};

	//template<typename T= unsigned >
	struct HashRandomValue {
		std::vector<unsigned long long> m_hash_randNum;
		void initialize(int num, ofec::Random *rnd, 
			unsigned long long from ,unsigned long long to) {
			m_hash_randNum.resize(num);
			for (auto& it : m_hash_randNum) {
				it = rnd->uniform.nextNonStd<unsigned long long>(from, to);
			}
		}
		
	};


	struct HashNode {

		bool m_divided_flag = false;
		//int m_div_idx;
		//unsigned m_hash_divide_num = 997;
		std::vector<std::pair<unsigned,std::vector<int>>> m_hash_ids;
		std::vector<HashNode> m_hash_nodes;

//		bool find(unsigned hashValue, int newId, )

	//	void initialize()
	};


	template<typename T> 
	class HashTable {

		int m_divide_threadhold = 1e3;
		unsigned m_hash_div_num;
		std::vector<HashNode> m_hash_nodes;

		std::vector<int> m_prime_range;
		//std::vector<unsigned> m_hash_values;
		std::vector<T> m_hash_info;

		int m_maxPrimeNumber = 2e6;
		

	public:
		void initiazlize() {

			std::vector<int> primeNumber;
			UTILITY::getPrimeNumber(m_maxPrimeNumber, primeNumber);
		//	std::cout << primeNumber.back() << std::endl;
			m_hash_div_num = primeNumber.back();

			int minNum(1e3), maxNum(1e3+50);
			auto fromIter = std::lower_bound(primeNumber.begin(), primeNumber.end(), minNum);
			auto toIter = std::lower_bound(primeNumber.begin(), primeNumber.end(), maxNum);
			
			m_prime_range.clear();
			for (auto cur = fromIter; cur != toIter; ++cur) {
				m_prime_range.push_back(*cur);
			}
			m_hash_nodes.resize(m_hash_div_num);
	//		int stop = -1;
		}

		bool insertValue(const T& val, unsigned hashVal,int& id) {
			HashNode* curNode = &m_hash_nodes[hashVal % m_hash_div_num];
			int curDivIdx(-1);
			while (curNode->m_divided_flag && curDivIdx+1 < m_prime_range.size()) {
				++curDivIdx;
				curNode = &curNode->m_hash_nodes[hashVal % m_prime_range[curDivIdx]];
				//++curDivIdx;
			}
			std::vector<int>* nearIds = nullptr;
			for (auto& nearInfo : curNode->m_hash_ids) {
				if (nearInfo.first == hashVal) {
					nearIds = &nearInfo.second;
				}
			}
			id = -1;
			if (nearIds != nullptr) {
				for (const auto& nearId : *nearIds) {
					if (m_hash_info[nearId] == val) {
						id = nearId;
					}
				}
			}

			bool flag_exist = true;
			if (id == -1) {
				id = m_hash_info.size();
				m_hash_info.push_back(val);
//				m_hash_info[id] = val;
				if (nearIds == nullptr) {
					std::pair<unsigned, std::vector<int>> cur_info;
					cur_info.first = hashVal;
					cur_info.second.push_back(id);
					curNode->m_hash_ids.emplace_back(cur_info);
				}
				else {
					nearIds->push_back(id);
				}
				flag_exist = false;
			}


			if (!curNode->m_divided_flag&&curNode->m_hash_ids.size() > m_divide_threadhold && curDivIdx + 1 < m_prime_range.size()) {
				++curDivIdx;
				auto curDiv = m_prime_range[curDivIdx];
				curNode->m_divided_flag = true;
				curNode->m_hash_nodes.resize(curDiv);
				for (auto& curInfo : curNode->m_hash_ids) {
					int sonId = curInfo.first % curDiv;
					curNode->m_hash_nodes[sonId].m_hash_ids.emplace_back(curInfo);
				}
			}

			return flag_exist;
		}

		

	};
}

#endif // ! OFEC_HASH_TABLE_H
