#ifndef MIN_HASH_LSH_H
#define MIN_HASH_LSH_H

#include<vector>
#include<map>
#include<array>

namespace ofec {


	template<typename TInfo>
	class LSH_info {
		int m_fre = 0;
		double m_obj = 0;
		long long  m_com_times = 0;
		TInfo m_info;

		void init(const TInfo & info) {
			m_fre = 0;
			m_obj = 0;
			m_com_times = 0;
			m_info = info;
		}
	};

    template<typename TInfo, size_t TNum=1>
	class MinHashLSH {

		std::vector<LSH_info<TInfo>> m_infos;
		std::vector<int> m_empty_idxs;
		std::vector<std::map<std::array<int,TNum>, LSH_info<TInfo>>> m_buckets;

	protected:
		int assignIds() {
			if (!m_empty_idxs.empty()) {
				int cur = m_empty_idxs.back();
				m_empty_idxs.pop_back();
				return cur;
			}
			else {
				m_infos.push_back(LSH_info<TInfo>());
				return m_infos.size() - 1;
			}
		}

		void deleteInfo() {
//			return 
		}
	public:
		//	void update(std::vector<std::array<int, TNum>>& num, TInfo& info);
		
		    
	};
	//template<typename TInfo, size_t TNum>
	//inline void MinHashLSH<TInfo, TNum>::update(std::vector<std::array<int, TNum>>& num, TInfo& info)
	//{
	//	
	//}
}

#endif