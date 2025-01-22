#ifndef OFEC_LOCAL_OPTIMA_NETWORK_H
#define OFEC_LOCAL_OPTIMA_NETWORK_H

#include "../../../../../core/instance.h"
#include "../../../../../core/problem/problem.h"
#include "../../../combination/LKH_origin/INCLUDE/LKH.h"
#include "../lon.h"





namespace ofec {

	class LocalOptimaNetworkTSP{


	//public:
	//	struct RecordInfo {
	//		int runId = 0;
	//		long long beforeFit = 0;
	//		int beforeId = 0;

	//		long long afterFit = 0;
	//		int afterId = 0;
	//	};
	protected:
		std::vector<lon::TraceEdge> m_historyEdge;
		std::vector<std::vector<int>> m_sols;
		std::vector<double> m_nodeFit;


		const int Runs = 1000;
		const int chainedLKLoop = 1e4;



		struct ThreadInfo {
			std::vector<lon::TraceEdge> m_historyEdge;
			std::vector<std::vector<int>> m_sols;
			std::vector<double> m_nodeFit;

			std::shared_ptr<ofec::Random> m_rnd;
			
			void initialize(Problem* pro,Random* rnd);
			void clear();
		
		};

		struct GlobalInfo {
			unsigned m_maxRun = 0;
			std::mutex m_mtx;
		};

		
		
	protected:
		void clear() {
			m_historyEdge.clear();
			m_sols.clear();
			m_nodeFit.clear();
		}
		//void addInputParameters() {}
		void samplingLKH(ThreadInfo& cur, GlobalInfo& globalInfo,  Problem* pro);

		void insertDatas(std::vector<ThreadInfo>& totalInfo, Problem* pro, Random* rnd);
		
	public:

		//void setPath(const std::string& filedir, const std::string& filename) {
		//	m_filedir = filedir;
		//	m_filename = filename;
		//	
		//}
		
		void sampleSingleThread(Problem* pro, Random* rnd);
		void sampleMultiThread(Problem* pro, Random* rnd);

		void trasferToLon(lon::LonInfo& lon);
		
		const std::vector<std::vector<int>>& sols()const {
			return m_sols;
		}
		const std::vector<double>& nodeFit()const {
			return m_nodeFit;
		}

		
	};
}


#endif //  OFEC_LOCAL_OPTIMA_NETWORK_H
