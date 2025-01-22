#ifndef GENERAL_MULTITHREAD_H
#define GENERAL_MULTITHREAD_H


#include <string>
#include <ostream>
#include <functional>
#include <mutex>
#include <iostream>
#include "../../core/exception.h"

#include <list>
#include <thread>
#include <vector>

namespace ofec {

	class GeneralMultiThreadInfo {
	public:
		std::string m_error_name;
		virtual ~GeneralMultiThreadInfo() {};
	};



	class GeneralMultiThread {
		std::mutex m_error_out_mutex;
		int m_numThread = 1;
		
		
	protected:
	
		void try_catch_fun(
			std::unique_ptr<GeneralMultiThreadInfo>& curInfo, 
			std::unique_ptr<GeneralMultiThreadInfo>& threadInfo, std::ostream& out) {
			//ms_fun(curInfo, ms_global_info);
			try {
				ms_fun(curInfo, threadInfo, ms_global_info);

			}
			catch (const std::exception& e) {
				std::string info = std::string(e.what()) + " at " +
					curInfo->m_error_name;
				std::lock_guard<std::mutex> guard(m_error_out_mutex);
				if (ms_cout_flag) std::cout << info << std::endl;
				out << info << std::endl;
			}
			//catch (const ofec::Exception& e) {
			//	std::string info = std::string(e.what()) + " at " +
			//		curInfo->m_error_name;
			//	std::lock_guard<std::mutex> guard(m_error_out_mutex);
			//	if (ms_cout_flag) std::cout << info << std::endl;
			//	out << info << std::endl;
			//}

			{
				std::lock_guard<std::mutex> guard(m_error_out_mutex);
				if (ms_cout_flag) std::cout << "finish task" << curInfo->m_error_name << std::endl;
			}
		}
		void try_catch_fun(
			std::unique_ptr<GeneralMultiThreadInfo>& curInfo,
			std::unique_ptr<GeneralMultiThreadInfo>& threadInfo) {
			try {
				ms_fun(curInfo, threadInfo, ms_global_info);

			}
			catch (const std::exception& e) {
				std::string info = std::string(e.what()) + " at " +
					curInfo->m_error_name;
				std::lock_guard<std::mutex> guard(m_error_out_mutex);
				if (ms_cout_flag) std::cout << info << std::endl;
			}
			//catch (const ofec::Exception& e) {
			//	std::string info = std::string(e.what()) + " at " +
			//		curInfo->m_error_name;
			//	std::lock_guard<std::mutex> guard(m_error_out_mutex);
			//	if (ms_cout_flag) std::cout << info << std::endl;
			//}

			{
				std::lock_guard<std::mutex> guard(m_error_out_mutex);
				if (ms_cout_flag) std::cout << "finish task" << curInfo->m_error_name << std::endl;
			}
		}

		void runTask(std::unique_ptr<GeneralMultiThreadInfo>& threadInfo,std::ostream& out) {
			bool flagInfo(false);
			std::unique_ptr<GeneralMultiThreadInfo> curInfo;
			while (true) {
				flagInfo = false;
				{
					std::lock_guard<std::mutex> guard(ms_info_mtx);
					if (!ms_info_buffer.empty()) {

						swap(curInfo, ms_info_buffer.front());
						if (ms_cout_flag)  std::cout << "running\t" << curInfo->m_error_name << std::endl;
						ms_info_buffer.pop_front();
						flagInfo = true;
					}
				}
				if (!flagInfo)break;
				try_catch_fun(curInfo, threadInfo, out);
			}
		}
		void runTaskNoOut(std::unique_ptr<GeneralMultiThreadInfo>& threadInfo) {
			bool flagInfo(false);
			std::unique_ptr<GeneralMultiThreadInfo> curInfo;
			while (true) {
				flagInfo = false;
				{
					std::lock_guard<std::mutex> guard(ms_info_mtx);
					if (!ms_info_buffer.empty()) {

						swap(curInfo, ms_info_buffer.front());
						if (ms_cout_flag)  std::cout << "running\t" << curInfo->m_error_name << std::endl;
						ms_info_buffer.pop_front();
						flagInfo = true;
					}
				}
				if (!flagInfo)break;
				try_catch_fun(curInfo, threadInfo);
			}
		}

		void runTaskNoException(std::unique_ptr<GeneralMultiThreadInfo>& threadInfo) {
			bool flagInfo(false);
			std::unique_ptr<GeneralMultiThreadInfo> curInfo;
			while (true) {
				flagInfo = false;
				{
					std::lock_guard<std::mutex> guard(ms_info_mtx);
					if (!ms_info_buffer.empty()) {

						swap(curInfo, ms_info_buffer.front());
						ms_info_buffer.pop_front();
						if (ms_cout_flag)  std::cout << "running\t" << curInfo->m_error_name << std::endl;
						flagInfo = true;
					}
				}
				if (!flagInfo)break;
				ms_fun(curInfo, threadInfo, ms_global_info);
			}
		}

	public:
		//GeneralMultiThread() = delete;
		bool ms_cout_flag = false;
		std::mutex ms_info_mtx;
		std::unique_ptr<GeneralMultiThreadInfo> ms_global_info;
		std::vector<std::unique_ptr<GeneralMultiThreadInfo>> ms_theadInfo;
		std::list<std::unique_ptr<GeneralMultiThreadInfo>> ms_info_buffer;

		std::function<void(
			std::unique_ptr<GeneralMultiThreadInfo>& curInfo,
			std::unique_ptr<GeneralMultiThreadInfo>& threadInfo,
			std::unique_ptr< GeneralMultiThreadInfo>& globalInfo)> ms_fun;




		GeneralMultiThread() {
			setNumberThread(std::thread::hardware_concurrency());
		}


		void setNumberThread(int numThread) {
			m_numThread = numThread;
			ms_theadInfo.resize(m_numThread);
		}
		void setDefaultNumberThread() {
			setNumberThread(std::thread::hardware_concurrency());
		}
		int numberThread()const { return m_numThread; }



		void run(std::ostream& out) {

			std::vector<std::thread> thrds;
			int num_task = m_numThread;
			//	num_task = 1;
			for (size_t i = 0; i < num_task; ++i) {
				thrds.push_back(std::thread(
					&GeneralMultiThread::runTask, this, std::ref(ms_theadInfo[i]), std::ref(out)));
			}
			for (auto& thrd : thrds)
				thrd.join();
		}


		void run() {

			std::vector<std::thread> thrds;
			int num_task = m_numThread;
			//	num_task = 1;
			for (size_t i = 0; i < num_task; ++i) {
				thrds.push_back(std::thread(
					&GeneralMultiThread::runTaskNoOut, this ,std::ref(ms_theadInfo[i])));
			}
			for (auto& thrd : thrds)
				thrd.join();
		}


		void runNoException() {

			std::vector<std::thread> thrds;
			int num_task = m_numThread;
			//	num_task = 1;
			for (size_t i = 0; i < num_task; ++i) {
				thrds.push_back(std::thread(
					&GeneralMultiThread::runTaskNoException, this, std::ref(ms_theadInfo[i])));
			}
			for (auto& thrd : thrds)
				thrd.join();
		}
		static void assignThreads(int num_samples, int num_task, std::vector<int>& ids);
	};
}



#endif