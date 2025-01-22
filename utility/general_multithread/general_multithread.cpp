#include "general_multithread.h"
//
//std::mutex utility::GeneralMultiThread::m_error_out_mutex;
//bool utility::GeneralMultiThread::ms_cout_flag = false;
////std::unique_ptr<utility::GeneralMultiThreadInfo> utility::GeneralMultiThread::ms_global_info = nullptr;
//std::list<std::unique_ptr<utility::GeneralMultiThreadInfo>> utility::GeneralMultiThread::ms_info_buffer;
//std::mutex utility::GeneralMultiThread::ms_info_mtx;
//
//std::function<void(std::unique_ptr<utility::GeneralMultiThreadInfo>& curInfo)> utility::GeneralMultiThread::ms_fun = nullptr;
//


void ofec::GeneralMultiThread::assignThreads(int num_samples, int num_task, std::vector<int>& ids) {
	ids.clear();
	int rest = num_samples % num_task;
	int id1 = 0, id2 = 0;
	std::pair<int, int> from_to;
	ids.push_back(id2);
	for (size_t i = 0; i < num_task; ++i) {
		id1 = id2;
		id2 = id1 + num_samples / num_task + (rest-- > 0 ? 1 : 0);
		ids.push_back(id2);
	}
}