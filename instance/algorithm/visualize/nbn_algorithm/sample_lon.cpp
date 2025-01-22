#include"sample_lon.h"

const int ofec_demo::SampleLON::mc_num_sample = 2e3;
const int ofec_demo::SampleLON::mc_num_exploration = 50;
const int ofec_demo::SampleLON::mc_num_exploration_sample=2e3;
const double ofec_demo::SampleLON::mc_num_exploration_radius = 0.3;
const int ofec_demo::SampleLON::NNnetworkQueue::mc_queue_size = 50;
void ofec_demo::SampleLON::toIdx(int idx, std::pair<int, int>& idxs)
{
	idxs.first = idx / mc_num_exploration_sample;
	idxs.second = idx % mc_num_exploration_sample;
}

void ofec_demo::SampleLON::toId(const std::pair<int, int>& idxs, int& id)
{
	id = idxs.first * mc_num_exploration_sample + idxs.second;
}
