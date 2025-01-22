#include "data_set.h"

using namespace torch;

namespace ofec {

	extern Device g_device;

	Dataset::Dataset(int num_data, const torch::Tensor &x, const torch::Tensor &y) :
		m_num_data(num_data),
		m_inputs(x),
		m_outputs(y) {}

	data::Example<> Dataset::get(size_t index) {
		return { m_inputs[index], m_outputs[index] };
	}

	optional<size_t> Dataset::size() const {
		return m_num_data;
	}
}
