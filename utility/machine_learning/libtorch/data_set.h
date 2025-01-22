#ifndef OFEC_DATA_SET_H
#define OFEC_DATA_SET_H


#ifdef  OFEC_DEMO_LITTORCH
#include <torch/torch.h>

namespace ofec {

	extern torch::Device g_device;

	class Dataset : public torch::data::Dataset<Dataset> {
	private:
		torch::Tensor m_inputs, m_outputs;
		int m_num_data;

	public:
		Dataset(int num_data, const torch::Tensor &x, const torch::Tensor &y);
		torch::data::Example<> get(size_t index) override;
		torch::optional<size_t> size() const override;
	};

}
#endif // !OFEC_DATA_SET_H

#endif
