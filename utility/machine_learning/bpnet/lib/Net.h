/**
 * @author  Gavin
 * @date    2022/2/10
 * @Email   gavinsun0921@foxmail.com
 */

#ifndef OFEC_NET_H
#define OFEC_NET_H

#include "../../../random/newran.h"

namespace ofec {
    class Sample {
    private:
        std::vector<Real> m_feature, m_label;//feature is input, label is output
    public:
        Sample();

        Sample(const std::vector<Real>& feature, const std::vector<Real>& label);

        void display();
        std::vector<Real> getSampleFeatures() const { return m_feature; }
        std::vector<Real> getSampleLabels() const { return m_label; }
        const Real getSampleFeature(size_t inx) const { return m_feature[inx]; }
        const Real getSampleLabel(size_t inx) const { return m_label[inx]; }

        void setSampleFeatures(std::vector<Real> v) {  m_feature=v; }
        void setSampleLabels(std::vector<Real> v) { m_label=v; }
        void setSampleFeature(size_t inx,Real v) {  m_feature[inx]=v; }
        void setSampleLabel(size_t inx, Real v) { m_label[inx]=v; }
    
    };

    class Node {
    private:
        Real m_value, m_bias, m_bias_delta;
        std::vector<Real> m_weight, m_weight_delta;
    public:
        explicit Node(size_t nextLayerSize);
        std::vector<Real> getNodeWeight() { return m_weight; }
        std::vector<Real> getNodeWeightDelta() { return m_weight_delta; }
        Real getNodeWeight(size_t inx) const{ return m_weight[inx]; }
        Real getNodeWeightDelta(size_t inx) const { return m_weight_delta[inx]; }
        void setNodeWeight(size_t inx,Real v) { m_weight[inx] = v; }
        void setNodeWeightDelta(size_t inx, Real v) { m_weight_delta[inx] = v; }

        Real getNodeBias() const { return m_bias; }
        Real getNodeValue() const { return m_value; }
        Real getNodeBiasDelta() const { return m_bias_delta; }

        void setNodeBias(Real v) {m_bias=v; }
        void setNodeValue(Real v) { m_value=v; }
        void setNodeBiasDelta(Real v) { m_bias_delta=v; }
    };

    class Net {
    private:
        size_t m_num_input_node;
        size_t m_num_hiden_node;
        size_t m_num_output_node;

        Real m_learning_rate;
        Real m_threshold;
        size_t m_max_epoch;

        std::vector<Real> m_his_loss;
        /*Real m_learning_rate = 0.8;
        Real m_threshold = 1e-4;
        size_t m_max_epoch = 1e5;*/

        std::vector<std::unique_ptr<Node>> inputLayer;
        std::vector<std::unique_ptr<Node>> hiddenLayer;
        std::vector<std::unique_ptr<Node>> outputLayer;

        /**
         * Clear all gradient accumulation
         *
         * Set 'weight_delta'(the weight correction value) and
         * 'bias_delta'(the bias correction value) to 0 of nodes
         */
        void grad_zero();

        /**
         * Forward propagation
         */
        void forward();

        /**
         * Calculate the value of Loss Function
         * @param label the label of sample (std::vector / numeric)
         * @return loss
         */
        Real calculateLoss(const std::vector<Real>& label);

        /**
         * update learning rate according to delta of loss
         * @param label the label of sample (Real / numeric)
         * @return loss
         */
        void updateLearningRate(Real loss);

        /**
         * Back propagation
         * @param label label of sample (std::vector / numeric)
         */
        void backward(const std::vector<Real>& label);

        /**
         * Revise 'weight' and 'bias according to
         * 'weight_delta'(the weight correction value) and
         * 'bias_weight'(the bias correction value)
         * @param batch_size
         */
        void revise(size_t batch_size);

    public:

        Net(size_t input_layer, size_t hiden_layer, size_t output_layer,Real lr,Real threshold,size_t epoch, Random* rnd);

        /**
         * Training network with training data
         * @param trainDataSet The sample set
         * @return Training convergence
         */
        bool train(const std::vector<Sample>& trainDataSet);

        /**
         * Using network to predict sample
         * @param feature The feature of sample (std::vector)
         * @return Sample with 'feature' and 'label'(predicted)
         */
        Sample predict(const std::vector<Real>& feature);

        /**
         * Using network to predict the sample set
         * @param predictDataSet The sample set
         * @return The sample set, in which each sample has 'feature' and 'label'(predicted)
         */
        std::vector<Sample> predict(const std::vector<Sample>& predictDataSet);

        
        Real getLearningRate() { return m_learning_rate; }
        Real getThreshold() { return m_threshold; }
        size_t getMaxEpoch() { return m_max_epoch; }
        std::vector<Real> getHisLoss() { return m_his_loss; }
        void setLearningRate(Real v) { m_learning_rate=v; }
        void setThreshold(Real v) {  m_threshold=v; }
        void setMaxEpoch(size_t v) {  m_max_epoch=v; }
    };
}

#endif // !OFEC_NET_H
