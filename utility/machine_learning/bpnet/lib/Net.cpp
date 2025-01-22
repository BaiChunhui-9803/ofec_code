/**
 * @author  Gavin
 * @date    2022/2/10
 * @Email   gavinsun0921@foxmail.com
 */

#include "Net.h"
#include "Utils.h"

namespace ofec {
    Net::Net(size_t num_input_node, size_t num_hiden_node, size_t num_output_node, Real lr, Real threshold, size_t epoch, Random* rnd) :m_num_input_node(num_input_node), \
        m_num_hiden_node(num_hiden_node), m_num_output_node(num_output_node),m_learning_rate(lr),m_threshold(threshold),m_max_epoch(epoch){

        /*std::mt19937 rd;
        rd.seed(std::random_device()());

        std::uniform_real_distribution<Real> distribution(-1, 1);*/

        /**
         * Initialize input layer
         */
        for (size_t i = 0; i < num_input_node; ++i) {
            inputLayer.emplace_back(std::make_unique<Node>(num_hiden_node));

            for (size_t j = 0; j < num_hiden_node; ++j) {

                // Initialize 'weight'(the weight value)
                // from the i-th node in the input layer to the j-th node in the hidden layer
                Real temp = 2 * rnd->uniform.next() - 1;
                inputLayer.back()->setNodeWeight(j,temp);

                // Initialize 'weight_delta'(the weight correction value)
                // from the i-th node in the input layer to the j-th node in the hidden layer
                inputLayer.back()->setNodeWeightDelta(j,0.);
            }
        }

        /**
         * Initialize hidden layer
         */
        for (size_t j = 0; j < num_hiden_node; ++j) {
            hiddenLayer.emplace_back(std::make_unique<Node>(num_output_node));

            // Initialize 'bias'(the bias value)
            // of the j-th node in the hidden layer
            hiddenLayer.back()->setNodeBias(2 * rnd->uniform.next() - 1);

            // Initialize 'bias_delta'(the bias correction value)
            // of the j-th node in the hidden layer
            hiddenLayer.back()->setNodeBiasDelta(0.);
            for (size_t k = 0; k < num_output_node; ++k) {

                // Initialize 'weight'(the weight value)
                // from the j-th node in the hidden layer to the k-th node in the output layer
                hiddenLayer.back()->setNodeWeight(k, 2 * rnd->uniform.next() - 1);

                // Initialize 'weight_delta'(the weight correction value)
                // from the j-th node in the hidden layer to the k-th node in the output layer
                hiddenLayer.back()->setNodeWeightDelta(k, 0.);
            }
        }

        // Initialize output layer
        for (size_t k = 0; k < num_output_node; ++k) {
            outputLayer.emplace_back(std::make_unique<Node>(0));

            // Initialize 'bias'(the bias value)
            // of the k-th node in the output layer
            outputLayer.back()->setNodeBias(2 * rnd->uniform.next() - 1);

            // Initialize 'bias_delta'(the bias correction value)
            // of the k-th node in the output layer
            outputLayer.back()->setNodeBiasDelta(0.);
        }
    }

    void Net::grad_zero() {

        // Clear 'weight_delta'(the weight correction value)
        // of all nodes in the input layer
        for (auto& nodeOfInputLayer : inputLayer) {
            nodeOfInputLayer->getNodeWeightDelta().assign(nodeOfInputLayer->getNodeWeightDelta().size(), 0.);
        }

        // Clear 'weight_delta'(the weight correction value) and 'bias_delta'(the bias correction value)
        // of all nodes in the hidden layer
        for (auto& nodeOfHiddenLayer : hiddenLayer) {
            nodeOfHiddenLayer->setNodeBiasDelta(0.);
            nodeOfHiddenLayer->getNodeWeightDelta().assign(nodeOfHiddenLayer->getNodeWeightDelta().size(), 0.);
        }

        // Clear 'bias_delta'(the bias correction value)
        // of all nodes in the hidden layer
        for (auto& nodeOfOutputLayer : outputLayer) {
            nodeOfOutputLayer->setNodeBiasDelta(0.);
        }
    }

    void Net::forward() {

        /**
         * The input layer propagate forward to the hidden layer.
         * MathJax formula: h_j = \sigma( \sum_i x_i w_{ij} - \beta_j )
         */
        for (size_t j = 0; j < m_num_hiden_node; ++j) {
            Real sum = 0;
            for (size_t i = 0; i < m_num_input_node; ++i) {
                sum += inputLayer[i]->getNodeValue() * inputLayer[i]->getNodeWeight(j);
            }
            sum -= (hiddenLayer[j]->getNodeBias());

            hiddenLayer[j]->setNodeValue(ofec::Utils::sigmoid(sum));
        }

        /**
         * The hidden layer propagate forward to the output layer.
         * MathJax formula: \hat{y_k} = \sigma( \sum_j h_j v_{jk} - \lambda_k )
         */
        for (size_t k = 0; k < m_num_output_node; ++k) {
            Real sum = 0;
            for (size_t j = 0; j < m_num_hiden_node; ++j) {
                sum += (hiddenLayer[j]->getNodeValue() * hiddenLayer[j]->getNodeWeight(k));
            }
            sum -= (outputLayer[k]->getNodeBias());

            outputLayer[k]->setNodeValue(ofec::Utils::sigmoid(sum));
        }
    }

    Real Net::calculateLoss(const std::vector<Real>& label) {
        Real loss = 0.;

        /**
         * MathJax formula: Loss = \frac{1}{2}\sum_k ( y_k - \hat{y_k} )^2
         */
        for (size_t k = 0; k < m_num_output_node; ++k) {
            Real tmp = std::fabs(outputLayer[k]->getNodeValue() - label[k]);
            loss += (tmp * tmp / 2);
        }

        return loss;
    }

    void Net::updateLearningRate(Real loss) {
        auto& his_loss = getHisLoss();
        if (his_loss.size() > 0) {
            Real last_loss = his_loss.back();
            auto lr = getLearningRate();
            if (loss < last_loss) {
                if (lr < 1) {
                    setLearningRate(getLearningRate() * 1.01);
                }
                else {
                    setLearningRate(0.5);
                }
            }
            else if (loss >= last_loss) {
                if (lr > 0.1) {
                    setLearningRate(getLearningRate() * 0.99);
                }
                else {
                    setLearningRate(0.5);
                }
            }
        }
    }

    void Net::backward(const std::vector<Real>& label) {

        /**
         * Calculate 'bias_delta'(the bias correction value)
         * of the k-th node in the output layer
         * MathJax formula: \Delta \lambda_k = - \eta (y_k - \hat{y_k}) \hat{y_k} (1 - \hat{y_k})
         */
        for (size_t k = 0; k < m_num_output_node; ++k) {
            Real bias_delta =-(label[k] - outputLayer[k]->getNodeValue())* outputLayer[k]->getNodeValue() * (1.0 - outputLayer[k]->getNodeValue());

            outputLayer[k]->setNodeBiasDelta(outputLayer[k]->getNodeBiasDelta() + bias_delta);
        }

        /**
         * Calculate 'weight_delta'(the weight correction value)
         * from the j-th node in the hidden layer to the k-th node in the output layer
         * MathJax formula: \Delta v_{jk} = \eta ( y_k - \hat{y_k} ) \hat{y_k} ( 1 - \hat{y_k} ) h_j
         */
        for (size_t j = 0; j < m_num_hiden_node; ++j) {
            for (size_t k = 0; k < m_num_output_node; ++k) {
                Real weight_delta =(label[k] - outputLayer[k]->getNodeValue())* outputLayer[k]->getNodeValue() * (1.0 - outputLayer[k]->getNodeValue())* hiddenLayer[j]->getNodeValue();

                hiddenLayer[j]->setNodeWeightDelta(k,hiddenLayer[j]->getNodeWeightDelta(k) + weight_delta);
            }
        }

        /**
         * Calculate 'bias_delta'(the bias correction value)
         * of the j-th node in the hidden layer
         * MathJax formula: \Delta \beta_j = - \eta \sum_k ( y_k - \hat{y_k} ) \hat{y_k} ( 1 - \hat{y_k} ) v_{jk} h_j ( 1 - h_j )
         */
        for (size_t j = 0; j < m_num_hiden_node; ++j) {
            Real bias_delta = 0.;
            for (size_t k = 0; k < m_num_output_node; ++k) {
                bias_delta +=
                    ( - (label[k] - outputLayer[k]->getNodeValue())
                    * outputLayer[k]->getNodeValue() * (1.0 - outputLayer[k]->getNodeValue())
                    * hiddenLayer[j]->getNodeWeight(k));
            }
            bias_delta *=(hiddenLayer[j]->getNodeValue() * (1.0 - hiddenLayer[j]->getNodeValue()));

            hiddenLayer[j]->setNodeBiasDelta(hiddenLayer[j]->getNodeBiasDelta() + bias_delta);
        }

        /**
         * Calculate 'weight_delta'(the weight correction value)
         * from the i-th node in the input layer to the j-th node in the hidden layer
         * MathJax formula: \Delta w_{ij} = \eta \sum_k ( y_k - \hat{y_k} ) \hat{y_k} ( 1 - \hat{y_k} ) v_{jk} h_j ( 1 - h_j ) x_i
         */
        for (size_t i = 0; i < m_num_input_node; ++i) {
            for (size_t j = 0; j < m_num_hiden_node; ++j) {
                Real weight_delta = 0.;
                for (size_t k = 0; k < m_num_output_node; ++k) {
                    weight_delta +=
                        ((label[k] - outputLayer[k]->getNodeValue())
                        * outputLayer[k]->getNodeValue() * (1.0 - outputLayer[k]->getNodeValue())
                        * hiddenLayer[j]->getNodeWeight(k));
                }
                weight_delta *=
                    (hiddenLayer[j]->getNodeValue() * (1.0 - hiddenLayer[j]->getNodeValue())
                    * inputLayer[i]->getNodeValue());

                inputLayer[i]->setNodeWeightDelta(j,inputLayer[i]->getNodeWeightDelta(j) + weight_delta);
            }
        }
    }

    bool Net::train(const std::vector<Sample>& trainDataSet) {
        for (size_t epoch = 0; epoch <= m_max_epoch; ++epoch) {

            grad_zero();

            Real max_loss = 0.;

            for (const auto& trainSample : trainDataSet) {

                // Load trainSample's feature into the network
                for (size_t i = 0; i < m_num_input_node; ++i) {
                    inputLayer[i]->setNodeValue(trainSample.getSampleFeature(i));
                }
                    
                

                // Forward propagation
                forward();

                // Calculate 'loss'
                Real loss = calculateLoss(trainSample.getSampleLabels());
                max_loss = std::max(max_loss, loss);

                // Back propagation
                backward(trainSample.getSampleLabels());

            }

            updateLearningRate(max_loss);

            m_his_loss.push_back(max_loss);

            // Deciding whether to stop training
            if (max_loss < m_threshold) {
                std::cout << "Training SUCCESS in  "  << epoch << " epochs." << std::endl;
                std::cout << "Final maximum error(loss): " << max_loss << std::endl;
                return true;
            }
            else if (epoch % 1000 == 0) {
                std::cout << "#epoch - max_loss: " << epoch << "  " << max_loss << "  " << "learning rate :"<<getLearningRate() << std::endl;;
            }

            // Revise 'weight' and 'bias' of each node
            revise(trainDataSet.size());
        }

        std::cout << "Failed within " << m_max_epoch << " epoch." << std::endl;

        return false;
    }

    void Net::revise(size_t batch_size) {

        auto batch_size_Real = (Real)batch_size;

        for (size_t i = 0; i < m_num_input_node; ++i) {
            for (size_t j = 0; j < m_num_hiden_node; ++j) {

                /**
                 * Revise 'weight' according to 'weight_delta'(the weight correction value)
                 * from the i-th node in the input layer to the j-th node in the hidden layer
                 */
                inputLayer[i]->setNodeWeight(j,inputLayer[i]->getNodeWeight(j) + (m_learning_rate * inputLayer[i]->getNodeWeightDelta(j) / batch_size_Real));

            }
        }

        for (size_t j = 0; j < m_num_hiden_node; ++j) {

            /**
             * Revise 'bias' according to 'bias_delta'(the bias correction value)
             * of the j-th node in the hidden layer
             */
            hiddenLayer[j]->setNodeBias(hiddenLayer[j]->getNodeBias() + (m_learning_rate * hiddenLayer[j]->getNodeBiasDelta() / batch_size_Real));

            for (size_t k = 0; k < m_num_output_node; ++k) {

                /**
                 * Revise 'weight' according to 'weight_delta'(the weight correction value)
                 * from the j-th node in the hidden layer to the k-th node in the output layer
                 */
                hiddenLayer[j]->setNodeWeight(k,hiddenLayer[j]->getNodeWeight(k) + (m_learning_rate * hiddenLayer[j]->getNodeWeightDelta(k) / batch_size_Real));

            }
        }

        for (size_t k = 0; k < m_num_output_node; ++k) {

            /**
             * Revise 'bias' according to 'bias_weight'(the bias correction value)
             * of the k-th node in the output layer
             */
            outputLayer[k]->setNodeBias(outputLayer[k]->getNodeBias() + (m_learning_rate * outputLayer[k]->getNodeBiasDelta() / batch_size_Real));

        }
    }

    Sample Net::predict(const std::vector<Real>& feature) {

        // load sample into the network
        for (size_t i = 0; i < m_num_input_node; ++i)
            inputLayer[i]->setNodeValue( feature[i]);

        forward();

        std::vector<Real> label(m_num_output_node);
        for (size_t k = 0; k < m_num_output_node; ++k)
            label[k] = outputLayer[k]->getNodeValue();

        Sample pred = Sample(feature, label);
        return pred;
    }

    std::vector<Sample> Net::predict(const std::vector<Sample>& predictDataSet) {
        std::vector<Sample> predSet;

        for (auto& sample : predictDataSet) {
            Sample pred = predict(sample.getSampleFeatures());
            predSet.push_back(pred);
        }

        return predSet;
    }

    Node::Node(size_t nextLayerSize) {
        m_weight.resize(nextLayerSize);
        m_weight_delta.resize(nextLayerSize);
    }

    Sample::Sample() = default;

    Sample::Sample(const std::vector<Real>& feature, const std::vector<Real>& label) :m_feature(feature), m_label(label) {

    }

    void Sample::display() {
        std::cout << "input : ";
        for (auto& x : m_feature) {
            std::cout << x << " ";
        }
        std::cout << std::endl;
        std::cout << "output: ";
        for (auto& y : m_label) {
            std::cout << y;
        }
        std::cout << std::endl;
    }
}

