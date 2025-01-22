// Copyright 2017 Kakao Corp. <http://www.kakaocorp.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once



//#include "../"
#include "hnsw_node.h"
#include "../functional.h"
#include "../../core/problem/solution.h"
#include "../../core/problem/problem.h"
namespace n2
{




    class DataSol : public Data {
    public:

        void bindData(ofec::SolBase* sol) {
            m_sol = sol;
            //m_pro = pro;
        }
        const ofec::SolBase* getSol()const {
            return m_sol;
        }



    protected:
        ofec::SolBase* m_sol = nullptr;

    };

    class SolDistance {
    protected:
        //	ofec::Problem* m_pro = nullptr;
    public:
        //void bindPro(ofec::Problem* pro) {
        //	m_pro = pro;
        //}
        //inline float operator()(const float* v1, const float* v2, size_t qty) const {
        //	std::vector<double> a1(qty);
        //	std::vector<double> a2(qty);
        //	for (int idx(0); idx < qty; ++idx) {
        //		a1[idx] = v1[idx];
        //		a2[idx] = v2[idx];
        //	}
        //	return ofec::euclideanDistance(a1.begin(), a1.end(), a2.begin());
        //	//   Eigen::Map<const Eigen::VectorXf, Eigen::Unaligned> p(v1, qty, 1), q(v2, qty, 1);
        //	   //return (p - q).squaredNorm();
        //}
        //inline float operator()(const HnswNode* n1, const HnswNode* n2, size_t qty) const {
        //	auto data1 = dynamic_cast<const DataSol*>(n1->getOriginData());
        //	auto sol1 = data1->getSol();
        //	auto sol2 = dynamic_cast<const DataSol*>(n2->getOriginData())->getSol();
        //	return sol1->variableDistance(*sol2, m_pro);
        ////	return (*this)(n1->GetData(), n2->GetData(), qty);
        //}



        inline float operator()(const HnswNode* n1, const HnswNode* n2, size_t qty) const {
            return 1.0;
            //        return (*this)(n1->GetData(), n2->GetData(), qty);
        }
    };
    class SolDis {

    protected:
        ofec::Problem* m_pro = nullptr;
    public:

        void bindPro(ofec::Problem* pro) {
            m_pro = pro;
        }

        //inline float operator()(const float* v1, const float* v2, size_t qty) const {
        //    std::vector<double> a1(qty);
        //    std::vector<double> a2(qty);
        //    for (int idx(0); idx < qty; ++idx) {
        //        a1[idx] = v1[idx];
        //        a2[idx] = v2[idx];
        //    }
        //    return ofec::euclideanDistance(a1.begin(), a1.end(), a2.begin());
        // //   Eigen::Map<const Eigen::VectorXf, Eigen::Unaligned> p(v1, qty, 1), q(v2, qty, 1);
        //    //return (p - q).squaredNorm();
        //}
        inline float operator()(const HnswNode* n1, const HnswNode* n2, size_t qty) const {
            return 1.0;
            //        return (*this)(n1->GetData(), n2->GetData(), qty);
        }
    };




    class L2Distance2 {

    protected:
        ofec::Problem* m_pro = nullptr;
    public:

        void bindPro(ofec::Problem* pro) {
            m_pro = pro;
        }

        inline float operator()(const float* v1, const float* v2, size_t qty) const {
            std::vector<double> a1(qty);
            std::vector<double> a2(qty);
            for (int idx(0); idx < qty; ++idx) {
                a1[idx] = v1[idx];
                a2[idx] = v2[idx];
            }
            return ofec::euclideanDistance(a1.begin(), a1.end(), a2.begin());
         //   Eigen::Map<const Eigen::VectorXf, Eigen::Unaligned> p(v1, qty, 1), q(v2, qty, 1);
            //return (p - q).squaredNorm();
        }
        //inline float operator()(const HnswNode* n1, const HnswNode* n2, size_t qty) const {
        //    return 1.0;
        //    //        return (*this)(n1->GetData(), n2->GetData(), qty);
        //}
    };


    class DistanceBase {
    public: 
        virtual void bindPro(ofec::Problem* pro) {}
    };

    class L2Distance: public DistanceBase {

    protected:
        ofec::Problem* m_pro = nullptr;
    public:

        virtual void bindPro(ofec::Problem* pro) override {
            m_pro = pro;
        }

        //inline float operator()(const float* v1, const float* v2, size_t qty) const {
        //    std::vector<double> a1(qty);
        //    std::vector<double> a2(qty);
        //    for (int idx(0); idx < qty; ++idx) {
        //        a1[idx] = v1[idx];
        //        a2[idx] = v2[idx];
        //    }
        //    return ofec::euclideanDistance(a1.begin(), a1.end(), a2.begin());
        // //   Eigen::Map<const Eigen::VectorXf, Eigen::Unaligned> p(v1, qty, 1), q(v2, qty, 1);
        //    //return (p - q).squaredNorm();
        //}
        //inline float operator()(const HnswNode* n1, const HnswNode* n2, size_t qty) const {
        //    return 1.0;
        //    //        return (*this)(n1->GetData(), n2->GetData(), qty);
        //}

              inline float operator()(const HnswNode* n1, const HnswNode* n2, size_t qty) const {
        	auto data1 = dynamic_cast<const DataSol*>(n1->getOriginData());
        	auto sol1 = data1->getSol();
        	auto sol2 = dynamic_cast<const DataSol*>(n2->getOriginData())->getSol();
        	return sol1->variableDistance(*sol2, m_pro);
        //	return (*this)(n1->GetData(), n2->GetData(), qty);
        }

    };

    class AngularDistance : public DistanceBase {
    public:
        inline float operator()(const float* v1, const float* v2, size_t qty) const {
            double value(0);
            // std::vector<double> a1(qty);
             //std::vector<double> a2(qty);
            for (int idx(0); idx < qty; ++idx) {
                value += v1[idx] * v2[idx];
            }

            return 1.0 - value;

            //  Eigen::Map<const Eigen::VectorXf, Eigen::Unaligned> p(v1, qty, 1), q(v2, qty, 1);
            //return 1.0 - p.dot(q);
        }
        inline float operator()(const HnswNode* n1, const HnswNode* n2, size_t qty) const {
            return (*this)(n1->GetData(), n2->GetData(), qty);
        }
    };

    class DotDistance : public DistanceBase {
    public:
        inline float operator()(const float* v1, const float* v2, size_t qty) const {
            double value(0);
            // std::vector<double> a1(qty);
             //std::vector<double> a2(qty);
            for (int idx(0); idx < qty; ++idx) {
                value += v1[idx] * v2[idx];
            }

            return -value;
            //   Eigen::Map<const Eigen::VectorXf, Eigen::Unaligned> p(v1, qty, 1), q(v2, qty, 1);
         //   return -p.dot(q);
        }
        inline float operator()(const HnswNode* n1, const HnswNode* n2, size_t qty) const {
            return (*this)(n1->GetData(), n2->GetData(), qty);
        }
    };

} // namespace n2
