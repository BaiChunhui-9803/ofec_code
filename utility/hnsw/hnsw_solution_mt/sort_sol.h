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

#include "hnsw_sol_node.h"

namespace n2 {

    class FurtherSolFirst:public EdgeInfo {
    public:
       // FurtherSolFirst(HnswSolutionNode* node, float distance) : node_(node), distance_(distance) {}
        FurtherSolFirst() = default;
        FurtherSolFirst(HnswSolutionNode* node, double dis) : EdgeInfo(node,dis){};
        FurtherSolFirst(const EdgeInfo& a) :EdgeInfo(a) {}
        bool operator< (const FurtherSolFirst& n) const {
            return (distance_ < n.distance_);
        }

    };

    class CloserSolFirst :public EdgeInfo {
    public:


        CloserSolFirst() = default;
        CloserSolFirst(HnswSolutionNode* node, double dis) : EdgeInfo(node, dis) {};
        CloserSolFirst(const EdgeInfo& a) :EdgeInfo(a) {}
        //CloserSolFirst(HnswSolutionNode* node, float distance) : node_(node), distance_(distance) {}
        //inline double GetDistance() const { return distance_; }
        //inline HnswSolutionNode* GetNode() const { return node_; }
        bool operator< (const CloserSolFirst& n) const {
            return (distance_ > n.distance_);
        }
    //private:
    //    HnswSolutionNode* node_;
    //    double distance_;
    };

} // namespace n2
