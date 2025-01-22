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

#include <cstddef>
#include <memory>
#include <mutex>
#include <queue>


#include "sort_sol.h"


namespace n2 {

    class BaseNeighborSelectingPoliciesSol {
    public:
        BaseNeighborSelectingPoliciesSol() {}
        virtual ~BaseNeighborSelectingPoliciesSol() {}

        virtual void Select(size_t m, 
            bool select_nn, 
            ofec::Environment *env,
            std::priority_queue<FurtherSolFirst>& result) = 0;
    };

    class NaiveNeighborSelectingPoliciesSol : public BaseNeighborSelectingPoliciesSol {
    public:
        NaiveNeighborSelectingPoliciesSol() {}
        ~NaiveNeighborSelectingPoliciesSol() override {}
        void Select(size_t m, bool select_nn,
            ofec::Environment *env,
            std::priority_queue<FurtherSolFirst>& result) override;
    };



    class HeuristicNeighborSelectingPoliciesSol : public BaseNeighborSelectingPoliciesSol {
    public:
        HeuristicNeighborSelectingPoliciesSol() : save_remains_(false) {}
        HeuristicNeighborSelectingPoliciesSol(bool save_remain) : save_remains_(save_remain) {}
        ~HeuristicNeighborSelectingPoliciesSol() override {}
        /**
         * Returns selected neighbors to result
         * (analagous to SELECT-NEIGHBORS-HEURISTIC in Yu. A. Malkov's paper.)
         *
         * select_nn: if true, select 0.25 * m nearest neighbors to result without applying the heuristic algorithm
         */
        void Select(size_t m,  bool select_nn,
            ofec::Environment *env,
            std::priority_queue<FurtherSolFirst>& result) override;


    private:
        bool save_remains_;
    public:

    };


    extern void selectNeighbors(std::vector<EdgeInfo>& neighbors,
        int level, int maxM, int maxM0, bool is_naive_, 
        BaseNeighborSelectingPoliciesSol* select_policy, ofec::Environment *env);
    



} // namespace n2
