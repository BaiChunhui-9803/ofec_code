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

#include <mutex>
#include <vector>
#include <unordered_set>
#include <shared_mutex>
#include "../../../core/problem/solution.h"
#include "../../../core/environment/environment.h"


namespace n2 {

    class BaseNeighborSelectingPoliciesSol;
    class HnswSolutionNode;

    class EdgeInfo {

    protected:
        HnswSolutionNode* node_;
        ofec::Real distance_;

    public:
        //HnswSolutionNode* m_node = nullptr;
        //double m_distance = 0;


        inline ofec::Real GetDistance() const { return distance_; }
        inline HnswSolutionNode* GetNode() const { return node_; }

        void set(HnswSolutionNode* node, ofec::Real dis) {
            node_ = node;
            distance_ = dis;
        }


        EdgeInfo() = default;
        EdgeInfo(HnswSolutionNode* node, double dis) : node_(node), distance_(dis){};
        EdgeInfo(const EdgeInfo& other) {
            node_ = other.node_;
            distance_ = other.distance_;
        };

        EdgeInfo& operator=(const EdgeInfo& other) {
            node_ = other.node_;
            distance_ = other.distance_;
            return *this;
        };

        bool operator==(const EdgeInfo& other) const {
            if (node_ == other.node_)
                return true;
            return false;
        };

        bool operator<(const EdgeInfo& other) {
            if (node_ < other.node_)
                return true;
            else if (node_ == other.node_)
                return true;

            return false;
        };

    };


    class HnswSolutionNode {
    public:

     //   using EdgeInfo = std::pair<HnswSolutionNode*, double> ;

        explicit HnswSolutionNode(int id, ofec::SolutionBase* data, int level, size_t max_m, size_t max_m0)
            : id_(id), level_(level), max_m_(max_m), max_m0_(max_m0), friends_at_layer_(level + 1) {
            
            if (data == nullptr) {
                data_ = nullptr;
            }
            else {
                data_ = data->shared_from_this();
            }
            for (int i = 1; i <= level; ++i) {
                friends_at_layer_[i].reserve(max_m_ + 1);
            }
            friends_at_layer_[0].reserve(max_m0_ + 1);
        }

        void setSolution(ofec::SolutionBase* sol) {
            data_ = sol->shared_from_this();
        }


    //    void CopyHigherLevelLinksToOptIndex(char* mem_offset, uint64_t memory_per_node_higher_level) const;
   //     void CopyDataAndLevel0LinksToOptIndex(char* mem_offset, int higher_level_offset) const;

        inline int GetId() const { return id_; }
        inline int GetLevel() const { return level_; }
        inline int& level() { return level_; }
        inline size_t GetMaxM() const { return max_m_; }
        inline size_t GetMaxM0() const { return max_m0_; }


        const ofec::SolutionBase& getOriginData() const{
            return *data_.get();
        }
        ofec::SolutionBase& getOriginData() {
            return *data_.get();
        }

        inline void getFriends(std::vector<EdgeInfo>& neighbors, int level)const {
            std::shared_lock lock(mutex_);
            neighbors = friends_at_layer_[level];
        }

        inline const std::vector<EdgeInfo>& getFriends(int level) {
           //std::shared_lock lock(mutex_);
            return friends_at_layer_[level];
        }


        
        void link(HnswSolutionNode* target, double dis, 
            int level, 
            bool is_naive_,
            BaseNeighborSelectingPoliciesSol* select_policy,
            ofec::Environment *env);
        void linkSingleThread(HnswSolutionNode* target, double dis, int level, 
            bool is_naive_,
            BaseNeighborSelectingPoliciesSol* select_policy,
            ofec::Environment *env);

        const std::vector<std::vector<EdgeInfo>>& getTotalFriends()const {
            return friends_at_layer_;
        }
        std::vector<std::vector<EdgeInfo>>& getTotalFriends() {
            return friends_at_layer_;
        }


        inline void SetFriends(int level,const std::vector<EdgeInfo>& new_friends) {
            std::unique_lock lock(mutex_);
            friends_at_layer_[level] = new_friends;
        }
        inline void setFriendsSingleThread(int level, const std::vector<EdgeInfo>& new_friends) {
           // std::unique_lock lock(mutex_);
            friends_at_layer_[level] = new_friends;
        }
     //   void CopyLinksToOptIndex(char* mem_offset, int level) const;


        void getNeighborIds(std::vector<std::vector<int>>& neighbors) {
            neighbors.resize(friends_at_layer_.size());
            for (int idx(0); idx < neighbors.size(); ++idx) {
                neighbors[idx].resize(friends_at_layer_[idx].size());
                for (int idy(0); idy < neighbors[idx].size(); ++idy) {
                    neighbors[idx][idy] = friends_at_layer_[idx][idy].GetNode()->GetId();
                }
            }
        }


        std::shared_mutex& mutex() {
            return mutex_;
        }
    private:
        int id_;
        std::shared_ptr<ofec::SolutionBase> data_;
        int level_;
        const size_t max_m_;
        const size_t max_m0_;

        std::vector<std::vector<EdgeInfo>> friends_at_layer_;
        mutable std::shared_mutex mutex_;
    };

} // namespace n2
namespace std {
    template<> struct hash<n2::EdgeInfo>
    {
        std::size_t operator()(const n2::EdgeInfo& p) const noexcept
        {
            return hash<int>{}(p.GetNode()->GetId());
        }

        //  hash() = default;
    };
}
