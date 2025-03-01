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

#include "hnsw_node.h"

namespace n2 {

HnswNode::HnswNode(int id, const Data* data, int level, size_t max_m, size_t max_m0)
        : id_(id), data_(data), level_(level), max_m_(max_m), max_m0_(max_m0), friends_at_layer_(level+1) {
    for (int i = 1; i <= level; ++i) {
        friends_at_layer_[i].reserve(max_m_ + 1);
    }
    friends_at_layer_[0].reserve(max_m0_ + 1);
}

void HnswNode::CopyHigherLevelLinksToOptIndex(char* mem_offset, uint64_t memory_per_node_higher_level) const {
    char* mem_data = mem_offset;
    for (int level = 1; level <= level_; ++level) {
        CopyLinksToOptIndex(mem_data, level);
        mem_data += memory_per_node_higher_level;
    }
}

void HnswNode::CopyDataAndLevel0LinksToOptIndex(char* mem_offset, int higher_level_offset) const {
    char* mem_data = mem_offset;
    *((int*)(mem_data)) = higher_level_offset;
    mem_data += sizeof(int);
    CopyLinksToOptIndex(mem_data, 0);
    mem_data += (sizeof(int) + sizeof(int)*max_m0_);
    auto& data = data_->GetData();
    for (size_t i = 0; i < data.size(); ++i) {
        *((float*)(mem_data)) = (float)data[i];
        mem_data += sizeof(float);
    }
}

void HnswNode::CopyLinksToOptIndex(char* mem_offset, int level) const {
    char* mem_data = mem_offset;
    const auto& neighbors = friends_at_layer_[level];
    *((int*)(mem_data)) = (int)(neighbors.size());
    mem_data += sizeof(int);
    for (size_t i = 0; i < neighbors.size(); ++i) {
        *((int*)(mem_data)) = (int)neighbors[i]->GetId();
        mem_data += sizeof(int);
    }
}

} // namespace n2
