
#pragma once
#include <vector>

namespace custom_fun {

    class VisitedLazyCountList {
    public:
        VisitedLazyCountList() = default;
        VisitedLazyCountList(int size) :visited(size, 0), totalSize(size){}
        ~VisitedLazyCountList() = default;

        void initialize() {
            mark_ = 1;
            visited.clear();
        }
        void resize(int numSize) {
            visited.resize(numSize);
            totalSize = visited.size();
        }
        void insert(int number = 1) {
            visited.resize(visited.size() + number, 0);
            totalSize = visited.size();
        }
        inline bool Visited(unsigned int index) const { return visited[index] >= mark_; }
        inline bool NotVisited(unsigned int index) const { return visited[index] < mark_; }
        inline void Mark(unsigned int index) {
            visited[index] = mark_;
        }
        inline void VisitedOnce(unsigned int index)const {
            return visited[index] == mark_;
        }
        inline void increaseValue(unsigned int index) {
            ++visited[index];
        }
        inline void Reset() {
            if (++mark_+ totalSize == 0) {
                mark_ = 1;
                std::fill(visited.begin(), visited.end(), 0);
                //memset(visited_, 0, sizeof(unsigned int) * size_);
            }
        }
        //   inline unsigned int* GetVisited() { return visited_; }
        inline unsigned long long GetVisitMark() { return mark_; }

    private:
        std::vector<unsigned long long> visited;
        unsigned long long mark_ = 0;
        unsigned long long totalSize = 0;
        
    };

} // namespace n2
