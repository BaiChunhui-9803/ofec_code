#ifndef OFEF_MIN_MAX_HEAP
#define OFEC_MIN_MAX_HEAP

#include <queue>

namespace ofec {
	template<typename T>
	class MedianHeap {
	private:
		std::priority_queue<T, std::vector<T>, std::greater<T>> m_min_heap;
		std::priority_queue<T, std::vector<T>, std::less<T>> m_max_heap;

		void fixChaos();

	public:
		bool empty() const { return m_min_heap.empty() && m_max_heap.empty(); }
		T median();
		void insert(const T &val);
		void clear();
	};

	template<typename T>
	T MedianHeap<T>::median() {
		if (m_min_heap.size() == m_max_heap.size())
			return (m_min_heap.top() + m_max_heap.top()) / 2;
		else if (m_max_heap.size() > m_min_heap.size())
			return m_max_heap.top();
		else
			return m_min_heap.top();
	}

	template<typename T>
	void MedianHeap<T>::insert(const T &val) {
		if (empty())
			m_min_heap.push(val);
		else {
			if (val < median())
				m_max_heap.push(val);
			else
				m_min_heap.push(val);
		}
		fixChaos();
	}

	template<typename T>
	void MedianHeap<T>::fixChaos() {
		if (m_max_heap.size() > m_min_heap.size() + 1) {
			m_min_heap.push(m_max_heap.top());
			m_max_heap.pop();
		}
		else if (m_min_heap.size() > m_max_heap.size() + 1){
			m_max_heap.push(m_min_heap.top());
			m_min_heap.pop();
		}
	}

	template<typename T>
	void MedianHeap<T>::clear() {
		m_min_heap = std::priority_queue<T, std::vector<T>, std::greater<T>>();
		m_max_heap = std::priority_queue<T, std::vector<T>, std::less<T>>();
	}
}

#endif // !OFEF_MIN_MAX_HEAP
