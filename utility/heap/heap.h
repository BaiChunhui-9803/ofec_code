#ifndef OFEC_MIN_MAX_HEAP_H
#define OFEC_MIN_MAX_HEAP_H

#include <vector>

namespace ofec {
	template <typename T>
	class Heap {
	private:
		std::vector<T> m_data;
		std::vector<size_t> m_seq;
		std::vector<int> m_pos;
	public:
		Heap(const std::vector<T> &data);
		void updateValue(size_t id, const T &val);
		void remove(size_t id);
		void add(const T &val);
		size_t top();
	private:
		void build();
		void checkUp(size_t id);
		void checkDown(size_t id);
	};
	
	template<typename T>
	void Heap<T>::build() {
		m_pos.resize(m_data.size());
		m_seq.clear();
		for (size_t i = 0; i < m_data.size(); ++i) {
			m_pos[i] = m_seq.size();
			m_seq.push_back(i);
			checkUp(i);
		}
	}
	
	template<typename T>
	Heap<T>::Heap(const std::vector<T> &data) : m_data(data) {
		build();
	} 

	template<typename T>
	void Heap<T>::updateValue(size_t id, const T &val) {
		T pre_val = m_data[id];
		m_data[id] = val;
		if (pre_val > val)
			checkUp(id);
		else if (pre_val < val) 
			checkDown(id);
	}
	
	template<typename T>
	void Heap<T>::remove(size_t id) {
		if (m_pos[id] == -1)
			return;
		size_t id_back = m_seq.back();
		m_seq[m_pos[id]] = id_back;
		m_pos[id_back] = m_pos[id];
		m_pos[id] = -1;
		m_seq.pop_back();
		checkDown(id_back);
	}

	template<typename T>
	void Heap<T>::add(const T &val) {
		m_pos.push_back(m_seq.size());
		m_seq.push_back(m_data.size());
		m_data.push_back(val);
		checkUp(m_data.size() - 1);
	}
	
	template<typename T>
	size_t Heap<T>::top() {
		return m_seq.front();
	}
	
	template<typename T>
	inline void Heap<T>::checkUp(size_t id) {
		int child = m_pos[id];
		int parent = (child - 1) / 2;
		while (child > 0) {
			if (m_data[m_seq[parent]] <= m_data[m_seq[child]])
				break;
			auto temp = m_seq[parent];
			m_seq[parent] = m_seq[child];
			m_seq[child] = temp;
			m_pos[m_seq[child]] = child;
			m_pos[m_seq[parent]] = parent;
			child = parent;
			parent = (child - 1) / 2;
		}
	}
	
	template<typename T>
	inline void Heap<T>::checkDown(size_t id) {
		int parent = m_pos[id];
		int child = 2 * parent + 1;
		while (child < m_seq.size()) {
			if (child + 1 < m_seq.size() && m_data[m_seq[child]] > m_data[m_seq[child + 1]])
				child++;
			if (m_data[m_seq[parent]] <= m_data[m_seq[child]])
				break;
			auto temp = m_seq[parent];
			m_seq[parent] = m_seq[child];
			m_seq[child] = temp;
			m_pos[m_seq[child]] = child;
			m_pos[m_seq[parent]] = parent;
			parent = child;
			child = 2 * parent + 1;
		}
	}
}

#endif // !OFEC_MIN_MAX_HEAP_H
