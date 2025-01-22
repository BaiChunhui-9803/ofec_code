#include"memory_record.h"


//void ofec::MemeryRecord::initialize(Problem *pro, double eps_distance)
//{
//	m_memory.clear();
//	m_problem.get() = pro;
//	m_eps_distance = eps_distance;
//}

void ofec::MemeryRecord::inputFile(const std::string& filename)
{
	std::unique_ptr<SolBase> sol;
	std::ifstream in(filename);
	if (in.is_open()) {
		//while (m_problem->inSolution(in, sol)) {
		//	m_memory.emplace_back(sol.release());
		//}
	}

	in.close();
}

void ofec::MemeryRecord::outputFile(const std::string& filename)
{
	std::ofstream out(filename);
	//for (auto& it : m_memory) {
	//	m_problem->outSolution(out, it);
	//}
	out.close();
}



void ofec::MemeryRecord::append(std::unique_ptr<SolBase>& p) {
	m_com_relation.resize(m_memory.size());
	std::fill(m_com_relation.begin(), m_com_relation.end(), 0);
	bool p_active(true);
	p->objective().resize(1);
	// store fitness value
	p->objective()[0] = m_fitness_fun(m_problem.get(), *p);
	//	m_memory.emplace_back(p.release());
	for (int mid(0); mid < m_memory.size(); ++mid) {
		if (m_memory[mid]->objective()[0] > p->objective()[0]) {
			if (m_memory[mid]->variableDistance(*p, m_problem.get()) < m_eps_distance) {
				p_active = false;
			}
			m_com_relation[mid] = 1;
		}
	}
	if (p_active) {
		bool erase_memory(false);
		//		m_active.resize(m_memory.size());
		//		std::fill(m_active.begin(), m_active.end(), true);
		for (int mid(0); mid < m_memory.size(); ++mid) {
			if (m_com_relation[mid] == 0) {
				if (p->objective()[0] > m_memory[mid]->objective()[0]) {
					m_com_relation[mid] = -1;
					if (m_memory[mid]->variableDistance(*p, m_problem.get()) < m_eps_distance) {
						m_com_relation[mid] = -2;
						erase_memory = true;
					}
				}
			}
		}
		if (erase_memory) {
			int active_idx(0);
			for (int mid(0); mid < m_memory.size(); ++mid) {
				if (m_com_relation[mid] == -2) {
					if (mid != active_idx) {
						swap(m_memory[mid], m_memory[active_idx]);
					}
					++active_idx;
				}
			}
			swap(m_memory[active_idx++], p);
			m_memory.resize(active_idx);
		}
		else {
			m_memory.emplace_back(p.release());
		}
	}
}