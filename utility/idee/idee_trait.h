#ifndef   OFEC_UTILITY_IDEE_TRAIT_H
#define   OFEC_UTILITY_IDEE_TRAIT_H


#include "../../utility/function/custom_function.h"



namespace ofec {
	

	struct EAX_TSP_trait {


		static void outputIDEEsolution(std::ostream& out,const std::vector<std::vector<int>>& sols) {
			for (int idx(0); idx < sols.size(); ++idx) {
				out << idx << "\t";
				bool first = true;
				for (auto& it : sols[idx]) {
					if (first)
						first = false;
					else out << ",";
					out << it;
				}
				out << std::endl;
			}
		}

		int m_numPop = 100;

		std::vector<std::vector<int>> m_GenPopSolIds;

		std::vector<std::pair<int, int>> m_runId_genIds;

		std::vector<int> m_fitness;
		std::vector<std::vector<int>> m_sols;

		void clear() {
			m_GenPopSolIds.clear();
			m_runId_genIds.clear();
			m_fitness.clear();
			m_sols.clear();
		}



		void inputPopSolIds(const std::string& filepath) {
			std::stringstream buf;
			UTILITY::readFileBuffer(filepath, buf);
			std::string headline;
			std::getline(buf, headline);
			int generationId(0);
			std::vector<int>  curPop;
			while (buf >> generationId) {
				curPop.resize(m_numPop);
				for (auto& it : curPop) {
					buf >> it;
				}
				m_GenPopSolIds.emplace_back(std::move(curPop));
			}
		}


		void inputRunGenPopSolIds(const std::string& filepath) {
			std::stringstream buf;
			UTILITY::readFileBuffer(filepath, buf);
			std::string headline;
			std::getline(buf, headline);
			std::pair<int, int> curRunGen;
			std::vector<int>  curPop;
			while (buf >> curRunGen.first >> curRunGen.second) {
				curPop.resize(m_numPop);
				for (auto& it : curPop) {
					buf >> it;
				}
				m_runId_genIds.push_back(curRunGen);
				m_GenPopSolIds.emplace_back(std::move(curPop));
			}
		}

		void inputSolInfo(const std::string& filepath) {
			std::stringstream buf;
			UTILITY::readFileBuffer(filepath, buf);

			std::string headline;
			std::getline(buf, headline);
			int solId(0);
			int fitness(0);

			while (buf >> solId >> fitness) {
				m_fitness.push_back(fitness);
			}
		}

		void inputSols(const std::string& filepath) {
			std::stringstream buf;
			UTILITY::readFileBuffer(filepath, buf);
			std::string curLine;
			std::stringstream linebuf;
			std::getline(buf, curLine);
			int solId(0);
			int solutionVec(0);
			std::vector<int> curSol;

			while (buf >> solId) {
				std::getline(buf, curLine);
				linebuf.clear();
				linebuf << curLine;
				curSol.clear();
				while (linebuf >> solutionVec) {
					curSol.push_back(solutionVec);
				}

				m_sols.emplace_back(std::move(curSol));
			}

			/*
			* 			double curVal(0);
				std::string curLine;
				std::getline(buf, curLine);
				std::stringstream tmpbuf;
				tmpbuf << curLine;
				while (tmpbuf >> curVal) {
					cur_data.push_back(curVal);
				}
			*/
		}


		void outputPopSolIds(const std::string& filepath) {
			std::ofstream out(filepath);
			out << "RunId\tGenId\tSolIds" << std::endl;
			for (int idx(0); idx < m_GenPopSolIds.size(); ++idx) {
				out << m_runId_genIds[idx].first << "\t" << m_runId_genIds[idx].second;
				for (auto& it : m_GenPopSolIds[idx]) {
					out << "\t" << it;
				}
				out << std::endl;
			}
			out.close();
		}
		void outputSolInfo(const std::string& filepath) {

			std::ofstream out(filepath);
			out << "ID\tFitness" << std::endl;
			for (int idx(0); idx < m_fitness.size(); ++idx) {
				out << idx << "\t" << m_fitness[idx] << std::endl;
			}
			out.close();
		}


		void outputSols(const std::string& filepath) {
			std::ofstream out(filepath);
			out << "ID\tSolution" << std::endl;
			for (int idx(0); idx < m_sols.size(); ++idx) {
				out << idx;
				for (auto& it : m_sols[idx]) {
					out << "\t" << it;
				}
				out << std::endl;
			}

			out.close();
		}


		void outputIDEESols(const std::string& filepath) {
			std::ofstream out(filepath);
			for (int idx(0); idx < m_sols.size(); ++idx) {
				out << idx << "\t";
				bool first = true;
				for (auto& it : m_sols[idx]) {
					if (first)
						first = false;
					else out << ",";
					out << it;
				}
				out << std::endl;
			}

			out.close();
		}


	};

}


#endif // !  OFEC_UTILITY_IDEE_TRAIT_H
