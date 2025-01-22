#ifndef OFEC_DPC_H
#define OFEC_DPC_H

// Clustering by fast search and find of density peaks
// Science 27 June 2014: 
// Vol. 344 no. 6191 pp. 1492-1496 
// DOI: 10.1126/science.1242072
// http://www.sciencemag.org/content/344/6191/1492.full// 
//
// Code Author: Eric Yuan
// Blog: http://eric-yuan.me
// You are FREE to use the following code for ANY purpose.
//
// Have fun with it

#include "../../core/algorithm/population.h"

namespace ofec {
	class DPC {
	public:
		struct Result {
			struct Infor {
				int clst_no;	//cluster No. 
				bool is_center, is_core;
			};
			std::vector<Infor> info;
			int num_clst;
		};
		enum LocalDensity { kGaussianKernel = 0, kCutOffKernel = 1 };

		template<typename TInd>
		void setData(Population<TInd> &pop, Problem *pro, int space = 0);
		template<typename TInd>
		void setData(const std::vector<std::unique_ptr<TInd>> &inds, Problem *pro, int space = 0);
		template<typename TInd>
		void setData(const std::vector<TInd> &inds, Problem *pro, int space = 0);

		void setMinRhoAndMinDelta(double, double);
		void setDisMatrix(const std::vector<std::vector<double> > &);
		void clustering(Result &result);
		void clustering();
		void clustering(double ratio_rho, double ratio_delta);
		std::vector<std::vector<size_t>> clusters() const;
		const std::vector<double>& getRho() const { return mv_rho; }
		const std::vector<double>& getDelta() const { return mv_delta; }

	private:
		void getdc();
		void getLocalDensity(LocalDensity LocalDensityVersion = kGaussianKernel);
		void getDistanceToHigherDensity();
		void getClusterNum();
		void assign();
		void getHalo(Result &result);
		void autoSetMinRhondMinDelta(double val);
		void setMinRhoNdMinDelta(double ratio_rho, double ratio_delta);

	private:
		int m_num_samples;
		std::vector<double> mv_rho;
		std::vector<double> mv_delta;
		std::vector<int> mv_cl;
		std::vector<int> mv_nneigh;
		std::vector<int> mv_icl;
		std::vector<int> mv_halo;
		double m_dc, m_rho_min, m_deltamin;
		int m_NCLUST;
		std::vector<std::vector<double> > mv_dis;
	};
	
	template<typename TInd>
	void DPC::setData(Population<TInd> &pop, Problem *pro, int space) {
		m_num_samples = pop.size();
		mv_rho.assign(m_num_samples, 0);
		mv_delta.assign(m_num_samples, std::numeric_limits<double>::max());
		mv_cl.assign(m_num_samples, -1);
		mv_nneigh.assign(m_num_samples, -1);
		mv_halo.assign(m_num_samples, 0);
		mv_icl.clear();
		m_NCLUST = 0;
		m_dc = 0.0;
		m_rho_min = 1;
		m_deltamin = 0.1;
		mv_dis.resize(m_num_samples);
		for (size_t i = 0; i < m_num_samples; ++i) {
			mv_dis[i].resize(m_num_samples);
		}
		for (size_t i = 0; i < m_num_samples; ++i) {
			mv_dis[i][i] = 0;
			for (size_t j = i + 1; j < m_num_samples; ++j) {
				if (space == 0) {
					mv_dis[i][j] = mv_dis[j][i] = pop[i].variableDistance(pop[j], pro);
				}
				else {
					mv_dis[i][j] = mv_dis[j][i] = pop[i].objectiveDistance(pop[j]);
				}
			}
		}
	}

	template<typename TInd>
	void DPC::setData(const std::vector<std::unique_ptr<TInd>> &inds, Problem *pro, int space) {
		m_num_samples = inds.size();
		mv_rho.assign(m_num_samples, 0);
		mv_delta.assign(m_num_samples, std::numeric_limits<double>::max());
		mv_cl.assign(m_num_samples, -1);
		mv_nneigh.assign(m_num_samples, -1);
		mv_halo.assign(m_num_samples, 0);
		mv_icl.clear();
		m_NCLUST = 0;
		m_dc = 0.0;
		m_rho_min = 1;
		m_deltamin = 0.1;
		mv_dis.resize(m_num_samples);
		for (size_t i = 0; i < m_num_samples; ++i) {
			mv_dis[i].resize(m_num_samples);
		}
		for (size_t i = 0; i < m_num_samples; ++i) {
			mv_dis[i][i] = 0;
			for (size_t j = i + 1; j < m_num_samples; ++j) {
				if (space == 0) {
					mv_dis[i][j] = mv_dis[j][i] = inds[i]->variableDistance(*inds[j], pro);
				}
				else {
					mv_dis[i][j] = mv_dis[j][i] = inds[i]->objectiveDistance(*inds[j]);
				}
			}
		}
	}

	template<typename TInd>
	void DPC::setData(const std::vector<TInd> &inds, Problem *pro, int space) {
		m_num_samples = inds.size();
		mv_rho.assign(m_num_samples, 0);
		mv_delta.assign(m_num_samples, std::numeric_limits<double>::max());
		mv_cl.assign(m_num_samples, -1);
		mv_nneigh.assign(m_num_samples, -1);
		mv_halo.assign(m_num_samples, 0);
		mv_icl.clear();
		m_NCLUST = 0;
		m_dc = 0.0;
		m_rho_min = 1;
		m_deltamin = 0.1;
		mv_dis.resize(m_num_samples);
		for (size_t i = 0; i < m_num_samples; ++i) {
			mv_dis[i].resize(m_num_samples);
		}
		for (size_t i = 0; i < m_num_samples; ++i) {
			mv_dis[i][i] = 0;
			for (size_t j = i + 1; j < m_num_samples; ++j) {
				if (space == 0) {
					mv_dis[i][j] = mv_dis[j][i] = inds[i].variableDistance(inds[j], pro);
				}
				else {
					mv_dis[i][j] = mv_dis[j][i] = inds[i].objectiveDistance(inds[j]);
				}
			}
		}
	}
}

#endif // !OFEC_DPC_H