#pragma once
#include"peak.h"
#include <numeric> 
#include "../../../../../core/definition.h"

namespace ofec {

	class MultiplePeaksModel
	{
	private:
		std::random_device rd;
		std::mt19937 gen;
		std::uniform_real_distribution<double> uni_domain;
		std::normal_distribution<double> nor_domain;
		std::uniform_real_distribution<double> uni_rand_shapeRange;
		std::uniform_real_distribution<double> uni_rand_radiusRange;
		std::uniform_real_distribution<double> uni_rand_heightRange;


	public:
		int num_var;
		int num_peaks;
		int num_optimal;
		std::string m_topology;
		int shapeHeightCorrelation;
		std::string m_peakShape;
		bool m_rotation;

		std::vector<std::unique_ptr<Peak>> m_peaks;

		MultiplePeaksModel(const int vars = 2, const int peaksnum = 4, const int optimanum = 1, const std::string topology = "random", const int correlation = 0, const std::string peakShape = "ellipse", const bool rotation = true);

		void randomUniformPeaks(int num_peaks);
		void createInstance();
		void makeFunnel();
		void clusteredPeaks(int num_peaks);

		double G(Real* x, const std::unique_ptr<Peak>& peak);
		int getActivePeak(Real* x);//õxڵķ

		double getGvalue(Real* x);

		int getLocalOptima();

		void createInstanceWithExactNumberOfOptima();// int numberOptima = 1, std::string topology = "random", int shapeHeightCorrelation = 0);

	};

}