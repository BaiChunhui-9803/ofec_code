
#ifndef MATLAB_NBN_VISUALIZATION_H
#define MATLAB_NBN_VISUALIZATION_H

#include <vector>
#include <string>
#include <array>


namespace ofec {
	struct MatlabNBNVisualization {

		enum  FigureType {
			kContour2D = 1 << 1,
			kContour2DNBN = 1 << 2,
			kContour3DNBN = 1 << 3,
			kGraph3DNBN = 1 << 4,
			kGraph3DNBNAxis = 1 << 5
		};

		

		std::vector<std::array<double, 3>> m_opt3D_pts_colors;
		std::vector<std::array<double, 3>> m_opt3D_pts;
		std::vector<std::pair<int, int>> m_opt3D_edges;
		std::vector<int> m_marker_idxs;
		int m_markerSize = 10;
		int m_cut = 1000;
		int contour3Level = 500;
		int m_fontSize = 20;
		int m_figure_type = 0;

		int m_color_label = 0;
		int m_axis_label = 1;
		std::string m_filetype = "eps";

		std::string m_filepath;


		void plot(int figureId)const;


		static void DrawMatlabFiguresOutput(int figureId,
			const std::vector<std::array<double, 3>>& opt3D_pts_colors,
			const std::vector<std::array<double, 3>>& opt3D_pts,
			const std::vector<std::pair<int, int>>& opt3D_edges,
			const std::vector<int>& marker_idxs,
			int m_markerSize,
			int m_cut,
			int contour3Level,
			const std::string& filepath,
			const std::string& filetype,
			int fontSize,
			int type);

	};


	struct colorMapper {
		std::vector<std::array<double, 3>> m_dcolor_table= { 
			{0, 0, 0} , {0, 0, 100},{0, 170, 170}, {170, 170, 0}, {255, 255, 255} };
		std::vector<double> m_value_div = {0,0.25,0.5,0.75,1.0};
		
		void init();
		void getColor(double val, std::array<double, 3>& dcolor);

		//fsample[0] = 0.0->csample[0] = (0, 0, 0)
		//	fsample[1] = 0.25->csample[1] = (0, 0, 100)
		//	fsample[2] = 0.5->csample[2] = (0, 170, 170)
		//	fsample[3] = 0.75->csample[3] = (170, 170, 0)
		//	fsample[4] = 1.0->csample[4] = (255, 255, 255)
	};





}

#endif