#include "matlab_nbn_visualization.h"
#ifdef  OFEC_DEMO_MATLAB_PLOT

#ifdef  OFEC_MSVC
//#include "NBN_figure_out.h"
#include "nbnFigureOut.h"
#endif //  OFEC_MSVC


//#include "figures_ExportFigure.h"
//#include "Graph3DColorMarkerAxis_ExportFigure.h"
#ifdef  OFEC_UNIX
#include <libNBNfigureOut.h>
#endif //  OFEC_UNIX
#endif //  



#include "../myexcept.h"

#include <cmath>


void ofec::MatlabNBNVisualization::plot(int figureId)
const {
	DrawMatlabFiguresOutput(figureId,
		m_opt3D_pts_colors, m_opt3D_pts, m_opt3D_edges,
		m_marker_idxs, m_markerSize, m_cut, contour3Level,
		m_filepath, m_filetype, m_fontSize, m_figure_type
	);
}


void ofec::MatlabNBNVisualization::DrawMatlabFiguresOutput(int figureId,
	const std::vector<std::array<double, 3>>& opt3D_pts_colors,
	const std::vector<std::array<double, 3>>& opt3D_pts,
	const std::vector<std::pair<int, int>>& opt3D_edges,
	const std::vector<int>& marker_idxs, int markerSize, int m_cut,
	int contour3Level,
	const std::string& filepath, const std::string& filetype, int fontSize,
	int type)
{
	std::string new_filepath = filepath;

//.	new_filepath += std::string("-sample_size-") + std::to_string(opt3D_pts.size()) + std::string("-");


	using namespace ofec;
#ifdef OFEC_DEMO_MATLAB_PLOT

	std::cout << "OFEC_DEMO_MATLAB_PLOT ok" << std::endl;

	mxDouble nInData = figureId;
	mwArray nIn(nInData);


	int ptrSize(opt3D_pts.size());

	mwArray XIn(ptrSize, 1, mxDOUBLE_CLASS);
	mwArray YIn(ptrSize, 1, mxDOUBLE_CLASS);
	mwArray ZIn(ptrSize, 1, mxDOUBLE_CLASS);


	int numEdge(opt3D_edges.size());
	mwArray EdgeBIn(1, numEdge, mxUINT32_CLASS);
	mwArray EdgeAIn(1, numEdge, mxUINT32_CLASS);


	{
		std::vector<mxUint32> dataInt;
		dataInt.resize(numEdge);
		int idx(0);
		if (numEdge > 0) {
			for (auto& it : opt3D_edges) {
				dataInt[idx++] = it.first + 1;
			}
			EdgeBIn.SetData(&dataInt[0], dataInt.size());
			idx = 0;
			for (auto& it : opt3D_edges) {
				dataInt[idx++] = it.second + 1;
			}
			EdgeAIn.SetData(&dataInt[0], dataInt.size());

		}

	}
	int colorSize(opt3D_pts_colors.size());

	mwArray colorRIn(colorSize, 1, mxDOUBLE_CLASS);
	mwArray colorGIn(colorSize, 1, mxDOUBLE_CLASS);
	mwArray colorBIn(colorSize, 1, mxDOUBLE_CLASS);


	{
		std::vector<mxDouble> dataDouble;
		dataDouble.resize(ptrSize);
		int idx(0);
		for (auto& it : opt3D_pts) {
			dataDouble[idx++] = it[0];
		}
		XIn.SetData(&dataDouble[0], dataDouble.size());


		idx = 0;
		for (auto& it : opt3D_pts) {
			dataDouble[idx++] = it[1];
		}
		YIn.SetData(&dataDouble[0], dataDouble.size());


		idx = 0;
		for (auto& it : opt3D_pts) {
			dataDouble[idx++] = it[2];
		}
		ZIn.SetData(&dataDouble[0], dataDouble.size());



		dataDouble.resize(colorSize);

		idx = 0;
		for (auto& it : opt3D_pts_colors) {
			dataDouble[idx++] = it[0];
		}
		colorRIn.SetData(&dataDouble[0], dataDouble.size());



		idx = 0;
		for (auto& it : opt3D_pts_colors) {
			dataDouble[idx++] = it[1];
		}
		colorGIn.SetData(&dataDouble[0], dataDouble.size());


		idx = 0;
		for (auto& it : opt3D_pts_colors) {
			dataDouble[idx++] = it[2];
		}
		colorBIn.SetData(&dataDouble[0], dataDouble.size());
	}

	mwArray markerPtrsIn(1, marker_idxs.size(), mxUINT32_CLASS);
	if (!marker_idxs.empty()) {
		std::vector<mxUint32> dataInt;
		dataInt.resize(marker_idxs.size());
		for (int idx(0); idx < dataInt.size(); ++idx) {
			dataInt[idx] = marker_idxs[idx] + 1;
		}
		markerPtrsIn.SetData(&dataInt[0], dataInt.size());
	}
	nInData = markerSize;
	mwArray PsizeIn(nInData);

	mxDouble cutInData = m_cut;
	mwArray cutIn(cutInData);



	mxDouble contour3LevelInData = contour3Level;
	mwArray contour3LevelIn(contour3LevelInData);




	const char* filetypeInData = filetype.c_str();
	mwArray filetypeIn(filetypeInData);


	mxDouble fontSizeInData = fontSize;
	mwArray fontSizeIn(fontSizeInData);

	try
	{
		if (type & kContour2D) {

			std::string cur_filepath = new_filepath + "_concour2D";
			const char* filenameInData = cur_filepath.c_str();
			mwArray filenameIn(filenameInData);
			concour2D_ExportFigure(nIn, XIn, YIn, ZIn, cutIn, contour3LevelIn, filenameIn, filetypeIn, fontSizeIn);
		}
		if (type & kContour2DNBN) {


			std::vector<double> dataDouble(ptrSize, 0);
			mwArray Z0In(ptrSize, 1, mxDOUBLE_CLASS);
			Z0In.SetData(&dataDouble[0], dataDouble.size());

			std::string cur_filepath = new_filepath + "_concour2DNBN";
			const char* filenameInData = cur_filepath.c_str();
			mwArray filenameIn(filenameInData);
			concour2DNBN_ExportFigure(nIn, EdgeBIn, EdgeAIn, XIn, YIn, Z0In, colorRIn, colorGIn, colorBIn, cutIn, markerPtrsIn, PsizeIn, contour3LevelIn, filenameIn, filetypeIn, fontSizeIn);


		}
		if (type & kContour3DNBN) {


			std::string cur_filepath = new_filepath + "_Contour3DNBN";
			const char* filenameInData = cur_filepath.c_str();
			mwArray filenameIn(filenameInData);
			concour3DNBN_ExportFigure(nIn, EdgeBIn, EdgeAIn, XIn, YIn, ZIn, colorRIn, colorGIn, colorBIn, cutIn, markerPtrsIn, PsizeIn, contour3LevelIn, filenameIn, filetypeIn, fontSizeIn);

		}
		//if (type & PlotFiguresOutput::kGraph3DNBN) {
		//	std::string cur_filepath = new_filepath + "_Graph3DNBN";
		////	std::string cur_filepath = "name_Graph3DNBN";
		//	const char* filenameInData = cur_filepath.c_str();
		//	mwArray filenameIn(filenameInData);
		//	Graph3DColorMarker_ExportFigure(nIn, EdgeBIn, EdgeAIn, XIn, YIn, ZIn, colorRIn, colorGIn, colorBIn, markerPtrsIn, PsizeIn, filenameIn, filetypeIn, fontSizeIn);
		//}


		if (type & kGraph3DNBNAxis) {
			std::string cur_filepath = new_filepath + "_Graph3DNBNAxis";

			//	std::string cur_filepath = "name_Graph3DNBN";
			const char* filenameInData = cur_filepath.c_str();
			mwArray filenameIn(filenameInData);


			mxDouble colorbar_labelInData = 0;
			mwArray colorbar_labelIn(colorbar_labelInData);
			mxDouble axis_labelInData = 1;
			mwArray axis_labelIn(axis_labelInData);

			Graph3DColorMarkerAxis_ExportFigure(nIn, EdgeBIn, EdgeAIn, XIn, YIn, ZIn, colorRIn, colorGIn, colorBIn, markerPtrsIn, PsizeIn, filenameIn, filetypeIn, fontSizeIn, colorbar_labelIn, axis_labelIn);

			//	Graph3DColorMarker_ExportFigure(nIn, EdgeBIn, EdgeAIn, XIn, YIn, ZIn, colorRIn, colorGIn, colorBIn, markerPtrsIn, PsizeIn, filenameIn, filetypeIn, fontSizeIn);

		}

		//concour2DNBN_ExportFigure(nIn, EdgeBIn, EdgeAIn, xiIn, yiIn, ziIn, colorRIn, colorGIn, colorBIn, cutIn, markerPtrsIn, PsizeIn, contour3LevelIn, filenameIn, filetypeIn, fontSizeIn);

	//	concour3D(nIn, XIn, YIn, ZIn, cutIn, contour3LevelIn);

		//concour3DNBN(nIn, EdgeBIn, EdgeAIn, XIn, YIn, ZIn, colorRIn, colorGIn, colorBIn, cutIn, markerPtrsIn, PsizeIn, contour3LevelIn);
	}
	catch (const mwException& e)
	{
		std::string error_info = "concour3D.cpp: PlotMatlabGraphColor3Marker\t" + std::string(e.what());
	//	THROW(error_info.c_str());
		//	THROW(error_info.c_str());
	}
	catch (...)
	{
		std::string error_info = "concour3D.cpp: PlotMatlabGraphColor3Marker\t" + std::string("Unexpected error thrown");
	//	THROW(error_info.c_str());
	}



#endif // OFEC_DEMO_MATLAB_PLOT

}

void ofec::colorMapper::init()
{
	for (auto& it : m_dcolor_table) {
		for (auto& it2 : it) {
			it2 /= 255.;
		}
	}


}

void ofec::colorMapper::getColor(double val, std::array<double, 3>& dcolor)
{
	int from = std::floor(val * 4);
	if (from >= 4) {
		dcolor = m_dcolor_table[from];
	}
	else {
		int to = from + 1;
		for (int ridx(0); ridx < 3; ++ridx) {
			dcolor[ridx] = m_dcolor_table[from][ridx] * (m_value_div[to] - val) + m_dcolor_table[to][ridx] * (val - m_value_div[from]);
		}
	}
}
