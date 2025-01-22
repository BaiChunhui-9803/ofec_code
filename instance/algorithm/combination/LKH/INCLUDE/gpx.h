#ifndef _GPX_H
#define _GPX_H

//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
#include "LKH.h"
namespace LKH {




#ifdef __cplusplus
	extern "C" {
#endif
#define new_int(n) ((int *) calloc(n, sizeof(int)))
#define new_tour(n) ((tour *) calloc(n, sizeof(tour)))
#define new_gate_structure(n) ((gate_structure *)\
                               calloc(n, sizeof(gate_structure)))

		typedef struct Adj
		{
			int vertex;
			struct Adj *nextAdj;
		} Adj;

		typedef struct Graph
		{
			int numVertices;
			Adj **firstAdj, **lastAdj;
		} Graph;

		Graph *new_Graph(int n);
		void insertEdge(Graph *g, int v1, int v2);
		void freeGraph(Graph *g);
		void compCon(Graph *g, int *vector_comp);

		typedef struct {
			int num;
			int time;
		} gate_structure;

		typedef struct {
			gate_structure* inputs = nullptr;
			gate_structure *outputs = nullptr;
			gate_structure first_entry ;
			gate_structure last_exit ;
			GainType fitness = 0;
		} tour;



		struct gpxTour {
			int* id = nullptr;;                 // vector with the id(candidate component) of each node
			int* size_ = nullptr;;               // vector
			int* n_inputs = nullptr;;
			int* n_outputs = nullptr;;
			tour* blue = nullptr, * red = nullptr;;
			int n_=0;                   // size of the tours 
			int** M_neigh = nullptr;;           // neighbourhood matrix: first collumn indicates the
			// number of neighbours, and the collumns 2 and 3
			// indicate the index of the neighbours
			int** M_neigh2 = nullptr;;          // neighbourhood matrix: the collumns indicate the
			// number of i conections to the neighbours indicated
			// in collumns 2 and 3

			int* test = nullptr;;                      // test of the candidtes: 1 - true component; 0 - otherwise

		};





		GainType gpx(int *solution_blue, int *solution_red, int *offspring, LKH::LKHAlg *Alg);
		void new_candidates(int *vector_comp, int n_new, gpxTour& curTour, LKHAlg::gpxInfo& info);
		void free_candidates(gpxTour& curTour, LKHAlg::gpxInfo& info);
		void findInputs(int *sol_blue, int *sol_red, int cand, gpxTour& curTour);
		void testComp(int cand, gpxTour& curTour);
		int testUnfeasibleComp(int *sol_blue, int cand,gpxTour& curTour);
		void fusion(int *sol_blue, int *sol_red, int cand, gpxTour& curTour);
		void fusionB(int *sol_blue, int *sol_red, int cand, gpxTour& curTour);
		//void fusionB_v2(int *sol_blue, int *sol_red);
		GainType off_gen(int *sol_blue, int *sol_red, int *offspring,
			int *label_list, gpxTour& curTour, LKH::LKHAlg *Alg);

		int *alloc_vectori(int lines);
		int **alloc_matrixi(int lines, int collums);
		void dealloc_matrixi(int **Matrix, int lines);
		int weight(int i, int j, LKH::LKHAlg *Alg);
		int d4_vertices_id(int *solution_blue, int *solution_red,
			int *d4_vertices, int *common_edges_blue,
			int *common_edges_red, LKHAlg::gpxInfo& info);
		void insert_ghost(int *solution, int *solution_p2, int *d4_vertices,
			int *label_list_inv, LKHAlg::gpxInfo& info);
		void tourTable(int *solution_blue_p2, int *solution_red_p2,
			int *solution_red, int *label_list, int *label_list_inv,
			int *vector_comp, int n_new, int *common_edges_blue_p2,
			int *common_edges_red_p2, LKHAlg::gpxInfo& info);
#ifdef __cplusplus
	}
#endif
}
#endif
