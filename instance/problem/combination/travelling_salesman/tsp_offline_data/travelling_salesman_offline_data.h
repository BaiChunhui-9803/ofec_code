/********* Begin Register Information **********
{
	"name": "TSP-offlineData",
	"identifier": "TravellingSalesmanOfflineData",
	"tags": [  "travelling salesman problem", "single-objective" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li & Yiya Diao
* Email: changhe.lw@google.com  Or diaoyiyacug@gmail.com
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/
// Created: 7 May 2024

#ifndef OFEC_TSP_OFFLINE_DATA_H
#define OFEC_TSP_OFFLINE_DATA_H

#include "../travelling_salesman.h"
#include "../../basin_divisioin/init_pop_basin.h"
#include "../../../../../utility/hnsw/hnsw_nbn/hnsw_model.h"




namespace ofec {
#define CAST_TSP_OD(pro) dynamic_cast<TravellingSalesmanOfflineData*>(pro)
	class Environment;

	class TravellingSalesmanOfflineData : public InitNBN_BasinBase, public TravellingSalesman {
		OFEC_CONCRETE_INSTANCE(TravellingSalesmanOfflineData)
	public:
		enum class ConstructionType {
			file, sampling, empty
		};
	protected:

		ConstructionType m_constructionType = ConstructionType::file;
		std::string m_file_dir;
		std::function<void(ofec::SolutionBase& sol, ofec::Environment* env)> m_eval_fun;

	public:
		void addInputParameters();


	protected:
		void initialize_(Environment* env) override;

	public:
		static void calculateNBNonline(
			const std::string& filedir,
			const std::string& filename,
			Hash& hash,
			nbn::HnswModel& model,
			std::vector<std::shared_ptr<SolutionBase>>& solbases,
			const std::function<void(ofec::SolutionBase& sol, ofec::Environment* env)>& eval_fun,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			ofec::Environment* env
		);
		static void insertNBNfromFile(
			const std::string& filedir,
			const std::string& filename,
			Hash& hash,
			nbn::HnswModel& model,
			std::vector<std::shared_ptr<SolutionBase>>& solbases,
			const std::function<void(ofec::SolutionBase& sol, ofec::Environment* env)>& eval_fun,
			std::vector<int>& belong,
			std::vector<double>& dis2parent,
			ofec::Environment* env
		);

	

		static void insertSolutions(
			const std::string& filedir,
			const std::string& filename,
			std::vector<std::shared_ptr<SolutionBase>>& solbases,
			const std::function<void(ofec::SolutionBase& sol, ofec::Environment* env)>& eval_fun,
			ofec::Environment* env
			);
		static void insertAlgorithmsPath(
			const std::string& filedir,
			const std::string& filename,
			std::vector<std::vector<std::vector<int>>>& solIds,
			ofec::Environment* env);
	};
}

#endif