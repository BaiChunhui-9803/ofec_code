/********* Begin Register Information **********
{
	"name": "OneMaxWModel",
	"identifier": "OneMaxWModel",
	"problem tags": [ "onemax", "single-objective" ],
	"dependency on libraries": [ "fmt", "Mynlohmann" ]
}
*********** End Register Information **********/

/*************************************************************************
* Project: Library of Open Frameworks for Evolutionary Computation (OFEC)
*************************************************************************
* Author: Changhe Li 
* Email: changhe.lw@google.com  
* Language: C++
*************************************************************************
*  This file is part of OFEC. This library is free software;
*  you can redistribute it and/or modify it under the terms of the
*  GNU General Public License as published by the Free Software
*  Foundation; either version 2, or (at your option) any later version.
*************************************************************************/



#ifndef WMODEL_ONE_MAX_H
#define WMODEL_ONE_MAX_H


#include "../one_max/one_max.h"
#include "../../../../utility/ioh_experimenter/include/ioh/problem/wmodel.hpp"


//https://iohprofiler.github.io/IOHexperimenter/cpp/ioh_problem_wmodel.html
//https://github.com/thomasWeise/BBDOB_W_Model


		/**
		 * @brief Construct a new WModel object
		 *
		 * @param problem_id the problem id
		 * @param name the name of the problem
		 * @param instance instance id
		 * @param n_variables the dimension of the problem
		 * @param dummy_select_rate select rate
		 * @param epistasis_block_size block size
		 * @param neutrality_mu neutrality parameter
		 * @param ruggedness_gamma ruggedness parameter
		 */

namespace ofec {
#define GET_OneMaxWModel(pro) dynamic_cast<OneMaxWModel*>(pro)

	class OneMaxWModel : public OneMax {
		OFEC_CONCRETE_INSTANCE(OneMaxWModel)
	public:
		enum class WModelProblemType { kWModelOneMax = 0, kWModelLeadingOnes };


	protected:


		const std::vector<std::string> m_wmodel_names = {"WModelOneMax","WModelLeadingOnes"};
		WModelProblemType m_type = WModelProblemType::kWModelOneMax;
		int m_instance_id = 0;
		double m_dummy_select_rate = 0.0;
		int m_epistasis_block_size = 0;
		int m_neutrality_mu = 0;
		int m_ruggedness_gamma = 0;
		std::unique_ptr<ioh::problem::WModel> m_model;

		enum class ConstructType {InputFile, Parameters};
		ConstructType m_constructType = ConstructType::Parameters;

		
	public:

		virtual void evaluate(const VariableBase& vars, std::vector<Real>& objs, std::vector<Real>& cons)override;


		std::string getProblemTotalName() {
			std::string proFileName = m_wmodel_names[static_cast<int>(m_type)] 
				+ "_Dim_" + std::to_string(m_number_variables)
				+"_InstanceId_" + std::to_string(m_instance_id)
				+"_DummySelectRate_" + std::to_string(int(m_dummy_select_rate*1000))
				+ "_EpistasisBlockSize_" + std::to_string(m_epistasis_block_size)
				+"_NeutralityMu_" +  std::to_string(m_neutrality_mu)
				+"_RuggednessGamma_" + std::to_string(m_ruggedness_gamma)
				;
			return std::move(proFileName);
		}


		void outputParams(std::ostream& out) {
			out << "type" << "\t" << static_cast<int>(m_type) << std::endl;
			out << "_Dim_" << "\t" << std::to_string(m_number_variables) << std::endl;
			out << "_InstanceId_" << "\t" << m_instance_id << std::endl;
				out << "_DummySelectRate_" << "\t" << m_dummy_select_rate  << std::endl;
				out << "_EpistasisBlockSize_" << "\t" << m_epistasis_block_size << std::endl;
				out << "_NeutralityMu_" << "\t" << m_neutrality_mu << std::endl;
				out << "_RuggednessGamma_" << "\t" << m_ruggedness_gamma << std::endl;
		}

		void inputParams(std::ifstream& in) {
			int ctype = 0;
			std::string info;
			in >> info >> ctype;
			m_type = static_cast<WModelProblemType>(ctype);
			in >> info >> m_number_variables;
			in >> info >> m_instance_id;
			in >> info >> m_dummy_select_rate;
			in >> info >> m_epistasis_block_size;
			in >> info >> m_neutrality_mu;
			in >> info >> m_ruggedness_gamma;
		/*	out << "_Dim_" << "\t" << std::to_string(m_num_vars) << std::endl;
			out << "_InstanceId_" << "\t" << m_instance_id << std::endl;
			out << "_DummySelectRate_" << "\t" << m_dummy_select_rate << std::endl;
			out << "_EpistasisBlockSize_" << "\t" << m_epistasis_block_size << std::endl;
			out << "_NeutralityMu_" << "\t" << m_neutrality_mu << std::endl;
			out << "_RuggednessGamma_" << "\t" << m_ruggedness_gamma << std::endl;*/
		}

		static void inputParams(std::istream& in, ofec::ParameterMap& param) {
			int m_type;// = WModelProblemType::kWModelOneMax;
			int m_instance_id = 0;
			double m_dummy_select_rate = 0.0;
			int m_epistasis_block_size = 0;
			int m_neutrality_mu = 0;
			int m_ruggedness_gamma = 0;
			int m_num_vars;
			//std::unique_ptr<ioh::problem::WModel> m_model;

			std::string info;
			in >> info >> m_type;
			in >> info >> m_num_vars;
			in >> info >> m_instance_id;
			in >> info >> m_dummy_select_rate;
			in >> info >> m_epistasis_block_size;
			in >> info >> m_neutrality_mu;
			in >> info >> m_ruggedness_gamma;

			param["type"] = m_type;
			param["number of variables"] = m_num_vars;
			param["_InstanceId_"] = m_instance_id;
			param["_DummySelectRate_"] = m_dummy_select_rate;
			param["_EpistasisBlockSize_"] = m_epistasis_block_size;
			param["_NeutralityMu_"] = m_neutrality_mu;
			param["_RuggednessGamma_"] = m_ruggedness_gamma;
		}


		
		

		void setInstanceId(int id) {
			m_instance_id = id;
		}
		//double m_dummy_select_rate = 0.0;
		void setDummySelectRate(double rate) {
			m_dummy_select_rate = rate;
		}
		//int epistasis_block_size = 0;
		void setEpistasisBlockSize(int blockSize) {
			m_epistasis_block_size = blockSize;
		}
		//int neutrality_mu = 0;
		void setNeutralityMu(int mu) {
			m_neutrality_mu = mu;
		}
		//int ruggedness_gamma = 0;
		/*
		* The integer parameter gamma ranges from 0 (r_0 is the canonical permutation and leaves the objective values unchanged) to n(n-1)/2,
		*/
		void setRuggednessGamma(int gamma) {
			m_ruggedness_gamma = gamma;
		}

	protected:
		void addInputParameters();
		void initialize_(Environment* env)override;
		virtual void updateOptima(Environment* env)override;
	};

}

#endif // !ONE_MAX_H

