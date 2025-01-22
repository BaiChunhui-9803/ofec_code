#include "wmodel.h"



void ofec::OneMaxWModel::addInputParameters() {

	m_input_parameters.add("WModelProblemType", new Enumeration(
		m_type, {"WModelOneMax","WModelLeadingOnes"}, WModelProblemType::kWModelOneMax));

	//m_input_parameters.add("Construction Type", new Enumeration(
	//	m_constructType, {"InputFile","Parameters"}, ConstructType::Parameters));
	//m_input_parameters.add("dataFile1", new FileName(
	//	m_constructType, { "InputFile","Parameters" }, ConstructType::Parameters));
}
void ofec::OneMaxWModel::evaluate(const VariableBase& vars, std::vector<Real>& objs, std::vector<Real>& cons)
{
	auto& x = dynamic_cast<const VariableType&>(vars);
	double evaluate_obj = (*m_model)(x.vect());
	objs.front() = evaluate_obj;
}

void ofec::OneMaxWModel::updateOptima(Environment* env) {
	m_optima.reset(new Optima<VariableType>());
	if (m_model->optimum().exists()) {


		if (!m_model->optimum().x.empty()) {
			SolutionType sol(m_number_objectives, m_number_constraints, m_number_variables);
			sol.variable().vect() = m_model->optimum().x;
			sol.objective().front()= m_model->optimum().y;
			dynamic_cast<Optima<VariableType>&>(*m_optima).appendSolution(sol);
		}
		else {
			std::vector<Real> objs(m_number_objectives);
			objs.front() = m_model->optimum().y;
			m_optima->appendObjective(objs);
		}
	}

}

void  ofec::OneMaxWModel::initialize_(Environment* env) {
	OneMax::initialize_(env);
	//if (v.has("dataFile1")) {
	//	std::ifstream in(v.get<std::string>("dataFile1"));
	//	inputParams(in);
	//	in.close();
	//}
	//else {
	//	m_num_vars = v.get<int>("number of variables");
	//	//m_type = static_cast<WModelProblemType>(v.get<int>("type"));

	//	//m_instance_id = v.get<int>("_InstanceId_");
	//	//m_dummy_select_rate = v.get<double>("_DummySelectRate_");
	//	//m_epistasis_block_size = v.get<double>("_EpistasisBlockSize_");
	//	//m_neutrality_mu = v.get<double>("_NeutralityMu_");
	//	//m_ruggedness_gamma = v.get<double>("_RuggednessGamma_");

	//}

	resizeObjective(1);
	resizeConstraint(0);
	m_optimize_mode[0] = OptimizeMode::kMaximize;


	//m_optima.reset(new Optima<VarVec<int>>());
	if (m_type == WModelProblemType::kWModelLeadingOnes) {
		m_model.reset(new ioh::problem::wmodel::WModelLeadingOnes(
			m_instance_id,m_number_variables, m_dummy_select_rate, m_epistasis_block_size, m_neutrality_mu, m_ruggedness_gamma));
	}
	else if (m_type == WModelProblemType::kWModelOneMax) {
		m_model.reset(new ioh::problem::wmodel::WModelOneMax(
			m_instance_id, m_number_variables, m_dummy_select_rate, m_epistasis_block_size, m_neutrality_mu, m_ruggedness_gamma));

	}



}
