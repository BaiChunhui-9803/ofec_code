#include "travelling_salesman_offline_data.h"
#include "../../../../algorithm/visualize/sampling/sampling_multiThread.h"
#include "../../../../../utility/hnsw/hnsw_nbn/hnsw_nbn.h"
#include "../../../../../utility/function/custom_function.h"
#include "../../../../../core/global.h"
#include "../../../../../core/environment/environment.h"



void ofec::TravellingSalesmanOfflineData::addInputParameters()
{
	m_input_parameters.add("directory to nbn offline data", new ofec::DirectoryName(m_file_dir, g_working_directory));
	m_input_parameters.add("construction way", new ofec::Enumeration(m_constructionType,
		{
			"read from file",
			"sample online",
			"nothing"
		},
		ConstructionType::file));

}

void ofec::TravellingSalesmanOfflineData::initialize_(Environment* env)
{
	TravellingSalesman::initialize_(env);
	auto rnd = std::make_shared<ofec::Random>(0.5);
	m_hnswModel.initialize(env, rnd.get(), std::thread::hardware_concurrency());
	
	m_eval_fun =
		[&](ofec::SolutionBase& sol, ofec::Environment* env) {
		//using namespace ofec;
		evaluate(sol.variableBase(), sol.objective(), sol.constraint());
		//sol.evaluate(env, false);
		ofec::Real pos = env->problem()->optimizeMode(0) == ofec::OptimizeMode::kMaximize ? 1 : -1;
		sol.setFitness(pos * sol.objective(0));
		};

	std::vector<std::shared_ptr<SolutionBase>> sols;

	switch (m_constructionType)
	{
	case ofec::TravellingSalesmanOfflineData::ConstructionType::file:
		insertNBNfromFile(m_file_dir, m_file_name, m_hash,m_hnswModel, sols, m_eval_fun, m_belong, m_dis2parent, env);
		break;
	case ofec::TravellingSalesmanOfflineData::ConstructionType::sampling:
		calculateNBNonline(m_file_dir, m_file_name, m_hash, m_hnswModel, sols, m_eval_fun, m_belong, m_dis2parent, env);
		break;
	default:
		break;
	}

	setOneBasin();
}

void ofec::TravellingSalesmanOfflineData::calculateNBNonline(
	const std::string& filedir,
	const std::string& filename,
	Hash& hashStruct,
	nbn::HnswModel& model,
	std::vector<std::shared_ptr<SolutionBase>>& solbases,
	const std::function<void(ofec::SolutionBase& sol, ofec::Environment* env)>& eval_fun,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	ofec::Environment* env
) {


	//auto tsp_pro = CAST_TSP(env->problem());

	ofec::ParameterMap param;
	std::string algName = "EAX-TSP-sampling";
	param["dataFile1"] = filename;

	double env_seed = env->random() == nullptr ? 0.5 : env->random()->seed;
	double pro_seed = env->problem()->random() == nullptr ? 0.5 : env->problem()->random()->seed;


	solbases.clear();
	std::vector<std::shared_ptr<ofec::SolutionInfo>> solInfos;
	SamplingMutithread::runAlgMultiTask(solbases, solInfos, env_seed, pro_seed, env->problem()->name(), param, algName, param, 30);

	//auto sols = SamplingMutithread::runAlgMultiTask(
	//	env_seed, pro_seed,
	//	env->problem()->name(), param, algName, param, 30);

	{
		int totalSols = 0;
		std::set<unsigned long long> setSols;

		for (auto& it : solbases) {

			++totalSols;
			auto hash = hashStruct.calHash(*it);
			if (setSols.find(hash) == setSols.end()) {
				setSols.insert(hash);

			}
		}

		//for (auto& it : sols) {
		//	for (auto& it2 : it) {
		//		for (auto& it3 : it2) {

		//			++totalSols;
		//			auto hash = hashStruct.calHash(*it3);
		//			if (setSols.find(hash) == setSols.end()) {
		//				solbases.push_back(it3->getSharedPtr());
		//				setSols.insert(hash);

		//			}
		//		}
		//	}
		//}
	}

	std::vector<ofec::SolutionBase*> ptr_solbases;
	for (auto& it : solbases) {
		ptr_solbases.push_back(it.get());
	}
	//	ofec::TravellingSalesman* pro_t = dynamic_cast<ofec::TravellingSalesman*>(pro.get());
	//std::cout << "number of solutions\t" << solbases.size();

	UTILITY::evaluateRandomSolutionsMultiThreads(ptr_solbases, env, eval_fun);


	{
		std::vector<int> solIds;
		//	new_datum->m_hnswModel2.setNumberThreads(1);
		model.addDataMutliThread(ptr_solbases, solIds);

	}
	std::vector<double> vFitness;
	for (auto& it : solbases) {
		vFitness.push_back(it->fitness());
	}

	{
		nbn::HnswModelNBN model2;
		model2.copy(model);
		//model2.setNumberThreads(1);
		model2.updateFitness(vFitness);
		model2.calculateNBN(belong, dis2parent);
	}
}
void ofec::TravellingSalesmanOfflineData::insertNBNfromFile(
	const std::string& filedir,
	const std::string& filename,
	Hash& hashStruct,
	nbn::HnswModel& model,
	std::vector<std::shared_ptr<SolutionBase>>& solbases,
	const std::function<void(ofec::SolutionBase& sol, ofec::Environment* env)>& eval_fun,
	std::vector<int>& belong,
	std::vector<double>& dis2parent,
	ofec::Environment* env
) {
	//solbases.clear();

	//std::vector<std::shared_ptr<ofec::SolutionBase>> solbasesCopy;
	//{
	//	ofec::ParameterVariantStream paramsStream;
	//	std::stringstream buf;
	//	std::ifstream in(filedir + "/solVariables_" + filename + ".txt");
	//	if (!in) {
	//		std::cout << "not open file" << std::endl;
	//		std::cout << filedir + "/solVariables_" + filename + ".txt" << std::endl;
	//	}
	//	buf << in.rdbuf();
	//	in.close();
	//	ofec::variants_stream::stream2ParametersMutithreadLine(buf, paramsStream);


	//	size_t dataSize;
	//	paramsStream >> dataSize;
	//	solbasesCopy.resize(dataSize);
	//	for (auto& it : solbasesCopy) {
	//		it.reset(env->problem()->createSolution());
	//		auto& cursol = dynamic_cast<TravellingSalesman::SolutionType&>(*it);
	//		paramsStream >> cursol.variable().vect();
	//	}
	//}
	//for (auto& it : solbasesCopy) {
	//	solbases.push_back(it);
	//}

	//UTILITY::evaluateRandomSolutionsMultiThreads(ptr_solbases, env, eval_fun);



	insertSolutions(filedir, filename, solbases, eval_fun, env);


	std::vector<ofec::SolutionBase*> ptr_solbases;
	for (auto& it : solbases) {
		ptr_solbases.push_back(it.get());
	}


	model.setSolutions(ptr_solbases);

	//auto filedir = dynamic_cast<FileName*>(m_input_parameters.at("directory to data"));
	//auto m_file_dir = filedir->directoryPath();
	{
		ofec::ParameterVariantStream paramsStream;
		std::stringstream buf;
		std::ifstream in(filedir + "/hnsw_parameters_" + filename + ".txt");
		if (!in) {
			std::cout << "not open file" << std::endl;
			std::cout << filedir + "/hnsw_parameters_" + filename + ".txt" << std::endl;
			throw ofec::Exception("could not find hnsw model file\n");
		}
		buf << in.rdbuf();
		in.close();

		std::cout << "buf size\t" << buf.str().size() << std::endl;
		ofec::variants_stream::stringstream2parameterStream(buf, paramsStream);

		model.configsfromParameters(paramsStream);

	}


	{

		//	ofec::variants_stream::m_out.clear();


		ofec::ParameterVariantStream paramsStream;
		std::stringstream buf;
		std::ifstream in(filedir + "/hnsw_data_" + filename + ".txt");
		if (!in) {
			throw ofec::Exception("could not find hnsw model file\n");
		}
		buf << in.rdbuf();
		in.close();

		//   testBufReader(buf);
		ofec::variants_stream::stringstream2parameterStream(buf, paramsStream);
		//    compare(paramsStream_in, paramsStream_out);
		model.datasFromParameters(paramsStream);

	}


	std::vector<double> vFitness;
	for (auto& it : solbases) {
		vFitness.push_back(it->fitness());
	}

	{
		nbn::HnswModelNBN model2;
		model2.copy(model);
		//model2.setNumberThreads(1);
		model2.updateFitness(vFitness);
		model2.calculateNBN(belong, dis2parent);
	}

	std::cout << "finish calculating nbn" << std::endl;
}

void ofec::TravellingSalesmanOfflineData::insertSolutions(
	const std::string& filedir, 
	const std::string& filename, 
	std::vector<std::shared_ptr<SolutionBase>>& solbases, 
	const std::function<void(ofec::SolutionBase& sol, ofec::Environment* env)>& eval_fun,
	ofec::Environment* env
	)
{
	solbases.clear();
	std::vector<std::shared_ptr<ofec::SolutionBase>> solbasesCopy;

	{
		ofec::ParameterVariantStream paramsStream;
		std::stringstream buf;
		std::ifstream in(filedir + "/solVariables_" + filename + ".txt");
		if (!in) {
			std::cout << "not open file" << std::endl;
			std::cout << filedir + "/solVariables_" + filename + ".txt" << std::endl;
		}
		buf << in.rdbuf();
		in.close();
		ofec::variants_stream::stringstream2parameterStream(buf, paramsStream);


		size_t dataSize;
		paramsStream >> dataSize;
		solbasesCopy.resize(dataSize);
		for (auto& it : solbasesCopy) {
			it.reset(env->problem()->createSolution());
			auto& cursol = dynamic_cast<TravellingSalesman::SolutionType&>(*it);
			paramsStream >> cursol.variable().vect();
		}
	}
	for (auto& it : solbasesCopy) {
		solbases.push_back(it);
	}
	std::vector<ofec::SolutionBase*> ptr_solbases;
	for (auto& it : solbases) {
		ptr_solbases.push_back(it.get());
	}

	UTILITY::evaluateRandomSolutionsMultiThreads(ptr_solbases, env, eval_fun);
}

void ofec::TravellingSalesmanOfflineData::insertAlgorithmsPath(const std::string& filedir, 
	const std::string& filename, 
	std::vector<std::vector<std::vector<int>>>& solIds, 
	ofec::Environment* env)
{

	ofec::ParameterVariantStream paramsStream;

	std::stringstream buf;
	std::ifstream in(filedir + "/solIds_" + filename + ".txt");
	if (!in) {
		//   std::cout << "not open file" << std::endl;
		std::cout << filedir + "/solIds_" + filename + ".txt" << std::endl;
		throw Exception("SamplingTSP_Visualization:: initialize: could not open file\t" + filedir + "/solIds_" + filename + ".txt");
	}
	buf << in.rdbuf();
	in.close();
	variants_stream::stringstream2parameterStream(buf, paramsStream);

	size_t dataSize;
	paramsStream >> dataSize;
	solIds.resize(dataSize);
	for (auto& it : solIds) {
		paramsStream >> dataSize;
		it.resize(dataSize);
		for (auto& it2 : it) {
			paramsStream >> it2;
		}
	}
}

