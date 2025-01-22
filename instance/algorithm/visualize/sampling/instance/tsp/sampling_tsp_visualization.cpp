#include "sampling_tsp_visualization.h"
#include "../../../../../../core/global.h"
#include "../../../../../../core/parameter/variants_to_stream.h"
#include "../../../../../problem/combination/travelling_salesman/travelling_salesman.h"
#include "../../../../../../core/environment/environment.h"
#include "../../../../../../core/exception.h"


#ifdef OFEC_PLAYBACK
#include <buffer/datum_inclusion.h>
#endif // OFEC_PLAYBACK

#ifdef OFEC_STATISTICS
#include <record/datum_inclusion.h>
#endif // OFEC_STATISTICS


namespace ofec {
    void SamplingTSP_Visualization::addInputParameters() {
		m_input_parameters.add("directory to nbn offline data", new ofec::DirectoryName(m_file_dir,
			g_working_directory));
        m_input_parameters.add("numRun", new RangedInt(m_numRun,0,50, 0));

        
	}
	void SamplingTSP_Visualization::initialize_(Environment* env) {
		
        Algorithm::initialize_(env);
        auto tsp_pro = CAST_TSP(env->problem());
        auto tspfilename = tsp_pro->fileName();
        {

            ofec::ParameterVariantStream paramsStream;

            std::stringstream buf;
            std::ifstream in(m_file_dir + "/solIds_" + tspfilename + ".txt");
            if (!in) {
             //   std::cout << "not open file" << std::endl;
                std::cout << m_file_dir + "/solIds_" + tspfilename + ".txt" << std::endl;
                throw Exception("SamplingTSP_Visualization:: initialize: could not open file\t" + m_file_dir + "/solIds_" + tspfilename + ".txt");
            }
            buf << in.rdbuf();
            in.close();
            variants_stream::stringstream2parameterStream(buf, paramsStream);

            size_t dataSize;
            paramsStream >> dataSize;
            m_solIds.resize(dataSize);
            for (auto& it : m_solIds) {
                paramsStream >> dataSize;
                it.resize(dataSize);
                for (auto& it2 : it) {
                    paramsStream >> it2;
                }
            }

        }

        if (m_numRun >= m_solIds.size()) {
            throw Exception("Error at SamplingTSP_Visualization: numRun out of range");
        }

        m_numIter = 0;



#ifdef  OFEC_DATUM_OFFLINE_NBN_POP_DATA_H
        ofec::datum::g_curPopIds.working=true;
        auto& curPopIds(ofec::datum::g_curPopIds.m_totalIds);
        curPopIds.clear();
        for (auto& it : m_solIds[m_numRun]) {
            for (auto& it2 : it) {
                curPopIds.push_back(it2);
            }
        }
#endif
	}
	void SamplingTSP_Visualization::run_(Environment* env) {		


        while (!terminating()&&m_numIter< m_solIds[m_numRun].size()) {
#ifdef  OFEC_DATUM_OFFLINE_NBN_POP_DATA_H
            ofec::datum::g_curPopIds.m_popIds = m_solIds[m_numRun][m_numIter];
#endif
            
            datumUpdated(env);
            ++m_numIter;
        }
       
	}
}