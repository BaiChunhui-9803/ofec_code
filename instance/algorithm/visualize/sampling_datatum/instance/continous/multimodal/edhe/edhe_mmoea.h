///********* Begin Register Information **********
//{
//	"name": "EDHE-MMOEA",
//	"identifier": "EDHE_MMOEA",
//	"tags": [ "continuous", "single-objective" ],
//	"dependency on libraries": [ "LibTorch" ]
//}
//*********** End Register Information **********/
//
//#ifndef OFEC_EDHE_MMOEA_H
//#define OFEC_EDHE_MMOEA_H
//
//#include "../../../../template/framework/edhe/continuous/edhec.h"
//
//namespace ofec {
//	class EDHE_MMOEA : public EDHEC {
//		OFEC_CONCRETE_INSTANCE(EDHE_MMOEA)
//	protected:
//		std::string m_mmoea_name;
//		size_t m_pop_size;
//		bool m_god_mode;
//
//		void addInputParameters();
//		void run_(Environment *env) override;
//
//		void initHills(Environment *env);
//		void selectHill(Environment *env);
//		void runDEnrand1(Environment *env);
//		void runANDE(Environment *env);
//	};
//}
//
//#endif // !OFEC_EDHE_MMOEA_H