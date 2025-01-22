#ifndef _BIT_H
#define _BIT_H

/* 
 * This header specifies the interface for the use of Binay Indexed Trees. 
 */

#include "LKH.h"
namespace LKH {
#ifdef __cplusplus
	extern "C" {
#endif

		void BIT_Make(int Size, LKH::LKHAlg *Alg);
		void BIT_Update(LKH::LKHAlg *Alg);

		int BIT_LoadDiff3Opt(LKH::LKHAlg::Node *t1, LKH::LKHAlg::Node *t2, LKH::LKHAlg::Node *t3, LKH::LKHAlg::Node *t4,
			LKH::LKHAlg::Node *t5, LKH::LKHAlg::Node *t6, LKH::LKHAlg *Alg);
		int BIT_LoadDiff4Opt(LKH::LKHAlg::Node *t1, LKH::LKHAlg::Node *t2, LKH::LKHAlg::Node *t3, LKH::LKHAlg::Node *t4,
			LKH::LKHAlg::Node *t5, LKH::LKHAlg::Node *t6, LKH::LKHAlg::Node *t7, LKH::LKHAlg::Node *t8, LKH::LKHAlg *Alg);
		int BIT_LoadDiff5Opt(LKH::LKHAlg::Node *t1, LKH::LKHAlg::Node *t2, LKH::LKHAlg::Node *t3, LKH::LKHAlg::Node *t4,
			LKH::LKHAlg::Node *t5, LKH::LKHAlg::Node *t6, LKH::LKHAlg::Node *t7, LKH::LKHAlg::Node *t8,
			LKH::LKHAlg::Node *t9, LKH::LKHAlg::Node *t10, int Case10, LKH::LKHAlg *Alg);
		int BIT_LoadDiff6Opt(LKH::LKHAlg::Node *t1, LKH::LKHAlg::Node *t2, LKH::LKHAlg::Node *t3, LKH::LKHAlg::Node *t4,
			LKH::LKHAlg::Node *t5, LKH::LKHAlg::Node *t6, LKH::LKHAlg::Node *t7, LKH::LKHAlg::Node *t8,
			LKH::LKHAlg::Node *t9, LKH::LKHAlg::Node *t10, LKH::LKHAlg::Node *t11, LKH::LKHAlg::Node *t12, LKH::LKHAlg *Alg);
#ifdef __cplusplus
	}
#endif
}
#endif
