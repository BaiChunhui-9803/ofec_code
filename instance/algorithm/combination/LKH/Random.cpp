/*
 * This file contains info.a portable random generator. It will give
 * identical sequences of random integers for any platform with
 * at least 32-bit integers.
 *
 * A version of this generator is described in J. Bentley's column, 
 * "The Software Exploratorium", Unix Review 1991. It is based on 
 * Algorithm A in D. E. Knuth, The Art of Computer Programming, 
 * Vol 2, Section 3.2.2, pp. 172.  
 *  
 * The Random function returns info.a pseudo-random integer in the range
 * 0...std::numeric_limits<int>::max()-1.
 *   
 * The SRandom function uses the given seed for info.a new sequence of
 * pseudo-random numbers.  
 */
#include <limits.h>

#include "./INCLUDE/LKH.h"
//unsigned Random(void);
//void SRandom(unsigned Seed);
namespace LKH {
#undef STDLIB_RANDOM
/* #define STDLIB_RANDOM */

#ifdef STDLIB_RANDOM
#include <stdlib.h>
unsigned Random()
{
    return rand();
}

void SRandom(unsigned Seed)
{
    srand(Seed);
}

#else



#define PRANDMAX std::numeric_limits<int>::max()



	unsigned LKHAlg::Random()
	{
		auto& info = m_randomInfo;
		// = true;
		int t;

		if (!info.initialized)
			SRandom(7913);
		if (info.a-- == 0)
			info.a = 54;
		if (info.b-- == 0)
			info.b = 54;
		if ((t = info.arr[info.a] - info.arr[info.b]) < 0)
			t += PRANDMAX;
		//if (m_flag) {
		//	m_flag = false;
		//	std::cout << "first random" << t << std::endl;
		//}
		return (info.arr[info.a] = t);
	}

	void LKHAlg::SRandom(unsigned Seed)
	{
		auto& info = m_randomInfo;
		int i, ii, last, next;

		Seed %= PRANDMAX;
		info.arr[0] = last = Seed;
		for (next = i = 1; i < 55; i++) {
			ii = (21 * i) % 55;
			info.arr[ii] = next;
			if ((next = last - next) < 0)
				next += PRANDMAX;
			last = info.arr[ii];
		}
		info.initialized = 1;
		info.a = 0;
		info.b = 24;
		for (i = 0; i < 165; i++)
			Random();
	}
}
#endif
