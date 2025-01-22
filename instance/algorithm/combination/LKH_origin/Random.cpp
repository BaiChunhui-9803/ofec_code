/*
 * This file contains a portable random generator. It will give
 * identical sequences of random integers for any platform with
 * at least 32-bit integers.
 *
 * A version of this generator is described in J. Bentley's column, 
 * "The Software Exploratorium", Unix Review 1991. It is based on 
 * Algorithm A in D. E. Knuth, The Art of Computer Programming, 
 * Vol 2, Section 3.2.2, pp. 172.  
 *  
 * The Random function returns a pseudo-random integer in the range
 * 0...std::numeric_limits<int>::max()-1.
 *   
 * The SRandom function uses the given seed for a new sequence of
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

	static thread_local int a = 0, b = 24, arr[55], initialized = 0;

	unsigned LKHAlg::Random()
	{
		// = true;
		int t;

		if (!initialized)
			SRandom(7913);
		if (a-- == 0)
			a = 54;
		if (b-- == 0)
			b = 54;
		if ((t = arr[a] - arr[b]) < 0)
			t += PRANDMAX;
		//if (m_flag) {
		//	m_flag = false;
		//	std::cout << "first random" << t << std::endl;
		//}
		return (arr[a] = t);
	}

	void LKHAlg::SRandom(unsigned Seed)
	{
		int i, ii, last, next;

		Seed %= PRANDMAX;
		arr[0] = last = Seed;
		for (next = i = 1; i < 55; i++) {
			ii = (21 * i) % 55;
			arr[ii] = next;
			if ((next = last - next) < 0)
				next += PRANDMAX;
			last = arr[ii];
		}
		initialized = 1;
		a = 0;
		b = 24;
		for (i = 0; i < 165; i++)
			Random();
	}
}
#endif
