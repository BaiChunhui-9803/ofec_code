#include "tsp_benchmark_generator.h"
#include <cmath>


void ofec::PortGen::sprand(int seed )
{	int i,ii;
	int last,next;
	arr[0] = last = seed;
	next = 1;
	for (i=1;i<55;i++) {
		ii = (21*i)%55;
		arr[ii] = next;
		next = last - next;
		if (next < 0) next += PRANDMAX;
		last = arr[ii];
	}
	a = 0;
	b = 24;
	for (i=0;i<165;i++) last = lprand();

}

int ofec::PortGen::lprand(void)
{
	long t;
	if (a-- == 0) a = 54;
	if (b-- == 0) b = 54;
	t = arr[a] - arr[b];
	if (t < 0) t += PRANDMAX;
	arr[a] = t;
	return t;
}

void ofec::PortGen::generate(int N, int seed, std::ostream& out)
{
	int factor = PRANDMAX / MAXCOORD;
	int i;
	int x, y;


	/* initialize random number generator */

	sprand(seed);

	out << "NAME : portgen-" << N << "-" << seed << std::endl;
	printf("NAME : portgen-%d-%d\n", N, seed);
	out << "COMMENT : portgen N=" << N << ", seed=" << seed << std::endl;
	printf("COMMENT : portgen N=%d, seed=%d\n", N, seed);
	out << "TYPE : TSP" << seed << std::endl;
	printf("TYPE : TSP\n");
	out << "DIMENSION : " << N << std::endl;
	printf("DIMENSION : %d\n", N);
	out << "EDGE_WEIGHT_TYPE : EUC_2D"<< std::endl;
	printf("EDGE_WEIGHT_TYPE : EUC_2D\n");
	out << "NODE_COORD_SECTION" << std::endl;
	printf("NODE_COORD_SECTION\n");

	for (i = 1; i <= N; i++) {
		x = lprand() / factor;
		y = lprand() / factor;
		out << i<<" "<<x<<" "<<y<< std::endl;
		printf("%d %d %d\n", i, x, y);
	}
	
}

void ofec::PortCGen::generate(int n, int seed, std::ostream& out)
{
	int N = n;
	int c;
	int i, j;
	int x, y;
	int nbase;
	double scale;


	/* initialize random number generator */

	sprand(seed);

	nbase = N / CLUSTERFACTOR;
	scale = SCALEFACTOR / std::sqrt((double)N);


	out << "NAME : portcgen-" << N << "-" << seed << std::endl;
	printf("NAME : portcgen-%d-%d\n", N, seed);
	out << "COMMENT : portcgen N=" << N << ", seed=" << seed << std::endl;
	printf("COMMENT : portcgen N=%d, seed=%d\n", N, seed);
	out << "TYPE : TSP" << seed << std::endl;
	printf("TYPE : TSP\n");
	out << "DIMENSION : " << N << std::endl;
	printf("DIMENSION : %d\n", N);
	out << "EDGE_WEIGHT_TYPE : EUC_2D" << std::endl;
	printf("EDGE_WEIGHT_TYPE : EUC_2D\n");
	out << "NODE_COORD_SECTION" << std::endl;
	printf("NODE_COORD_SECTION\n");



	for (i = 1; i <= nbase; i++)
		for (j = 0; j <= 1; j++)
			center[i][j] = (int)(lprand() / PRANDMAX * MAXCOORD);

	for (i = 1; i <= N; i++) {
		c = (int)(lprand() / PRANDMAX * nbase) + 1;
		x = center[c][0] + (int)(normal() * scale * MAXCOORD);
		y = center[c][1] + (int)(normal() * scale * MAXCOORD);

		out << i << " " << x << " " << y << std::endl;
		printf("%d %d %d\n", i, x, y);
	}
}

double  ofec::PortCGen::normal(void)	/* Algorithm 3.4.1.P, p. 117, Knuth v. 2 */
{
	static int	goodstill = 0;
	static double	nextvar;
	double		s, t, v1, v2;

	if (goodstill) {
		goodstill = 0;
		return nextvar;
	}
	else {
		goodstill = 1;
		do {
			v1 = 2 * lprand() / PRANDMAX - 1.0;
			v2 = 2 * lprand() / PRANDMAX - 1.0;
			s = v1 * v1 + v2 * v2;
		} while (s >= 1.0);
		t = std::sqrt((-2.0 * std::log(s)) / s);
		nextvar = v1 * t;	/* Knuth's x1 */
		return v2 * t;		/* Knuth's x2 */
	}
}

void  ofec::PortCGen::sprand(int seed)
{
	int i, ii;
	int last, next;
	arr[0] = last = seed;
	next = 1;
	for (i = 1; i < 55; i++) {
		ii = (21 * i) % 55;
		arr[ii] = next;
		next = last - next;
		if (next < 0) next += PRANDMAX;
		last = arr[ii];
	}
	a = 0;
	b = 24;
	for (i = 0; i < 165; i++) last = lprand();
}

double ofec::PortCGen::lprand(void)
{
	long t;
	if (a-- == 0) a = 54;
	if (b-- == 0) b = 54;
	t = arr[a] - arr[b];
	if (t < 0) t += PRANDMAX;
	arr[a] = t;
	return (double)t;
}

