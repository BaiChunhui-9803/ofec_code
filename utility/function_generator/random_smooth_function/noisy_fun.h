
#ifndef NOISY_FUN_H

#define NOISY_FUN_H

namespace ofec {
	namespace noisy_fun {

		static unsigned char m_perm[512] = { 151, 160, 137, 91, 90, 15,
		131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
		190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
		88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
		77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
		102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
		135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
		5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
		223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
		129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
		251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
		49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
		138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180,
		151, 160, 137, 91, 90, 15,
		131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
		190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
		88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
		77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
		102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
		135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
		5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
		223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
		129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
		251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
		49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
		138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180
		};




		static double m_grad3(int hash, double x, double y, double z) {
			int h = hash & 15;     /* Convert low 4 bits of hash code into 12 simple */
			double u = h < 8 ? x : y; /* gradient directions, and compute dot product. */
			double v = h < 4 ? y : h == 12 || h == 14 ? x : z; /* Fix repeats at h = 12 to 15 */
			return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
		}

#define m_FASTFLOOR(x) ( ((x)>0) ? ((int)x) : (((int)x)-1) )



		/* 3D simplex noise */
		static double m_slang_library_noise3(double x, double y, double z)
		{
			/* Simple skewing factors for the 3D case */
#define m_F3 0.333333333f
#define m_G3 0.166666667f

			double n0, n1, n2, n3; /* Noise contributions from the four corners */

								   /* Skew the input space to determine which simplex cell we're in */
			double s = (x + y + z) * m_F3; /* Very nice and simple skew factor for 3D */
			double xs = x + s;
			double ys = y + s;
			double zs = z + s;
			int i = m_FASTFLOOR(xs);
			int j = m_FASTFLOOR(ys);
			int k = m_FASTFLOOR(zs);

			double t = (double)(i + j + k) * m_G3;
			double X0 = i - t; /* Unskew the cell origin back to (x,y,z) space */
			double Y0 = j - t;
			double Z0 = k - t;
			double x0 = x - X0; /* The x,y,z distances from the cell origin */
			double y0 = y - Y0;
			double z0 = z - Z0;

			double x1, y1, z1, x2, y2, z2, x3, y3, z3;
			int ii, jj, kk;
			double t0, t1, t2, t3;

			/* For the 3D case, the simplex shape is a slightly irregular tetrahedron. */
			/* Determine which simplex we are in. */
			int i1, j1, k1; /* Offsets for second corner of simplex in (i,j,k) coords */
			int i2, j2, k2; /* Offsets for third corner of simplex in (i,j,k) coords */

							/* This code would benefit from a backport from the GLSL version! */
			if (x0 >= y0) {
				if (y0 >= z0)
				{
					i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0;
				} /* X Y Z order */
				else if (x0 >= z0) { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1; } /* X Z Y order */
				else { i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1; } /* Z X Y order */
			}
			else { /* x0<y0 */
				if (y0 < z0) { i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1; } /* Z Y X order */
				else if (x0 < z0) { i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1; } /* Y Z X order */
				else { i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0; } /* Y X Z order */
			}

			/* A step of (1,0,0) in (i,j,k) means a step of (1-c,-c,-c) in (x,y,z), */
			/* a step of (0,1,0) in (i,j,k) means a step of (-c,1-c,-c) in (x,y,z), and */
			/* a step of (0,0,1) in (i,j,k) means a step of (-c,-c,1-c) in (x,y,z), where */
			/* c = 1/6. */

			x1 = x0 - i1 + m_G3; /* Offsets for second corner in (x,y,z) coords */
			y1 = y0 - j1 + m_G3;
			z1 = z0 - k1 + m_G3;
			x2 = x0 - i2 + 2.0f * m_G3; /* Offsets for third corner in (x,y,z) coords */
			y2 = y0 - j2 + 2.0f * m_G3;
			z2 = z0 - k2 + 2.0f * m_G3;
			x3 = x0 - 1.0f + 3.0f * m_G3; /* Offsets for last corner in (x,y,z) coords */
			y3 = y0 - 1.0f + 3.0f * m_G3;
			z3 = z0 - 1.0f + 3.0f * m_G3;

			/* Wrap the integer indices at 256, to avoid indexing perm[] out of bounds */
			ii = i % 256;
			jj = j % 256;
			kk = k % 256;

			/* Calculate the contribution from the four corners */
			t0 = 0.6f - x0 * x0 - y0 * y0 - z0 * z0;
			if (t0 < 0.0f) n0 = 0.0f;
			else {
				t0 *= t0;
				n0 = t0 * t0 * m_grad3(m_perm[ii + m_perm[jj + m_perm[kk]]], x0, y0, z0);
			}

			t1 = 0.6f - x1 * x1 - y1 * y1 - z1 * z1;
			if (t1 < 0.0f) n1 = 0.0f;
			else {
				t1 *= t1;
				n1 = t1 * t1 * m_grad3(m_perm[ii + i1 + m_perm[jj + j1 + m_perm[kk + k1]]], x1, y1, z1);
			}

			t2 = 0.6f - x2 * x2 - y2 * y2 - z2 * z2;
			if (t2 < 0.0f) n2 = 0.0f;
			else {
				t2 *= t2;
				n2 = t2 * t2 * m_grad3(m_perm[ii + i2 + m_perm[jj + j2 + m_perm[kk + k2]]], x2, y2, z2);
			}

			t3 = 0.6f - x3 * x3 - y3 * y3 - z3 * z3;
			if (t3 < 0.0f) n3 = 0.0f;
			else {
				t3 *= t3;
				n3 = t3 * t3 * m_grad3(m_perm[ii + 1 + m_perm[jj + 1 + m_perm[kk + 1]]], x3, y3, z3);
			}

			/* Add contributions from each corner to get the final noise value. */
			/* The result is scaled to stay just inside [-1,1] */
			return 32.0f * (n0 + n1 + n2 + n3); /* TODO: The scale factor is preliminary! */
		}


		/* The result is scaled to stay just inside [0,1] */
		extern double m_ofNoise(double x, double y, double z);

	}
}



#endif 