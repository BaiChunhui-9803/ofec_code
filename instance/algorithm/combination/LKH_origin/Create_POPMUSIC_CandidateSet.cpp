#include "./INCLUDE/LKH.h"
namespace LKH {
#define maxNeighbors POPMUSIC_MaxNeighbors
#define trials POPMUSIC_Trials
#define NB_RES POPMUSIC_Solutions
#define SAMPLE_SIZE POPMUSIC_SampleSize

#define d(a, b) (((a) == (b) ? 0 : D(a, b)\
                 - (a)->Pi - (b)->Pi) / Precision)

	static void build_path(int n, int *path, int nb_clust, LKHAlg *Alg);
	static void fast_POPMUSIC(int n, int *path, int R, LKHAlg *Alg);
	static void swap(int *a, int *b);
	static int unif(int low, int high, LKHAlg *Alg);
	static void shuffle(int n, int *path, LKHAlg *Alg);
	static GainType length_path(int n, int *path, LKHAlg *Alg);
	static void path_threeOpt(int N, int **D, int *best_sol,
		GainType * best_cost, LKHAlg *Alg);

	static thread_local LKHAlg::Node **node = 0;
	static thread_local LKHAlg::Node **node_path = 0;

	static void optimize_path(int n, int *path, LKHAlg *Alg);

	/*
	 * The Create_POPMUSIC_CandidateSet function creates for each node
	 * a candidate set of K edges using the POPMUSIC algorithm.
	 *
	 * The function is called from the CreateCandidateSet function.
	 *
	 * Programmed by Keld Helsgaun and Eric Taillard, 2018.
	 *
	 * References:
	 *
	 *     É. D. Taillard and K. Helsgaun,
	 *     POPMUSIC for the Travelling Salesman Problem.
	 *     European Journal of Operational Research, 272(2):420-429 (2019).
	 *
	 *     K. Helsgaun,
	 *     Using POPMUSIC for Candidate Set Generation in the
	 *     Lin-Kernighan-Helsgaun TSP Solver.
	 *     Technical Report, Computer Science, Roskilde University, 2018.
	 */

	void LKHAlg::Create_POPMUSIC_CandidateSet(int K)
	{
		int n, i, no_res, setInitialSuc, deleted = 0, d;
		int *solution;
		GainType cost, costSum = 0;
		GainType costMin = std::numeric_limits<GainType>::max(), costMax = std::numeric_limits<GainType>::min();
		Node *N;
		double startTime, entryTime;
		int InitialTourAlgorithmSaved = InitialTourAlgorithm;

		entryTime = GetTime();
		if (TraceLevel >= 2)
			printff("Creating POPMUSIC candidate set ...\n");
		AddTourCandidates();
		if (MaxCandidates == 0) {
			N = FirstNode;
			do {
				if (!N->CandidateSet)
					eprintf("MAX_CANDIDATES = 0: No candidates");
			} while ((N = N->Suc) != FirstNode);
			if (TraceLevel >= 2)
				printff("done\n");
			return;
		}

		/* Create a tour containing all fixed or common edges */
		InitialTourAlgorithm = WALK;
		ChooseInitialTour();
		InitialTourAlgorithm = InitialTourAlgorithmSaved;
		/* N->V == 1 iff N is going to be deleted */
		N = FirstNode;
		do {
			N->V = FixedOrCommon(N, N->Pred);
		} while ((N = N->Suc) != FirstNode);

		n = Dimension;
		solution = (int *)malloc((n + 1) * sizeof(int));
		node = (Node **)malloc((n + 1) * sizeof(Node *));
		node_path = (Node **)malloc((n + 1) * sizeof(Node *));

		for (no_res = 1; no_res <= NB_RES; no_res++) {
			/* Create set of non-deleted nodes */
			n = 0;
			N = FirstNode;
			do {
				if (!N->V) {
					node[solution[n] = n] = N;
					n++;
				}
				else
					deleted++;
			} while ((N = N->Suc) != FirstNode);
			shuffle(n, solution, this);
			solution[n] = solution[0];
			node[n] = node[solution[0]];
			startTime = GetTime();
			build_path(n, solution, SAMPLE_SIZE, this);
			if (deleted > 0) {
				/* Create a one-way list representing the built tour */
				for (i = 1; i <= n; i++)
					node[solution[i - 1]]->Next = node[solution[i]];
				/* Insert the deleted nodes in the one-way list */
				N = node[solution[0]];
				do {
					if (!N->V && N->Suc->V) {
						Node *OldNext = N->Next;
						do
							N = N->Next = N->Suc;
						while (N->Suc->V);
						N->Next = OldNext;
					}
				} while ((N = N->Suc) != node[solution[0]]);
				/* Convert the one-way list to a solution array */
				N = node[solution[0]];
				n = 0;
				do {
					node[solution[n] = n] = N;
					n++;
				} while ((N = N->Next) != node[solution[0]]);
			}
			solution[n] = solution[0];
			node[n] = node[solution[0]];
			cost = length_path(n, solution, this);
			if (TraceLevel >= 2) {
				printff("%d: Initial cost:  %lld, ", no_res, cost);
				if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
					printff("Gap = %0.2f%%, ",
						100.0 * (cost - Optimum) / Optimum);
				printff("Time: %0.2f sec.\n", GetTime() - startTime);
			}
			startTime = GetTime();
			fast_POPMUSIC(n, solution, SAMPLE_SIZE * SAMPLE_SIZE, this);
			solution[n] = solution[0];
			node[n] = node[solution[0]];
			cost = length_path(n, solution, this);
			if (TraceLevel >= 2) {
				printff("%d: Improved cost: %lld, ", no_res, cost);
				if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
					printff("Gap = %0.2f%%, ",
						100.0 * (cost - Optimum) / Optimum);
				printff("Time: %0.2f sec.\n", GetTime() - startTime);
			}
			costSum += cost;
			if (cost > costMax)
				costMax = cost;
			setInitialSuc = 0;
			if (cost < costMin) {
				costMin = cost;
				setInitialSuc = POPMUSIC_InitialTour && !FirstNode->InitialSuc;
			}
			for (i = 0; i < n; i++) {
				Node *a = node[solution[i]];
				Node *b = node[solution[i + 1]];
				d = (this->*D)(a, b);
				AddCandidate(a, b, d, 1);
				AddCandidate(b, a, d, 1);
				if (setInitialSuc)
					a->InitialSuc = b;
			}
		}
		if (TraceLevel >= 2) {
			printff
			("Cost.min = " GainFormat ", Cost.avg = %0.2f, Cost.max = "
				GainFormat "\n", costMin, (double)costSum / NB_RES, costMax);
			if (Optimum != std::numeric_limits<GainType>::min() && Optimum != 0)
				printff
				("Gap.min = %0.2f%%, Gap.avg = %0.2f%%, Gap.max = %0.2f%%\n",
					100.0 * (costMin - Optimum) / Optimum,
					100.0 * ((double)costSum / NB_RES - Optimum) / Optimum,
					100.0 * (costMax - Optimum) / Optimum);
		}
		free(solution);
		free(node);
		free(node_path);
		ResetCandidateSet();
		if (K > 0)
			TrimCandidateSet(K);
		AddTourCandidates();
		if (CandidateSetSymmetric)
			SymmetrizeCandidateSet();
		if (TraceLevel >= 2) {
			CandidateReport();
			printff("POPMUSIC Time = %0.2f sec.\n", GetTime() - entryTime);
			printff("done\n");
		}
	}

	/************************ Compute the length of a path ************************/
	static GainType length_path(int n, int *path, LKHAlg *Alg)
	{
		LKHAlg::Node *a, *b;
		GainType length = 0;
		int i;

		for (i = 1, a = node[path[0]]; i <= n; i++, a = b) {
			b = node[path[i]];
			length += (((a) == (b) ? 0 : (Alg->*(Alg->D))(a, b)\
				- (a)->Pi - (b)->Pi) / Alg->Precision);
		}
		return length;
	}

	/*************** Optimize the path between path[0] and path[n] ****************/
	static void optimize_path(int n, int *path, LKHAlg *Alg)
	{
		int i, j;
		int *order, **d;
		GainType length;
		LKHAlg::Node *a, *b;

		order = (int *)malloc((n + 1) * sizeof(int));
		d = (int **)malloc((n + 1) * sizeof(int *));
		for (i = 0; i <= n; i++) {
			order[i] = i;
			d[i] = (int *)malloc((n + 1) * sizeof(int));
		}
		for (i = 0; i < n; i++) {
			a = node[path[i]];
			for (j = i + 1; j <= n; j++) {
				b = node[path[j]];
				d[i][j] = d[j][i] = (((a) == (b) ? 0 : (Alg->*(Alg->D))(a, b)\
					- (a)->Pi - (b)->Pi) / Alg->Precision);
			}
		}
		d[0][n] = d[n][0] = 0;
		length = 0;
		for (i = 0; i < n; i++)
			length += d[i][i + 1];
		for (i = 0; i <= n; i++)
			node_path[i] = node[path[i]];
		path_threeOpt(n, d, order, &length, Alg);
		for (i = 0; i <= n; i++)
			free(d[i]);
		free(d);
		for (i = 0; i <= n; i++)
			order[i] = path[order[i]];
		for (i = 0; i <= n; i++)
			path[i] = order[i];
		free(order);
	}

	/*********** Build recursively a path between path[0] and path[n] *************/
	static void build_path(int n, int *path, int nb_clust, LKHAlg *Alg)
	{
		int i, j, d, dmin, closest, start, end;
		int *tmp_path, *sample, *assignment, *start_clust, *assigned;

		if (n <= 2)
			return;
		if (n <= nb_clust * nb_clust) {
			optimize_path(n, path, Alg);
			return;
		}

		/* Temporary path */
		tmp_path = (int *)malloc((n + 1) * sizeof(int));
		for (i = 0; i <= n; i++)
			tmp_path[i] = path[i];

		/* Set tmp_path[1] be the closest city to path[0] */
		dmin = (((node[tmp_path[1]]) == (node[path[0]]) ? 0 : (Alg->*(Alg->D))(node[tmp_path[1]], node[path[0]])\
			- (node[tmp_path[1]])->Pi - (node[path[0]])->Pi) / Alg->Precision);
		closest = 1;
		for (i = 2; i < n; i++) {
			d = (((node[tmp_path[i]]) == (node[path[0]]) ? 0 : (Alg->*(Alg->D))(node[tmp_path[i]], node[path[0]])\
				- (node[tmp_path[i]])->Pi - (node[path[0]])->Pi) / Alg->Precision);
			if (d < dmin) {
				dmin = d;
				closest = i;
			}
		}
		swap(tmp_path + 1, tmp_path + closest);

		/* Set tmp_path[2] be the closest city to path[n] */
		dmin = (((node[tmp_path[2]]) == (node[path[n]]) ? 0 : (Alg->*(Alg->D))(node[tmp_path[2]], node[path[n]])\
			- (node[tmp_path[2]])->Pi - (node[path[n]])->Pi) / Alg->Precision);
		closest = 2;
		for (i = 3; i < n; i++) {
			d = (((node[tmp_path[i]]) == (node[path[n]]) ? 0 : (Alg->*(Alg->D))(node[tmp_path[i]], node[path[n]])\
				- (node[tmp_path[i]])->Pi - (node[path[n]])->Pi) / Alg->Precision);
			if (d < dmin) {
				dmin = d;
				closest = i;
			}
		}
		swap(tmp_path + 2, tmp_path + closest);

		/* Choose a sample of nb_clust-2 random cities from tmp_path */
		for (i = 3; i <= nb_clust; i++)
			swap(tmp_path + i, tmp_path + unif(i, n - 1, Alg));
		swap(tmp_path + 2, tmp_path + nb_clust);
		sample = (int *)malloc((nb_clust + 2) * sizeof(int));
		for (i = 0; i <= nb_clust; i++)
			sample[i] = tmp_path[i];
		sample[nb_clust + 1] = sample[0];
		optimize_path(nb_clust + 1, sample, Alg);

		/* Assign each city of path to the closest of the sample */
		assignment = (int *)malloc((n + 1) * sizeof(int));
		for (i = 1; i < n; i++) {
			dmin = std::numeric_limits<int>::max();
			for (j = 1; j <= nb_clust; j++) {
				if (path[i] == sample[j]) {
					closest = j;
					break;
				}
				d = (((node[path[i]]) == (node[sample[j]]) ? 0 : (Alg->*(Alg->D))(node[path[i]], node[sample[j]])\
					- (node[path[i]])->Pi - (node[sample[j]])->Pi) / Alg->Precision);
				if (d < dmin) {
					dmin = d;
					closest = j;
				}
			}
			assignment[i] = closest;
		}

		/* Build clusters: ith cluster has (start_clust[i+1]-start_clust[i])
		   cities */
		start_clust = (int *)calloc(nb_clust + 1, sizeof(int));
		for (i = 1; i < n; i++)
			start_clust[assignment[i]]++;
		for (i = 1; i <= nb_clust; i++)
			start_clust[i] += start_clust[i - 1];

		/* Clusters are stored in tmp_path in the order given by sample */
		assigned = (int *)calloc(nb_clust + 1, sizeof(int));
		for (i = 1; i < n; i++)
			tmp_path[start_clust[assignment[i] - 1] +
			assigned[assignment[i]]++] = path[i];

		/* Reorder original path */
		for (i = 1; i < n; i++)
			path[i] = tmp_path[i - 1];
		free(tmp_path);
		free(sample);
		free(assigned);
		free(assignment);
		for (i = 0; i < nb_clust; i++) {
			start = start_clust[i] - 1;
			if (start < 0)
				start = 0;
			end = start_clust[i + 1] + 1;
			if (end > n)
				end = n;
			/* Recursively optimize sub-path corresponding to each cluster */
			if (end - start < n)
				build_path(end - start, path + start, nb_clust, Alg);
		}
		free(start_clust);
	}

	static void reverse(int *path, int i, int j)
	{
		while (i < j)
			swap(path + i++, path + j--);
	}

	static void circular_right_shift(int n, int *path, int positions)
	{
		reverse(path, 0, positions - 1);
		reverse(path, positions, n - 1);
		reverse(path, 0, n - 1);
	}

	/* Fast POPMUSIC: Optimize independent subpaths of R cities */
	static void fast_POPMUSIC(int n, int *path, int R, LKHAlg *Alg)
	{
		int scan, i;

		if (R > n)
			R = n;
		/* Optimize subpaths of R cities with R/2 overlap; 2 scans */
		for (scan = 1; scan <= 2; scan++) {
			if (scan == 2) {
				circular_right_shift(n, path, R / 2);
				path[n] = path[0];
				node[n] = node[path[0]];
			}
			for (i = 0; i < n / R; i++)
				optimize_path(R, path + R * i, Alg);
			if (n % R != 0) /* Optimize last portion of the path */
				optimize_path(R, path + n - R, Alg);
		}
	}

	/* Iterated 3-opt code  */

	static thread_local int n__;
	static thread_local int **dist;
	static thread_local int *tour_, *pos;
	static thread_local int **neighbor;
	static thread_local int *neighbors;
	static thread_local int reversed;
	static thread_local char *dontLook;
	static thread_local GainType tourLength;

	static void createNeighbors(LKHAlg *Alg);
	static void threeOpt(LKHAlg *Alg);
	static void doubleBridgeKick(LKHAlg *Alg);
	static int prev(int v);
	static int next(int v);
	static int PREV(int v);
	static int NEXT(int v);

	static void path_threeOpt(int N, int **D, int *best_sol,
		GainType * best_cost, LKHAlg *Alg)
	{
		int i, j, trial;
		int *bestTour;
		GainType bestTourLength;

		n__ = N + 1;
		tour_ = (int *)malloc(n__ * sizeof(int));
		pos = (int *)malloc(n__ * sizeof(int));
		bestTour = (int *)malloc(n__ * sizeof(int));
		for (i = 0; i < n__; i++)
			pos[bestTour[i] = tour_[i] = best_sol[i]] = i;
		dist = D;
		createNeighbors(Alg);
		dontLook = (char *)calloc(n__, sizeof(char));
		bestTourLength = tourLength = *best_cost;
		if (Alg->POPMUSIC_Trials == 0)
			Alg->POPMUSIC_Trials = n__;
		for (trial = 1; trial <= Alg->POPMUSIC_Trials; trial++) {
			threeOpt(Alg);
			if (tourLength < bestTourLength) {
				for (i = 0; i < n__; i++)
					bestTour[i] = tour_[i];
				bestTourLength = tourLength;
			}
			else {
				for (i = 0; i < n__; i++) {
					pos[tour_[i] = bestTour[i]] = i;
					tourLength = bestTourLength;
				}
			}
			if (n__ <= 5 || trial == Alg->POPMUSIC_Trials)
				break;
			doubleBridgeKick(Alg);
		}
		*best_cost = bestTourLength;
		for (i = 0; i < n__; i++)
			pos[tour_[i] = bestTour[i]] = i;
		reversed = next(0) == N;
		for (i = 0, j = 0; j < n__; i = NEXT(i), j++)
			best_sol[j] = i;
		free(tour_);
		free(pos);
		free(bestTour);
		free(neighbors);
		for (i = 0; i < n__; i++)
			free(neighbor[i]);
		free(neighbor);
		free(dontLook);
	}

	static int unif(int low, int high, LKHAlg *Alg)
	{
		return low + Alg->Random() % (high - low + 1);
	}

	static void shuffle(int n, int *path, LKHAlg *Alg)
	{
		int i;
		for (i = 1; i < n; i++)
			swap(path + i, path + unif(0, i, Alg));
	}

	static void swap(int *a, int *b)
	{
		int tmp = *a;
		*a = *b;
		*b = tmp;
	}

	static int fixed(int a, int b, LKHAlg *Alg)
	{
		return (a == 0 && b == n__ - 1) || (a == n__ - 1 && b == 0) ||
			//FixedOrCommon(node_path[a], node_path[b]);
			(LKH::LKHAlg::Fixed(node_path[a], node_path[b]) || Alg->IsCommonEdge(node_path[a], node_path[b]));
	}

	static int prev(int v)
	{
		return tour_[pos[v] > 0 ? pos[v] - 1 : n__ - 1];
	}

	static int next(int v)
	{
		return tour_[pos[v] < n__ - 1 ? pos[v] + 1 : 0];
	}

	static int between(int v1, int v2, int v3)
	{
		int a = pos[v1], b = pos[v2], c = pos[v3];
		return a <= c ? b >= a && b <= c : b <= c || b >= a;
	}

	static int PREV(int v)
	{
		return reversed ? next(v) : prev(v);
	}

	static int NEXT(int v)
	{
		return reversed ? prev(v) : next(v);
	}

	static int BETWEEN(int v1, int v2, int v3)
	{
		return !reversed ? between(v1, v2, v3) : between(v3, v2, v1);
	}

	static void flip(int from, int to)
	{
		int i, j, size, tmp;

		if (from == to)
			return;
		if (reversed) {
			tmp = from;
			from = to;
			to = tmp;
		}
		i = pos[from], j = pos[to];
		size = j - i;
		if (size < 0)
			size += n__;
		if (size >= n__ / 2) {
			tmp = i;
			i = ++j < n__ ? j : 0;
			j = --tmp >= 0 ? tmp : n__ - 1;
		}
		while (i != j) {
			tmp = tour_[i];
			pos[tour_[i] = tour_[j]] = i;
			pos[tour_[j] = tmp] = j;
			if (++i == n__)
				i = 0;
			if (i != j && --j < 0)
				j = n__ - 1;
		}
	}

	static void threeOpt(LKHAlg *Alg)
	{
		int improved = 1, a, b, c, d, e, f, xa, xc, xe, i, j;
		GainType g0, g1, g2, g3, gain;


		while (improved) {
			improved = 0;
			for (b = 0; b < n__; b++) {
				if (dontLook[b])
					continue;
				dontLook[b] = 1;
				for (xa = 1; xa <= 2; xa++, reversed = !reversed) {
					a = PREV(b);
					if (fixed(a, b, Alg))
						continue;
					g0 = dist[a][b];
					for (i = 0; i < neighbors[b]; i++) {
						c = neighbor[b][i];
						if (c == prev(b) || c == next(b))
							continue;
						g1 = g0 - dist[b][c];
						if (g1 <= 0)
							break;
						for (xc = 1; xc <= 2; xc++) {
							d = xc == 1 ? PREV(c) : NEXT(c);
							if (d == a || fixed(c, d, Alg))
								continue;
							g2 = g1 + dist[c][d];
							if (xc == 1) {
								gain = g2 - dist[d][a];
								if (gain > 0) {
									flip(b, d);
									tourLength -= gain;
									dontLook[a] = dontLook[b] = 0;
									dontLook[c] = dontLook[d] = 0;
									improved = 1;
									i = neighbors[b];
									break;
								}
							}
							for (j = 0; j < neighbors[d]; j++) {
								e = neighbor[d][j];
								if (e == prev(d) || e == next(d) ||
									(xc == 2 && !BETWEEN(b, e, c)))
									continue;
								g3 = g2 - dist[d][e];
								if (g3 <= 0)
									break;
								for (xe = 1; xe <= xc; xe++) {
									if (xc == 1)
										f = BETWEEN(b, e,
											c) ? NEXT(e) : PREV(e);
									else
										f = xe == 1 ? PREV(e) : NEXT(e);
									if (f == a || fixed(e, f, Alg))
										continue;
									gain = g3 + dist[e][f] - dist[f][a];
									if (gain > 0) {
										if (xc == 1) {
											flip(b, d);
											if (f == PREV(e))
												flip(e, a);
											else
												flip(a, e);
										}
										else if (xe == 1) {
											flip(e, c);
											if (b == NEXT(a))
												flip(b, f);
											else
												flip(f, b);
										}
										else {
											flip(d, a);
											if (f == NEXT(e))
												flip(f, c);
											else
												flip(c, f);
											if (b == NEXT(d))
												flip(b, e);
											else
												flip(e, b);
										}
										tourLength -= gain;
										dontLook[a] = dontLook[b] = 0;
										dontLook[c] = dontLook[d] = 0;
										dontLook[e] = dontLook[f] = 0;
										improved = 1;
										goto Next_b;
									}
								}
							}
						}
					}
				Next_b:;
				}
			}
		}
	}

	static void createNeighbors(LKHAlg *Alg)
	{
		int i, j, k, d;

		neighbor = (int **)malloc(n__ * sizeof(int *));
		for (i = 0; i < n__; i++)
			neighbor[i] = (int *)malloc((Alg->POPMUSIC_MaxNeighbors + 1) * sizeof(int));
		neighbors = (int *)calloc(n__, sizeof(int));
		for (i = 0; i < n__; i++) {
			for (j = 0; j < n__; j++) {
				if (i == j)
					continue;
				d = dist[i][j];
				k = neighbors[i] <
					Alg->POPMUSIC_MaxNeighbors ? neighbors[i]++ : Alg->POPMUSIC_MaxNeighbors;
				while (k > 0 && d < dist[i][neighbor[i][k - 1]]) {
					neighbor[i][k] = neighbor[i][k - 1];
					k--;
				}
				neighbor[i][k] = j;
			}
		}
	}

	static void doubleBridgeKick(LKHAlg *Alg)
	{
		int t[4], r[4], a, b, c, d, e, f, g, h, i, j;

		reversed = 0;
		for (i = 0; i <= 3; i++) {
			r[i] = unif(0, n__ - i - 1, Alg);
			for (j = 0; j < n__ - i; j++) {
				t[i] = tour_[r[i]];
				if (!fixed(t[i], next(t[i]), Alg)) {
					swap(tour_ + r[i], tour_ + n__ - 1 - i);
					break;
				}
				if (++r[i] == n__ - i)
					r[i] = 0;
			}
			if (j == n__ - i)
				return;
		}
		for (i = 3; i >= 0; i--)
			swap(tour_ + r[i], tour_ + n__ - 1 - i);
		if (pos[t[0]] > pos[t[1]])
			swap(t, t + 1);
		if (pos[t[2]] > pos[t[3]])
			swap(t + 2, t + 3);
		if (pos[t[0]] > pos[t[2]])
			swap(t, t + 2);
		if (pos[t[1]] > pos[t[3]])
			swap(t + 1, t + 3);
		if (pos[t[1]] > pos[t[2]])
			swap(t + 1, t + 2);
		a = t[0];
		b = next(a);
		c = t[2];
		d = next(c);
		e = t[1];
		f = next(e);
		g = t[3];
		h = next(g);
		flip(b, c);
		if (f == next(e))
			flip(f, h);
		else
			flip(h, f);
		if (b == next(d))
			flip(c, d);
		else
			flip(d, c);
		dontLook[a] = dontLook[b] = dontLook[c] = dontLook[d] = 0;
		dontLook[e] = dontLook[f] = dontLook[g] = dontLook[h] = 0;
		tourLength -= dist[a][b] - dist[b][c] +
			dist[c][d] - dist[d][a] +
			dist[e][f] - dist[f][g] + dist[g][h] - dist[h][e];
	}
}