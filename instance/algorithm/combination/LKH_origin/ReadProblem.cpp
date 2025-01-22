#include "./INCLUDE/LKH.h"
#include "./INCLUDE/Heap.h"
#include <stdio.h>
#include <string.h>
/*
 * The ReadProblem function reads the problem data in TSPLIB format from the
 * file specified in the parameter file (PROBLEM_FILE).
 *
 * The following description of the file format is extracted from the TSPLIB
 * documentation.
 *
 * The file consists of a specification part and a data part. The specification
 * part contains information on the file format and on its contents. The data
 * part contains explicit data.
 *
 * (1) The specification part
 *
 * All entries in this section are of the form <keyword> : <value>, where
 * <keyword> denotes an alphanumerical keyword and <value> denotes
 * alphanumerical or numerical data. The terms <string>, <integer> and <Real>
 * denote character string, integer or Real data, respectively. The order of
 * specification of the keywords in the data file is arbitrary (in principle),
 * but must be consistent, i.e., whenever a keyword is specified, all
 * necessary information for the correct interpretation of the keyword has to
 * be known.
 *
 * Below is given a list of all available keywords.
 *
 * NAME : <string>e
 * Identifies the data file.
 *
 * TYPE : <string>
 * Specifies the type of data. Possible types are
 * TSP          Data for a symmetric traveling salesman problem
 * ATSP         Data for an asymmetric traveling salesman problem
 * SOP          Data for a sequence ordering problem
 * HCP          Hamiltonian cycle problem data
 * HPP          Hamiltonian path problem data (not available in TSPLIB)
 * TSPTW        Data for a TSP instance with time windows
 * CCVRP        Data for a cumulative capacitated vehicle routing problem
 * CVRP         Data for a symmetric capacitated vehicle routing problem
 * ACVRP        Data for an asymmetric capacitated vehicle routing problem
 * CVRPTW       Data for a capacitated vehicle routing problem with
 *              time windows
 * VRPMPD       Data for a mixed pickup and delivery problem with backhauls
 * 1-PDTSP      Data for a one-commodity pickup-and-delivery traveling
 *              salesman problem
 * MLP          Data for a minimum latency problem
 * m-PDTSP      Data for a mulity-commodity pickup-and-delivery traveling
 *              salesman problem
 * m1-PDTSP     Data for a mulity-commodity one-to-one pickup-and-delivery
 *              traveling salesman problem
 * OVRP         Data for an open vehicle routing problem
 * PDTSP        Data for a pickup and delivery traveling salesman problem
 * PDTSPF       Data for a pickup and delivery traveling salesman problem
 *              with FIFO loading
 * PDTSPL       Data for a pickup and delivery traveling salesman problem
 *              with LIFO loading
 * PDPTW        Data for a pickup and delivery problem with time windows
 * RCTVRP       Data for a risk-constrained cash-in-transit vehicle
 *              routing problem
 * RCTVRPTW     Data for a risk-constrained cash-in-transit vehicle
 *              routing problem with time windows
 * TRP          Data for a traveling repairman problem
 * TSPDL        Dara for a traveling salesman problem with draft limits
 * TSPPD        Data for a pickup and delivery travling salesman problem
 * VRPB         Data for a vehicle routing problem with backhauls
 * VRPBTW       Data for a vehicle routing problem with backhauls and
 *              time windows
 * CTSP         Data for a colored traveling salesman problem
 *
 * COMMENT : <string>
 * Additional comments (usually the name of the contributor or the creator of
 * the problem instance is given here).
 *
 * DIMENSION : < integer>
 * The number of nodes.
 *
 * CAPACITY : <integer>
 * Specifies the truck capacity in a CVRP.
 *
 * DISTANCE : <Real>
 * The maximum length allowed for each route in a CVRP.
 *
 * EDGE_WEIGHT_TYPE : <string>
 * Specifies how the edge weights (or distances) are given. The values are:
 * ATT          Special distance function for problem att48 and att532
 * CEIL_2D      Weights are Euclidean distances in 2-D rounded up
 * CEIL_3D      Weights are Euclidean distances in 3-D rounded up
 * EUC_2D       Weights are Euclidean distances in 2-D
 * EUC_3D       Weights are Euclidean distances in 3-D
 * EXACT_2D     Weights are EUC_2D distances (SCALE = 1000 as default)
 * EXACT_3D     Weights are EUC_3D distances (SCALE = 1000 as default)
 * EXPLICIT     Weights are listed explicitly in the corresponding section
 * FLOOR_2D     Weights are Euclidean distances in 2-D rounded down
 * FLOOR_3D     Weights are Euclidean distances in 3-D rounded down
 * GEO          Weights are geographical distances in kilometers (TSPLIB)
 *              Coordinates are given in the form DDD.MM where DDD are the
 *              degrees and MM the minutes
 * GEOM         Weights are geographical distances in meters (used for the
 *              world TSP). Coordinates are given in decimal form
 * GEO_MEEUS    Weights are geographical distances in kilometers, computed
 *              according to Meeus' formula.  Coordinates are given in the
 *              form DDD.MM where DDD are the degrees and MM the minutes
 * GEOM_MEEUS   Weights are geographical distances, computed according to
 *              Meeus' formula. Coordinates are given in decimal form
 * MAN_2D       Weights are Manhattan distances in 2-D
 * MAN_3D       Weights are Manhattan distances in 3-D
 * MAX_2D       Weights are maximum distances in 2-D
 * MAX_3D       Weights are maximum distances in 3-D
 * TOR_2D       Wirghes are toroidal distances in 2-D
 * TOR_3D       Wirghes are toroidal distances in 3-D
 * XRAY1        Distance function for crystallography problems (Version 1)
 * XRAY2        Distance function for crystallography problems (Version 2)
 * SPECIAL      There is a special distance function implemented in
 *              the Distance_SPECIAL function.
 *
 * EDGE-WEIGHT_FORMAT : <string>
 * Describes the format of the edge weights if they are given explicitly.
 * The values are
 * FUNCTION         Weights are given by a function (see above)
 * FULL_MATRIX      Weights are given by a full matrix
 * UPPER_ROW        Upper triangular matrix
 *                      (row-wise without diagonal entries)
 * LOWER_ROW        Lower triangular matrix
 *                      (row-wise without diagonal entries)
 * UPPER_DIAG_ROW   Upper triangular matrix
 *                      (row-wise including diagonal entries)
 * LOWER_DIAG_ROW   Lower triangular matrix
 *                      (row-wise including diagonal entries)
 * UPPER_COL        Upper triangular matrix
 *                      (column-wise without diagonal entries)
 * LOWER_COL        Lower triangular matrix
 *                      (column-wise without diagonal entries)
 * UPPER_DIAG_COL   Upper triangular matrix
 *                      (column-wise including diagonal entries)
 * LOWER_DIAG_COL   Lower triangular matrix
 *                      (column-wise including diagonal entries)
 *
 * EDGE_DATA_FORMAT : <string>
 * Describes the format in which the edges of a graph are given, if the
 * graph is not complete. The values are
 * EDGE_LIST    The graph is given by an edge list
 * ADJ_LIST     The graph is given by an adjacency list
 *
 * NODE_COORD_TYPE : <string>
 * Specifies whether the coordinates are associated with each node
 * (which, for example may be used for either graphical display or
 * distance computations.
 * The values are
 * TWOD_COORDS      Nodes are specified by coordinates in 2-D
 * THREED_COORDS    Nodes are specified by coordinates in 3-D
 * NO_COORDS        The nodes do not have associated coordinates
 * The default value is NO_COORDS. In the current implementation, however,
 * the value has no significance.
 *
 * DISPLAY_DATA_TYPE : <string>
 * Specifies how a graphical display of the nodes can be obtained.
 * The values are
 * COORD_DISPLAY    Display is generated from the node coordinates
 * TWOD_DISPLAY     Explicit coordinates in 2-D are given
 * NO_DISPLAY       No graphical display is possible
 *
 * The default value is COORD_DISPLAY if node coordinates are specifies and
 * NO_DISPLAY otherwise. In the current implementation, however, the value
 * has no significance.
 *
 * DEMAND_DIMENSION : < integer>
 * The number of objects in an m1-PDTSP.
 *
 * GRID_SIZE : <Real>
 * The grid size for toroidal instances.
 * Default: 1000000.0
 *
 * RISK_THRESHOLD : <integer>
 * The maximum risk alllowed for each route in an RCTVRP or RCTVRPTW instance.
 *
 * SALESMEN : <integer>
 * VEHICLES : <integer>
 * The number of vehicles/salesmen in a CVRP.
 *
 * SCALE : <integer>
 * Scale factor. Distances are multiplied by this factor.
 *
 * EOF
 * Terminates input data. The entry is optional.
 *
 * (2) The data part
 *
 * Depending on the choice of specifications some additional data may be
 * required. These data are given corresponding data sections following the
 * specification part. Each data section begins with the corresponding
 * keyword. The length of the sectionis either explicitly known form the
 * format specification, or the section is terminated by an appropriate
 * end-of-section identifier.
 *
 * NODE_COORD_SECTION :
 * Node coordinates are given in this section. Each line is of the form
 *
 *      <integer> <Real> <Real>
 *
 * if NODE_COORD_TYPE is TWOD_COORDS, or
 *
 *      <integer> <Real> <Real> <Real>
 *
 * if NODE_COORD_TYPE is THREED_COORDS. The integers give the number of the
 * respective nodes. The Real numbers are the associated coordinates.
 *
 * EDGE_DATA_SECTION :
 * Edges of the graph are specified in either of the two formats allowed in
 * the EDGE_DATA_FORMAT entry. If the type is EDGE_LIST, then the edges are
 * given as a sequence of lines of one of the forms
 *
 *      <integer> <integer>
 *      <integer> <integer> <integer>
 *
 * each entry giving the terminal nodes of some edge, and if three integers are
 * given, the last one specifies its weight. The list is terminated y a -1.
 * If the type is ADJ_LIST, the section consists of adjacency lists for nodes.
 * The adjacency list of a node x is specified as
 *
 *      <integer> <integer> ... <integer> -1
 *
 * where the first integer gives the number of node x and the following
 * integers (terminated by -1) the numbers of the nodes adjacent to x.
 * The list of adjacency lists are terminated by an additional -1.
 *
 * FIXED_EDGES_SECTION :
 * In this section, edges are listed that are required to appear in each
 * solution to the problem. The edges to be fixed are given in the form
 * (per line)
 *
 *      <integer> <integer>
 *
 * meaning that the edge (arc) from the first node to the second node has
 * to be contained in a solution. This section is terminated by a -1.
 *
 * DISPLAY_DATA_SECTION :
 * If DISPLAY_DATA_TYPE is TWOD_DISPLAY, the 2-dimensional coordinates from
 * which a display can be generated are given in the form (per line)
 *
 *      <integer> <Real> <Real>
 *
 * The integers specify the respective nodes and the Real numbers give the
 * associated coordinates. The contents of this section, however, has no
 * significance in the current implementation.
 *
 * EDGE_WEIGHT_SECTION :
 * The edge weights are given in the format specifies by the EDGE_WEIGHT_FORMAT
 * entry. At present, all explicit data are integral and is given in one of the
 * (self-explanatory) matrix formats, with explicitly known lengths.
 *
 * TOUR_SECTION :
 * A tour is specified in this section. The tour is given by a list of
 * integers giving the sequence in which the nodes are visited in the tour.
 * The tour is terminated by a -1. Note: In contrast to the TSPLIB format,
 * only one tour can be given in this section. The tour is used to limit
 * the search (the last edge to be excluded in a non-gainful move must not
 * belong to the tour). In addition, the Alpha field of its edges is set to
 * -1.
 *
 * BACKHAUL_SECTION :
 * This section is used for specifying VRPB instances.
 * It contains a list of backhaul nodes. This list is terminated by a -1.
 *
 * CTSP_SET_SECTION :
 * This section is used for specifying CTSP instances.
 * Each entry has the following format:
 * c v1 v2 ... vk -1, where c is the color number (colors are numbered
 * from 1 to SALESMEN), and v1 v2 ... vk are vertices with color c
 * (vertices are numbered from 1 to Alg->Dimension).
 *
 * The first integer gives the number of the node. The last two integers
 * give the earliest and latest time for the node.
 *
 * DEMAND_SECTION :
 * The demands of all nodes of a CVRP are given in the form (per line)
 *
 *    <integer> <integer>
 *
 * The first integer spcifies a node number, the second its demand. The depot
 * nodes must also occcur in this section. Their demands are 0.
 *
 * DEPOT_SECTION :
 * Contains a list of possible alternate depot nodes. This list is terminated
 * by a -1. The current implementation allows only one depot.
 *
 * DRAFT_LIMIT_SECTION :
 * The draft limits of all nodes of a CVRP are give in the form (per line)
 *
 *    <integer> <integer>
 *
 * The first integer spcifies a node number, the second its draft limit.
 * The depot nodes must also occcur in this section. Their demands are 0.
 *
 * PICKUP_AND_DELIVERY_SECTION :
 * This section is used for specifying specifying pickup-and-delivery
 * instances. Each line is of the form
 *
 *     <integer> <integer> <Real> <Real> <Real> <integer> <integer>
 *
 * The first integer gives the number of the node.
 * The second integer gives its demand (ignored for PDTSPF, PDTSPL, VRPMPD
 * and VRPSPD instances).
 * The third and fourth number give the earliest and latest time for the node.
 * The fifth number specifies the service time for the node.
 * The last two integers are used to specify pickup and delivery. For a PDPTW,
 * PDTSP, PDTSPF and PDTSPL instance, the first of these integers gives the
 * index of the pickup sibling, whereas the second integer gives the index of
 * the delivery sibling. For a VRPMPD and VRPSPD instance, the two integers
 * simply give the size of the pickup and delivery for the node.
 *
 * SERVICE_TIME_SECTION :
 * The service times of all nodes of a CVRP are given in the form (per line)
 *
 *    <integer> <Real>
 *
 * The integer specifies a node number, the Real its service time.
 * The depot node must also occur in this section. Its service time is 0.
 *
 * TIME_WINDOW_SECTION :
 * Time windows are given in this section. Each line is of the form
 *
 *      <integer> <Real> <realr>
 *
 * The first integer specifies a node number. The two reals specify
 * earliest and latest arrival time for the node, respectively.
 */
namespace LKH {
static const thread_local char Delimiters[] = " :=\n\t\r\f\v\xef\xbb\xbf";
static void CheckSpecificationPart(LKHAlg *Alg);
static char *Copy(const char *S);

static void CreateNodes(LKHAlg *Alg);
static int FixEdge(LKHAlg::Node * Na, LKHAlg::Node * Nb);
static void Read_BACKHAUL_SECTION(LKHAlg *Alg);
static void Read_CAPACITY(LKHAlg *Alg);
static void Read_CTSP_SET_SECTION(LKHAlg *Alg);
static void Read_DEMAND_DIMENSION(LKHAlg *Alg);
static void Read_DEMAND_SECTION(LKHAlg *Alg);
static void Read_DEPOT_SECTION(LKHAlg *Alg);
static void Read_DIMENSION(LKHAlg *Alg);
static void Read_DISPLAY_DATA_SECTION(LKHAlg *Alg);
static void Read_DISPLAY_DATA_TYPE(LKHAlg *Alg);
static void Read_DISTANCE(LKHAlg *Alg);
static void Read_DRAFT_LIMIT_SECTION(LKHAlg *Alg);
static void Read_EDGE_DATA_FORMAT(LKHAlg *Alg);
static void Read_EDGE_DATA_SECTION(LKHAlg *Alg);
static void Read_EDGE_WEIGHT_FORMAT(LKHAlg *Alg);
static void Read_EDGE_WEIGHT_SECTION(LKHAlg *Alg);
static void Read_EDGE_WEIGHT_TYPE(LKHAlg *Alg);
static void Read_FIXED_EDGES_SECTION(LKHAlg *Alg);
static void Read_GRID_SIZE(LKHAlg *Alg);
static void Read_NAME(LKHAlg *Alg);
static void Read_NODE_COORD_SECTION(LKHAlg *Alg);
static void Read_NODE_COORD_TYPE(LKHAlg *Alg);
static void Read_PICKUP_AND_DELIVERY_SECTION(LKHAlg *Alg);
static void Read_RISK_THRESHOLD(LKHAlg *Alg);
static void Read_SALESMEN(LKHAlg *Alg);
static void Read_SCALE(LKHAlg *Alg);
static void Read_SERVICE_TIME(LKHAlg *Alg);
static void Read_SERVICE_TIME_SECTION(LKHAlg *Alg);
static void Read_TIME_WINDOW_SECTION(LKHAlg *Alg);
static void Read_TOUR_SECTION(FILE ** File, LKHAlg *Alg);
static void Read_TYPE(LKHAlg *Alg);
static int TwoDWeightType(LKHAlg *Alg);
static int ThreeDWeightType(LKHAlg *Alg);
static void Convert2FullMatrix(LKHAlg *Alg);

void LKHAlg::ReadProblem(const std::string& filepath)
{
    int i, j, K;
    char *Line, *Keyword;

	std::ostringstream FileName;
	FileName << filepath;
 
    
    if (!(ProblemFile = fopen(FileName.str().c_str(), "r")))
        eprintf("Cannot open PROBLEM_FILE: \"%s\"", ProblemFileName);
    if (TraceLevel >= 1)
        //printff("Reading PROBLEM_FILE: \"%s\" ... ", ProblemFileName);
    FreeStructures();
    FirstNode = 0;
    WeightType = WeightFormat = ProblemType = -1;
    CoordType = NO_COORDS;
    Name = Copy("Unnamed");
    Type = EdgeWeightType = EdgeWeightFormat = 0;
    EdgeDataFormat = NodeCoordType = DisplayDataType = 0;
    Distance = 0;
    C = 0;
    c = 0;
    DistanceLimit = std::numeric_limits<double>::max();
    m_filestream_p = 0;
    while ((Line = ReadLine(ProblemFile))) {

        if (!(Keyword = lkh_strtok(Line, Delimiters, &m_filestream_p)))
            continue;
        for (i = 0; i < (int)strlen(Keyword); i++) {
            Keyword[i] = (char)toupper(Keyword[i]);
        }

 
        if (!strcmp(Keyword, "COMMENT"));
        else if (!strcmp(Keyword, "BACKHAUL_SECTION"))
            Read_BACKHAUL_SECTION(this);
        else if (!strcmp(Keyword, "CAPACITY"))
            Read_CAPACITY(this);
        else if (!strcmp(Keyword, "CTSP_SET_SECTION"))
            Read_CTSP_SET_SECTION(this);
        else if (!strcmp(Keyword, "DEMAND_DIMENSION"))
            Read_DEMAND_DIMENSION(this);
        else if (!strcmp(Keyword, "DEMAND_SECTION"))
            Read_DEMAND_SECTION(this);
        else if (!strcmp(Keyword, "DEPOT_SECTION"))
            Read_DEPOT_SECTION(this);
        else if (!strcmp(Keyword, "DIMENSION"))
            Read_DIMENSION(this);
        else if (!strcmp(Keyword, "DISPLAY_DATA_SECTION"))
            Read_DISPLAY_DATA_SECTION(this);
        else if (!strcmp(Keyword, "DISPLAY_DATA_TYPE"))
            Read_DISPLAY_DATA_TYPE(this);
        else if (!strcmp(Keyword, "DISTANCE"))
            Read_DISTANCE(this);
        else if (!strcmp(Keyword, "DRAFT_LIMIT_SECTION"))
            Read_DRAFT_LIMIT_SECTION(this);
        else if (!strcmp(Keyword, "EDGE_DATA_FORMAT"))
            Read_EDGE_DATA_FORMAT(this);
        else if (!strcmp(Keyword, "EDGE_DATA_SECTION"))
            Read_EDGE_DATA_SECTION(this);
        else if (!strcmp(Keyword, "EDGE_WEIGHT_FORMAT"))
            Read_EDGE_WEIGHT_FORMAT(this);
        else if (!strcmp(Keyword, "EDGE_WEIGHT_SECTION"))
            Read_EDGE_WEIGHT_SECTION(this);
        else if (!strcmp(Keyword, "EDGE_WEIGHT_TYPE"))
            Read_EDGE_WEIGHT_TYPE(this);
        else if (!strcmp(Keyword, "EOF"))
            break;
        else if (!strcmp(Keyword, "FIXED_EDGES_SECTION"))
            Read_FIXED_EDGES_SECTION(this);
        else if (!strcmp(Keyword, "GRID_SIZE"))
            Read_GRID_SIZE(this);
        else if (!strcmp(Keyword, "NAME"))
            Read_NAME(this);
        else if (!strcmp(Keyword, "NODE_COORD_SECTION"))
            Read_NODE_COORD_SECTION(this);
        else if (!strcmp(Keyword, "NODE_COORD_TYPE"))
            Read_NODE_COORD_TYPE(this);
        else if (!strcmp(Keyword, "PICKUP_AND_DELIVERY_SECTION"))
            Read_PICKUP_AND_DELIVERY_SECTION(this);
        else if (!strcmp(Keyword, "RISK_THRESHOLD"))
            Read_RISK_THRESHOLD(this);
        else if (!strcmp(Keyword, "SALESMEN") ||
                 !strcmp(Keyword, "VEHICLES"))
            Read_SALESMEN(this);
        else if (!strcmp(Keyword, "SCALE"))
            Read_SCALE(this);
        else if (!strcmp(Keyword, "SERVICE_TIME"))
            Read_SERVICE_TIME(this);
        else if (!strcmp(Keyword, "SERVICE_TIME_SECTION"))
            Read_SERVICE_TIME_SECTION(this);
        else if (!strcmp(Keyword, "TIME_WINDOW_SECTION"))
            Read_TIME_WINDOW_SECTION(this);
        else if (!strcmp(Keyword, "TOUR_SECTION"))
            Read_TOUR_SECTION(&ProblemFile, this);
        else if (!strcmp(Keyword, "TYPE"))
            Read_TYPE(this);
        else
            eprintf("Unknown keyword: %s", Keyword);
    }
    Swaps = 0;

    /* Adjust parameters */
    if (Seed == 0)
        Seed = (unsigned) time(0);
    if (Precision == 0)
        Precision = 100;
    if (InitialStepSize == 0)
        InitialStepSize = 1;
    if (MaxSwaps < 0)
        MaxSwaps = Dimension;
    if (KickType > Dimension / 2)
        KickType = Dimension / 2;
    if (Runs == 0)
        Runs = 10;
    if (MaxCandidates > Dimension - 1)
        MaxCandidates = Dimension - 1;
    if (ExtraCandidates > Dimension - 1)
        ExtraCandidates = Dimension - 1;
    if (Scale < 1)
        Scale = 1;
    if (SubproblemSize >= Dimension)
        SubproblemSize = Dimension;
    else if (SubproblemSize == 0) {
        if (AscentCandidates > Dimension - 1)
            AscentCandidates = Dimension - 1;
        if (InitialPeriod < 0) {
            InitialPeriod = Dimension / 2;
            if (InitialPeriod < 100)
                InitialPeriod = 100;
        }
        if (Excess < 0)
            Excess = 1.0 / DimensionSaved * Salesmen;
        if (MaxTrials == -1)
            MaxTrials = Dimension;
        HeapMake(Dimension, this);
    }
    if (POPMUSIC_MaxNeighbors > Dimension - 1)
        POPMUSIC_MaxNeighbors = Dimension - 1;
    if (POPMUSIC_SampleSize > Dimension)
        POPMUSIC_SampleSize = Dimension;
    Depot = &NodeSet[MTSPDepot];
    if (ProblemType == CVRP) {
        Node *N;
        int MinSalesmen;
        if (Capacity <= 0)
            eprintf("CAPACITY not specified");
        TotalDemand = 0;
        N = FirstNode;
        do
            TotalDemand += N->Demand;
        while ((N = N->Suc) != FirstNode);
        MinSalesmen =
            TotalDemand / Capacity + (TotalDemand % Capacity != 0);
        if (Salesmen == 1) {
            Salesmen = MinSalesmen;
            if (Salesmen > Dimension)
                eprintf("CVRP: SALESMEN larger than DIMENSION");
        } else if (Salesmen < MinSalesmen)
            eprintf("CVRP: SALESMEN too small to meet demand");
        assert(Salesmen >= 1 && Salesmen <= Dimension);
        if (Salesmen == 1)
            ProblemType = TSP;
        Penalty = &LKHAlg::Penalty_CVRP;////
    } else if (ProblemType == SOP || ProblemType == M1_PDTSP) {
        Constraint *Con;
        Node *Ni, *Nj;
        int n, k;
        OldDistance = Distance;
        Distance = &LKHAlg::Distance_SOP;
        if (ProblemType == M1_PDTSP) {
            for (i = 2; i < Dim; i++) {
                Ni = &NodeSet[i];
                for (k = n = 0; k < DemandDimension; k++) {
                    n = Ni->M_Demand[k];
                    if (n >= 0)
                        continue;
                    for (j = 2; j < Dim; j++) {
                        if (j == i)
                            continue;
                        Nj = &NodeSet[j];
                        if (Nj->M_Demand[k] == -n) {
                            Ni->C[j] = -1;
                            break;
                        }
                    }
                }
            }
        }
        for (j = 2; j < Dim; j++) {
            Nj = &NodeSet[j];
            for (i = 2; i < Dim; i++) {
                if (i != j && Nj->C[i] == -1) {
                    Ni = &NodeSet[i];
                    assert(Con =
                           (Constraint *) malloc(sizeof(Constraint)));
                    Con->t1 = Ni;
                    Con->t2 = Nj;
                    Con->Suc = FirstConstraint;
                    FirstConstraint = Con;
                    Con->Next = Ni->FirstConstraint;
                    Ni->FirstConstraint = Con;
                }
            }
        }
        Salesmen = 1;
        Penalty = ProblemType == SOP ? &LKHAlg::Penalty_SOP : &LKHAlg::Penalty_M1_PDTSP;////
    }
    if (ProblemType == TSPTW) {
        Salesmen = 1;
        Penalty = &LKHAlg::Penalty_TSPTW;
    } else
        TSPTW_Makespan = 0;
    if (Salesmen > 1) {
        if (Salesmen > Dim && MTSPMinSize > 0)
            eprintf("Too many salesmen/vehicles (>= DIMENSION)");
        MTSP2TSP();
    }
    if (ProblemType == ACVRP)
        Penalty = &LKHAlg::Penalty_ACVRP;
    else if (ProblemType == CCVRP)
        Penalty = &LKHAlg::Penalty_CCVRP;
    else if (ProblemType == CTSP)
        Penalty = &LKHAlg::Penalty_CTSP;
    else if (ProblemType == CVRPTW)
        Penalty = &LKHAlg::Penalty_CVRPTW;
    else if (ProblemType == MLP)
        Penalty = &LKHAlg::Penalty_MLP;
    else if (ProblemType == OVRP)
        Penalty = &LKHAlg::Penalty_OVRP;
    else if (ProblemType == PDTSP)
        Penalty = &LKHAlg::Penalty_PDTSP;
    else if (ProblemType == PDTSPF)
        Penalty = &LKHAlg::Penalty_PDTSPF;
    else if (ProblemType == PDTSPL)
        Penalty = &LKHAlg::Penalty_PDTSPL;
    else if (ProblemType == PDPTW)
        Penalty = &LKHAlg::Penalty_PDPTW;
    else if (ProblemType == ONE_PDTSP)
        Penalty = &LKHAlg::Penalty_1_PDTSP;
    else if (ProblemType == M_PDTSP)
        Penalty = &LKHAlg::Penalty_M_PDTSP;
    else if (ProblemType == M1_PDTSP)
        Penalty = &LKHAlg::Penalty_M1_PDTSP;
    else if (ProblemType == RCTVRP || ProblemType == RCTVRPTW)
        Penalty = &LKHAlg::Penalty_RCTVRP;
    else if (ProblemType == TRP)
        Penalty = &LKHAlg::Penalty_TRP;
    else if (ProblemType == TSPDL)
        Penalty = &LKHAlg::Penalty_TSPDL;
    else if (ProblemType == TSPPD)
        Penalty = &LKHAlg::Penalty_TSPPD;
    if (ProblemType == VRPB)
        Penalty = &LKHAlg::Penalty_VRPB;
    else if (ProblemType == VRPBTW)
        Penalty = &LKHAlg::Penalty_VRPBTW;
    else if (ProblemType == VRPPD)
        Penalty = &LKHAlg::Penalty_VRPPD;
    if (BWTSP_B > 0) {
        if (Penalty)
            eprintf("BWTSP not compatible with problem type %s\n", Type);
        ProblemType = BWTSP;
        free(Type);
        Type = Copy("BWTSP");
        Penalty = &LKHAlg::Penalty_BWTSP;
        if (BWTSP_L != std::numeric_limits<int>::max())
            BWTSP_L *= Scale;
    }
    if (Penalty && (SubproblemSize > 0 || SubproblemTourFile))
        eprintf("Partitioning not implemented for constrained problems");
    Depot->DepotId = 1;
    for (i = Dim + 1; i <= DimensionSaved; i++)
        NodeSet[i].DepotId = i - Dim + 1;
    if (Dimension != DimensionSaved) {
        NodeSet[Depot->Id + DimensionSaved].DepotId = 1;
        for (i = Dim + 1; i <= DimensionSaved; i++)
            NodeSet[i + DimensionSaved].DepotId = i - Dim + 1;
    }
    if (Scale < 1)
        Scale = 1;
    else {
        Node *Ni = FirstNode;
        do {
            Ni->Earliest *= Scale;
            Ni->Latest *= Scale;
            Ni->ServiceTime *= Scale;
        } while ((Ni = Ni->Suc) != FirstNode);
        ServiceTime *= Scale;
        RiskThreshold *= Scale;
        if (DistanceLimit != std::numeric_limits<double>::max())
            DistanceLimit *= Scale;
    }
    if (ServiceTime != 0) {
        for (i = 1; i <= Dim; i++)
            NodeSet[i].ServiceTime = ServiceTime;
        Depot->ServiceTime = 0;
    }
    if (CostMatrix == 0 && Dimension <= MaxMatrixDimension &&
        Distance != 0 && Distance != &LKHAlg::Distance_1
        && Distance != &LKHAlg::Distance_LARGE && Distance != &LKHAlg::Distance_ATSP
        && Distance != &LKHAlg::Distance_MTSP && Distance != &LKHAlg::Distance_SPECIAL) {
        Node *Ni, *Nj;
        assert(CostMatrix =
               (int *) calloc((size_t) Dim * (Dim - 1) / 2, sizeof(int)));
        Ni = FirstNode->Suc;
        do {
            Ni->C =
                &CostMatrix[(size_t) (Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
            if (ProblemType != HPP || Ni->Id <= Dim)
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = Fixed(Ni, Nj) ? 0 : (this->*Distance)(Ni, Nj);
            else
                for (Nj = FirstNode; Nj != Ni; Nj = Nj->Suc)
                    Ni->C[Nj->Id] = 0;
        }
        while ((Ni = Ni->Suc) != FirstNode);
        c = 0;
        WeightType = EXPLICIT;
    }
    if (ProblemType == TSPTW ||
        ProblemType == CVRPTW || ProblemType == VRPBTW ||
        ProblemType == PDPTW || ProblemType == RCTVRPTW) {
        M = std::numeric_limits<int>::max() / 2 / Precision;
        for (i = 1; i <= Dim; i++) {
            Node *Ni = &NodeSet[i];
            for (j = 1; j <= Dim; j++) {
                Node *Nj = &NodeSet[j];
                if (Ni != Nj &&
                    Ni->Earliest + Ni->ServiceTime + Ni->C[j] > Nj->Latest)
                    Ni->C[j] = M;
            }
        }
    }
    C = WeightType == EXPLICIT ? &LKHAlg::C_EXPLICIT : &LKHAlg::C_FUNCTION;
    D = WeightType == EXPLICIT ? &LKHAlg::D_EXPLICIT : &LKHAlg::D_FUNCTION;
    if (ProblemType != CVRP && ProblemType != CVRPTW &&
        ProblemType != CTSP &&
        ProblemType != TSP && ProblemType != ATSP) {
        M = std::numeric_limits<int>::max() / 2 / Precision;
        for (i = Dim + 1; i <= DimensionSaved; i++) {
            for (j = 1; j <= DimensionSaved; j++) {
                if (j == i)
                    continue;
                if (j == MTSPDepot || j > Dim)
                    NodeSet[i].C[j] = NodeSet[MTSPDepot].C[j] = M;
                NodeSet[i].C[j] = NodeSet[MTSPDepot].C[j];
                NodeSet[j].C[i] = NodeSet[j].C[MTSPDepot];
            }
        }
        if (ProblemType == CCVRP || ProblemType == OVRP)
            for (i = 1; i <= Dim; i++)
                NodeSet[i].C[MTSPDepot] = 0;
    }
    if (Precision > 1 && CostMatrix) {
        for (i = 2; i <= Dim; i++) {
            Node *N = &NodeSet[i];
            for (j = 1; j < i; j++)
                if (N->C[j] * Precision / Precision != N->C[j])
                    eprintf("PRECISION (= %d) is too large", Precision);
        }
    }
    if (SubsequentMoveType == 0) {
        SubsequentMoveType = MoveType;
        SubsequentMoveTypeSpecial = MoveTypeSpecial;
    }
    K = MoveType >= SubsequentMoveType || !SubsequentPatching ?
        MoveType : SubsequentMoveType;
    if (PatchingC > K)
        PatchingC = K;
    if (PatchingA > 1 && PatchingA >= PatchingC)
        PatchingA = PatchingC > 2 ? PatchingC - 1 : 1;
    if (NonsequentialMoveType == -1 ||
        NonsequentialMoveType > K + PatchingC + PatchingA - 1)
        NonsequentialMoveType = K + PatchingC + PatchingA - 1;
    if (PatchingC >= 1) {
        BestMove = BestSubsequentMove = &LKHAlg::BestKOptMove;
        if (!SubsequentPatching && SubsequentMoveType <= 5) {
            MoveFunction BestOptMove[] =
                { 0, 0, &LKHAlg::Best2OptMove, &LKHAlg::Best3OptMove,
				&LKHAlg::Best4OptMove, &LKHAlg::Best5OptMove
            };
            BestSubsequentMove = BestOptMove[SubsequentMoveType];
        }
    } else {
        MoveFunction BestOptMove[] = { 0, 0, &LKHAlg::Best2OptMove, &LKHAlg::Best3OptMove,
			&LKHAlg::Best4OptMove, &LKHAlg::Best5OptMove
        };
        BestMove = MoveType <= 5 ? BestOptMove[MoveType] : &LKHAlg::BestKOptMove;
        BestSubsequentMove = SubsequentMoveType <= 5 ?
            BestOptMove[SubsequentMoveType] : &LKHAlg::BestKOptMove;
    }
    if (MoveTypeSpecial)
        BestMove = &LKHAlg::BestSpecialOptMove;
    if (SubsequentMoveTypeSpecial)
        BestSubsequentMove = &LKHAlg::BestSpecialOptMove;
    if (ProblemType == HCP || ProblemType == HPP)
        MaxCandidates = 0;
    /*if (TraceLevel >= 1) {
        printff("done\n");
        PrintParameters();
    } else
        printff("PROBLEM_FILE = %s\n",
                ProblemFileName ? ProblemFileName : "");*/
	fclose(ProblemFile);
    if (InitialTourFileName)
        ReadTour(InitialTourFileName, &InitialTourFile);
    if (InputTourFileName)
        ReadTour(InputTourFileName, &InputTourFile);
    if (SubproblemTourFileName && SubproblemSize > 0)
        ReadTour(SubproblemTourFileName, &SubproblemTourFile);
    if (MergeTourFiles >= 1) {
        free(MergeTourFile);
        assert(MergeTourFile =
               (FILE **) malloc(MergeTourFiles * sizeof(FILE *)));
        for (i = 0; i < MergeTourFiles; i++)
            ReadTour(MergeTourFileName[i], &MergeTourFile[i]);
    }
    free(LastLine);
    LastLine = 0;
}

static int TwoDWeightType(LKHAlg *Alg)
{
    if (Alg->Asymmetric)
        return 0;
    return Alg->WeightType == LKH::EUC_2D || Alg->WeightType == LKH::MAX_2D ||
        Alg->WeightType == LKH::MAN_2D || Alg->WeightType == LKH::CEIL_2D ||
        Alg->WeightType == LKH::FLOOR_2D ||
        Alg->WeightType == LKH::GEO || Alg->WeightType == LKH::GEOM ||
        Alg->WeightType == LKH::GEO_MEEUS || Alg->WeightType == LKH::GEOM_MEEUS ||
        Alg->WeightType == LKH::ATT || Alg->WeightType == LKH::TOR_2D ||
        (Alg->WeightType == LKH::SPECIAL && Alg->CoordType == LKH::TWOD_COORDS);
}

static int ThreeDWeightType(LKHAlg *Alg)
{
    if (Alg->Asymmetric)
        return 0;
    return Alg->WeightType == LKH::EUC_3D || Alg->WeightType == LKH::MAX_3D ||
        Alg->WeightType == LKH::MAN_3D || Alg->WeightType == LKH::CEIL_3D ||
        Alg->WeightType == LKH::FLOOR_3D || Alg->WeightType == LKH::TOR_3D ||
        Alg->WeightType == LKH::XRAY1 || Alg->WeightType == LKH::XRAY2 ||
        (Alg->WeightType == LKH::SPECIAL && Alg->CoordType == LKH::THREED_COORDS);
}

static void CheckSpecificationPart(LKHAlg *Alg)
{
    if (Alg->ProblemType == -1)
		Alg->eprintf("TYPE is missing");
    if (Alg->Dimension < 3)
		Alg->eprintf("DIMENSION < 3 or not specified");
    if (Alg->WeightType == -1 && !Alg->Asymmetric && Alg->ProblemType != LKH::HCP &&
		Alg->ProblemType != LKH::HPP && !Alg->EdgeWeightType)
		Alg->eprintf("EDGE_WEIGHT_TYPE is missing");
    if (Alg->WeightType == LKH::EXPLICIT && Alg->WeightFormat == -1 && !Alg->EdgeWeightFormat)
		Alg->eprintf("EDGE_WEIGHT_FORMAT is missing");
    if (Alg->WeightType == LKH::EXPLICIT && Alg->WeightFormat == LKH::FUNCTION)
		Alg->eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if (Alg->WeightType != LKH::EXPLICIT &&
        (Alg->WeightType != LKH::SPECIAL || Alg->CoordType != LKH::NO_COORDS) &&
        Alg->WeightType != -1 && Alg->WeightFormat != -1 && Alg->WeightFormat != LKH::FUNCTION)
		Alg->eprintf("Conflicting EDGE_WEIGHT_TYPE and EDGE_WEIGHT_FORMAT");
    if ((Alg->ProblemType == LKH::ATSP || Alg->ProblemType == LKH::SOP) &&
        Alg->WeightType != LKH::EXPLICIT && Alg->WeightType != -1)
		Alg->eprintf("Conflicting TYPE and EDGE_WEIGHT_TYPE");
    if (Alg->CandidateSetType == LKH::DELAUNAY && !TwoDWeightType(Alg) &&
        Alg->MaxCandidates > 0)
		Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = DELAUNAY");
    if (Alg->CandidateSetType == LKH::QUADRANT && !TwoDWeightType(Alg) &&
        !ThreeDWeightType(Alg) && Alg->MaxCandidates + Alg->ExtraCandidates > 0)
		Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for CANDIDATE_SET_TYPE = QUADRANT");
    if (Alg->ExtraCandidateSetType == LKH::QUADRANT && !TwoDWeightType(Alg) &&
        !ThreeDWeightType(Alg) && Alg->ExtraCandidates > 0)
		Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for EXTRA_CANDIDATE_SET_TYPE = "
             "QUADRANT");
    if (Alg->InitialTourAlgorithm == LKH::QUICK_BORUVKA && !TwoDWeightType(Alg) &&
        !ThreeDWeightType(Alg))
		Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "QUICK-BORUVKA");
    if (Alg->InitialTourAlgorithm == LKH::SIERPINSKI && !TwoDWeightType(Alg))
		Alg->eprintf
            ("Illegal EDGE_WEIGHT_TYPE for INITIAL_TOUR_ALGORITHM = "
             "SIERPINSKI");
    if (Alg->DelaunayPartitioning && !TwoDWeightType(Alg))
		Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for DELAUNAY specification");
    if (Alg->KarpPartitioning && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
		Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for KARP specification");
    if (Alg->KCenterPartitioning && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
		Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for K-CENTER specification");
    if (Alg->KMeansPartitioning && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
		Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for K-MEANS specification");
    if (Alg->MoorePartitioning && !TwoDWeightType(Alg))
		Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for MOORE specification");
    if (Alg->RohePartitioning && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
		Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for ROHE specification");
    if (Alg->SierpinskiPartitioning && !TwoDWeightType(Alg))
		Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for SIERPINSKI specification");
    if (Alg->SubproblemBorders && !TwoDWeightType(Alg) && !ThreeDWeightType(Alg))
		Alg->eprintf("Illegal EDGE_WEIGHT_TYPE for BORDERS specification");
    if (Alg->InitialTourAlgorithm == LKH::MTSP_ALG && Alg->Asymmetric)
		Alg->eprintf("INTIAL_TOUR_ALGORITHM = MTSP is not applicable for "
                "asymetric problems");
}

static char *Copy(const char *S)
{
    char *Buffer;

    if (!S || strlen(S) == 0)
        return 0;
    assert(Buffer = (char *) malloc(strlen(S) + 1));
    strcpy(Buffer, S);
    return Buffer;
}

static void CreateNodes(LKHAlg *Alg)
{
    LKHAlg::Node *Prev = 0, *N = 0;
    int i;

    if (Alg->Dimension <= 0)
		Alg->eprintf("DIMENSION is not positive (or not specified)");
    if (Alg->Asymmetric) {
		Alg->Dim = Alg->DimensionSaved;
		Alg->DimensionSaved = Alg->Dimension + Alg->Salesmen - 1;
		Alg->Dimension = 2 * Alg->DimensionSaved;
    } else if (Alg->ProblemType == LKH::HPP) {
		Alg->Dimension++;
        if (Alg->Dimension > Alg->MaxMatrixDimension)
			Alg->eprintf("DIMENSION too large in HPP problem");
    }
    assert(Alg->NodeSet = (LKHAlg::Node *) calloc(Alg->Dimension + 1, sizeof(LKHAlg::Node)));
    for (i = 1; i <= Alg->Dimension; i++, Prev = N) {
        N = &Alg->NodeSet[i];
        if (i == 1)
			Alg->FirstNode = N;
        else
            LKH::LKHAlg::Link(Prev, N);
        N->Id = i;
        if (Alg->MergeTourFiles >= 1)
            assert(N->MergeSuc =
                   (LKHAlg::Node **) calloc(Alg->MergeTourFiles, sizeof(LKHAlg::Node *)));
        N->Earliest = 0;
        N->Latest = std::numeric_limits<int>::max();
    }
    LKH::LKHAlg::Link(N, Alg->FirstNode);
}


static int FixEdge(LKHAlg::Node * Na, LKHAlg::Node * Nb)
{
    if (!Na->FixedTo1 || Na->FixedTo1 == Nb)
        Na->FixedTo1 = Nb;
    else if (!Na->FixedTo2 || Na->FixedTo2 == Nb)
        Na->FixedTo2 = Nb;
    else
        return 0;
    if (!Nb->FixedTo1 || Nb->FixedTo1 == Na)
        Nb->FixedTo1 = Na;
    else if (!Nb->FixedTo2 || Nb->FixedTo1 == Na)
        Nb->FixedTo2 = Na;
    else
        return 0;
    return 1;
}

static void Read_NAME(LKHAlg *Alg)
{
    if (!(Alg->Name = Copy(Alg->lkh_strtok(0, Delimiters, &Alg->m_filestream_p))))
		Alg->eprintf("NAME: string expected");
}

static void Read_BACKHAUL_SECTION(LKHAlg *Alg)
{
    int Id;

    while (Alg->fscanint(Alg->ProblemFile, &Id) && Id != -1) {
        if (Id <= 0 || Id > Alg->Dim)
			Alg->eprintf("BACKHAUL_SECTION: Node number out of range: %d", Id);
		Alg->NodeSet[Id].Backhaul = 1;
		Alg->NodeSet[Id + Alg->DimensionSaved].Backhaul = 1;
    }
}

static void Read_CAPACITY(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters, &Alg->m_filestream_p);

    if (!Token || !sscanf(Token, "%d", &Alg->Capacity))
		Alg->eprintf("CAPACITY: Integer expected");
}

static void Read_CTSP_SET_SECTION(LKHAlg *Alg)
{   
    LKHAlg::Node *N;
    int Id, n, *ColorUsed;
    
    N = Alg->FirstNode;
    do {
        N->Color = 0;
    } while ((N = N->Suc) != Alg->FirstNode);
    assert(ColorUsed = (int *) calloc(Alg->Salesmen + 1, sizeof(int)));
    while (fscanf(Alg->ProblemFile, "%d", &Id) > 0) {
        if (Id < 1 || Id > Alg->Salesmen)
			Alg->eprintf("(CTSP_SET_SECTION) Color number %d outside range", Id);
        if (ColorUsed[Id])
			Alg->eprintf("(CTSP_SET_SECTION) Color number %d used twice", Id);
        ColorUsed[Id] = 1;
        for (;;) {
            if (fscanf(Alg->ProblemFile, "%d", &n) != 1)
				Alg->eprintf("(CTSP_SET_SECTION) Missing -1");
            if (n == -1)
                break; 
             if (n < 1 || n > Alg->DimensionSaved)
				 Alg->eprintf("(CTSP_SET_SECTION) Node %d outside range", n);
             N = &Alg->NodeSet[n];
             if (N->Color != 0 && N->Color != Id) 
				 Alg->eprintf("(CTSP_SET_SECTION) Node %d in two sets", n);
             if (N == Alg->Depot)
				 Alg->eprintf("(CTSP_SET_SECTION) Depot %d occurs in set %d", n, Id);
             N->Color = Id;
        }
    }
    free(ColorUsed);
}

static void Read_DEMAND_DIMENSION(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p);

    if (!Token || !sscanf(Token, "%d", &(Alg->DemandDimension)))
		Alg->eprintf("DIMENSION_DIMENSION: Integer expected");
    if (Alg->DemandDimension < 0)
		Alg->eprintf("DIMENSION_DIMENSION: < 0");
}

static void Read_DEMAND_SECTION(LKHAlg *Alg)
{
    int Id, Demand, i, k;
    LKHAlg::Node *N;

    for (i = 1; i <= Alg->Dim; i++) {
        Alg->fscanint(Alg->ProblemFile, &Id);
        if (Id <= 0 || Id > Alg->Dim)
			Alg->eprintf("DEMAND_SECTION: Node number out of range: %d", Id);
        N = &Alg->NodeSet[Id];
        if (Alg->DemandDimension > 1) {
            assert(N->M_Demand =
                   (int *) malloc(Alg->DemandDimension * sizeof(int)));
            for (k = 0; k < Alg->DemandDimension; k++) {
                if (!Alg->fscanint(Alg->ProblemFile, &Demand))
					Alg->eprintf("DEMAND_SECTION: Missing demand for node %d",
                            Id);
                N->M_Demand[k] = Demand;
            }
        } else if (!Alg->fscanint(Alg->ProblemFile, &N->Demand))
			Alg->eprintf("DEMAND_SECTION: Missing demand for node %d", Id);
    }
}

static void Read_DEPOT_SECTION(LKHAlg *Alg)
{
    int i;
    if (!Alg->fscanint(Alg->ProblemFile, &(Alg->MTSPDepot)))
		Alg->eprintf("DEPOT_SECTION: Integer expected");
    else if (Alg->MTSPDepot <= 0)
		Alg->eprintf("DEPOT_SECTION: Positive value expected");
    if (Alg->fscanint(Alg->ProblemFile, &i) && i != -1)
		Alg->eprintf("DEPOT_SECTION: Only one depot allowed");
}

static void Read_DIMENSION(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p);

    if (!Token || !sscanf(Token, "%d", &Alg->Dimension))
		Alg->eprintf("DIMENSION: Integer expected");
    if (Alg->Dimension < 0)
		Alg->eprintf("DIMENSION: < 0");
	Alg->DimensionSaved = Alg->Dim = Alg->Dimension;
}

static void Read_DISPLAY_DATA_SECTION(LKHAlg *Alg)
{
    LKHAlg::Node *N;
    int Id, i;

    CheckSpecificationPart(Alg);
    if (Alg->ProblemType == LKH::HPP)
		Alg->Dimension--;
    if (!Alg->DisplayDataType || strcmp(Alg->DisplayDataType, "TWOD_DISPLAY"))
		Alg->eprintf
            ("DISPLAY_DATA_SECTION conflicts with DISPLAY_DATA_TYPE: %s",
				Alg->DisplayDataType);
    if (!Alg->FirstNode)
        CreateNodes(Alg);
    N = Alg->FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != Alg->FirstNode);
    for (i = 1; i <= Alg->Dim; i++) {
        if (!Alg->fscanint(Alg->ProblemFile, &Id))
			Alg->eprintf("DIPLAY_DATA_SECTION: Missing nodes");
        if (Id <= 0 || Id > Alg->Dimension)
			Alg->eprintf("DIPLAY_DATA_SECTION: Node number out of range: %d",
                    Id);
        N = &Alg->NodeSet[Id];
        if (N->V == 1)
			Alg->eprintf("DIPLAY_DATA_SECTION: Node number occurs twice: %d",
                    N->Id);
        N->V = 1;
        if (!fscanf(Alg->ProblemFile, "%lf", &N->X))
			Alg->eprintf("DIPLAY_DATA_SECTION: Missing X-coordinate");
        if (!fscanf(Alg->ProblemFile, "%lf", &N->Y))
			Alg->eprintf("DIPLAY_DATA_SECTION: Missing Y-coordinate");
    }
    N = Alg->FirstNode;
    do
        if (!N->V && N->Id <= Alg->Dim)
            break;
    while ((N = N->Suc) != Alg->FirstNode);
    if (!N->V)
		Alg->eprintf("DIPLAY_DATA_SECTION: No coordinates given for node %d",
                N->Id);
    if (Alg->ProblemType == LKH::HPP)
		Alg->Dimension++;
}

static void Read_DISPLAY_DATA_TYPE(LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->DisplayDataType = Copy(Alg->lkh_strtok(0, Delimiters, &Alg->m_filestream_p))))
		Alg->eprintf("DISPLAY_DATA_TYPE: string expected");
    for (i = 0; i < strlen(Alg->DisplayDataType); i++)
		Alg->DisplayDataType[i] = (char) toupper(Alg->DisplayDataType[i]);
    if (strcmp(Alg->DisplayDataType, "COORD_DISPLAY") &&
        strcmp(Alg->DisplayDataType, "TWOD_DISPLAY") &&
        strcmp(Alg->DisplayDataType, "NO_DISPLAY"))
		Alg->eprintf("Unknown DISPLAY_DATA_TYPE: %s", Alg->DisplayDataType);
}

static void Read_DISTANCE(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p);

    if (!Token || !sscanf(Token, "%lf", &Alg->DistanceLimit))
		Alg->eprintf("DISTANCE: Real expected");
}

static void Read_DRAFT_LIMIT_SECTION(LKHAlg *Alg)
{
    int Id, i;
    LKHAlg::Node *N;

    for (i = 1; i <= Alg->Dim; i++) {
        Alg->fscanint(Alg->ProblemFile, &Id);
        if (Id <= 0 || Id > Alg->Dim)
			Alg->eprintf("DRAFT_LIMIT_SECTION: Node number out of range: %d",
                    Id);
        N = &Alg->NodeSet[Id];
        if (!Alg->fscanint(Alg->ProblemFile, &N->DraftLimit))
			Alg->eprintf("DRAFT_LIMIT_SECTION: Missing draft limit for node %d",
                    Id);
    }
}

static void Read_EDGE_DATA_FORMAT(LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->EdgeDataFormat = Copy(Alg->lkh_strtok(0, Delimiters, &Alg->m_filestream_p))))
		Alg->eprintf("EDGE_DATA_FORMAT: string expected");
    for (i = 0; i < strlen(Alg->EdgeDataFormat); i++)
		Alg->EdgeDataFormat[i] = (char) toupper(Alg->EdgeDataFormat[i]);
    if (strcmp(Alg->EdgeDataFormat, "EDGE_LIST") &&
        strcmp(Alg->EdgeDataFormat, "ADJ_LIST"))
		Alg->eprintf("Unknown EDGE_DATA_FORMAT: %s", Alg->EdgeDataFormat);
    if (Alg->SubproblemTourFileName)
		Alg->eprintf("EDGE_DATA_FORMAT"
                " cannot be used together with SUBPROBLEM_TOUR_FILE");
}

static void Read_EDGE_DATA_SECTION(LKHAlg *Alg)
{
    LKHAlg::Node *Ni, *Nj;
    int i, j, W = 0, WithWeights = 0, FirstLine = 1;
    double w = 0;
    char *Line;

    CheckSpecificationPart(Alg);
    if (!Alg->EdgeDataFormat)
        Alg->eprintf("Missing EDGE_DATA_FORMAT specification");
    if (!Alg->FirstNode)
        CreateNodes(Alg);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    if (!strcmp(Alg->EdgeDataFormat, "EDGE_LIST")) {
        Line = Alg->ReadLine(Alg->ProblemFile);
        if (sscanf(Line, "%d %d %lf\n", &i, &j, &w) == 3)
            WithWeights = 1;
        W = round(Alg->Scale * w);
        while (i != -1) {
            if (i <= 0 ||
                i > (!Alg->Asymmetric ? Alg->Dimension : Alg->Dimension / 2))
                Alg->eprintf("(EDGE_DATA_SECTION) Node number out of range: %d", i);
            if (!FirstLine)
                Alg->fscanint(Alg->ProblemFile, &j);
            if (j <= 0
                || j > (!Alg->Asymmetric ? Alg->Dimension : Alg->Dimension / 2))
                Alg->eprintf("(EDGE_DATA_SECTION) Node number out of range: %d",
                        j);
            if (i == j)
                Alg->eprintf("(EDGE_DATA_SECTION) Illegal edge: %d to %d",
                        i, j);
            if (Alg->Asymmetric)
                j += Alg->Dimension / 2;
            Ni = &Alg->NodeSet[i];
            Nj = &Alg->NodeSet[j];
            if (WithWeights) {
                if (!FirstLine) {
                    fscanf(Alg->ProblemFile, "%lf", &w);
                    W = round(Alg->Scale * w);
                }
                W *= Alg->Precision;
            }
            Alg->AddCandidate(Ni, Nj, W, 1);
            Alg->AddCandidate(Nj, Ni, W, 1);
            FirstLine = 0;
            if (!Alg->fscanint(Alg->ProblemFile, &i))
                i = -1;
        }
    } else if (!strcmp(Alg->EdgeDataFormat, "ADJ_LIST")) {
        if (!Alg->fscanint(Alg->ProblemFile, &i))
            i = -1;
        while (i != -1) {
            if (i <= 0 ||
                (!Alg->Asymmetric ? Alg->Dimension : Alg->Dimension / 2))
                Alg->eprintf
                ("(EDGE_DATA_SECTION) Node number out of range: %d",
                 i);
            Ni = &Alg->NodeSet[i];
            Alg->fscanint(Alg->ProblemFile, &j);
            while (j != -1) {
                if (j <= 0 ||
                    (Alg->ProblemType != LKH::ATSP ? Alg->Dimension : Alg->Dimension / 2))
                    Alg->eprintf
                    ("(EDGE_DATA_SECTION) Node number out of range: %d",
                     j);
                if (i == j)
                    Alg->eprintf("(EDGE_DATA_SECTION) Illgal edge: %d to %d",
                            i, j);
                if (Alg->Asymmetric)
                    j += Alg->Dimension / 2;
                Nj = &Alg->NodeSet[j];
                Alg->AddCandidate(Ni, Nj, 0, 1);
                Alg->AddCandidate(Nj, Ni, 0, 1);
                Alg->fscanint(Alg->ProblemFile, &j);
            }
            Alg->fscanint(Alg->ProblemFile, &i);
        }
    } else
        Alg->eprintf("(EDGE_DATA_SECTION) No EDGE_DATA_FORMAT specified");
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
    if (Alg->Asymmetric) {
        for (i = 1; i <= Alg->DimensionSaved; i++)
            FixEdge(&Alg->NodeSet[i], &Alg->NodeSet[i + Alg->DimensionSaved]);
    }
    Alg->WeightType = 1;
    Alg->MaxCandidates = Alg->ExtraCandidates = 0;
    Alg->Distance = WithWeights ? &LKHAlg::Distance_LARGE : &LKHAlg::Distance_1;
}

static void Read_EDGE_WEIGHT_FORMAT(LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->EdgeWeightFormat = Copy(Alg->lkh_strtok(0, Delimiters, &Alg->m_filestream_p))))
        Alg->eprintf("EDGE_WEIGHT_FORMAT: string expected");
    for (i = 0; i < strlen(Alg->EdgeWeightFormat); i++)
        Alg->EdgeWeightFormat[i] = (char) toupper(Alg->EdgeWeightFormat[i]);
    if (!strcmp(Alg->EdgeWeightFormat, "FUNCTION"))
        Alg->WeightFormat = LKH::FUNCTION;
    else if (!strcmp(Alg->EdgeWeightFormat, "FULL_MATRIX"))
        Alg->WeightFormat = LKH::FULL_MATRIX;
    else if (!strcmp(Alg->EdgeWeightFormat, "UPPER_ROW"))
        Alg->WeightFormat = LKH::UPPER_ROW;
    else if (!strcmp(Alg->EdgeWeightFormat, "LOWER_ROW"))
        Alg->WeightFormat = LKH::LOWER_ROW;
    else if (!strcmp(Alg->EdgeWeightFormat, "UPPER_DIAG_ROW"))
        Alg->WeightFormat = LKH::UPPER_DIAG_ROW;
    else if (!strcmp(Alg->EdgeWeightFormat, "LOWER_DIAG_ROW"))
        Alg->WeightFormat = LKH::LOWER_DIAG_ROW;
    else if (!strcmp(Alg->EdgeWeightFormat, "UPPER_COL"))
        Alg->WeightFormat = LKH::UPPER_COL;
    else if (!strcmp(Alg->EdgeWeightFormat, "LOWER_COL"))
        Alg->WeightFormat = LKH::LOWER_COL;
    else if (!strcmp(Alg->EdgeWeightFormat, "UPPER_DIAG_COL"))
        Alg->WeightFormat = LKH::UPPER_DIAG_COL;
    else if (!strcmp(Alg->EdgeWeightFormat, "LOWER_DIAG_COL"))
        Alg->WeightFormat = LKH::LOWER_DIAG_COL;
    else
        Alg->eprintf("Unknown EDGE_WEIGHT_FORMAT: %s", Alg->EdgeWeightFormat);
}

static void Read_EDGE_WEIGHT_SECTION(LKHAlg *Alg)
{
	LKHAlg::Node *Ni;
    int i, j, n, W;
    double w;

    if (Alg->ProblemType == LKH::SOP && Alg->ProblemType != LKH::M1_PDTSP) {
        Alg->fscanint(Alg->ProblemFile, &n);
        if (n != Alg->Dimension)
            Alg->eprintf("SOP: DIMENSION != n (%d != %d)", Alg->Dimension, n);
    } else
        n = Alg->Dimension;
    CheckSpecificationPart(Alg);
    if (!Alg->FirstNode)
        CreateNodes(Alg);
    if (!Alg->Asymmetric) {
        assert(Alg->CostMatrix =
               (int *) calloc((size_t) Alg->Dimension * (Alg->Dimension - 1) / 2,
                              sizeof(int)));
        Ni = Alg->FirstNode->Suc;
        do {
            Ni->C =
                &Alg->CostMatrix[(size_t) (Ni->Id - 1) * (Ni->Id - 2) / 2] - 1;
        }
        while ((Ni = Ni->Suc) != Alg->FirstNode);
    } else {
        n = Alg->Dimension / 2;
        assert(Alg->CostMatrix = (int *) calloc((size_t) n * n, sizeof(int)));
        for (Ni = Alg->FirstNode; Ni->Id <= n; Ni = Ni->Suc)
            Ni->C = &Alg->CostMatrix[(size_t) (Ni->Id - 1) * n] - 1;
    }
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    switch (Alg->WeightFormat) {
	case LKH::FULL_MATRIX:
        for (i = 1; i <= Alg->Dim; i++) {
            Ni = &Alg->NodeSet[i];
            for (j = 1; j <= Alg->Dim; j++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                W = round(Alg->Scale * w);
                if (W > std::numeric_limits<int>::max() / 2 / Alg->Precision)
                    W = std::numeric_limits<int>::max() / 2 / Alg->Precision;
                if (Alg->Salesmen > 1 && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                if (Alg->Asymmetric) {
                    Ni->C[j] = W;
                    if (j != i && W > Alg->M)
						Alg->M = W;
                } else if (j < i)
                    Ni->C[j] = W;
            }
        }
        break;
    case LKH::UPPER_ROW:
        for (i = 1; i < Alg->Dim; i++) {
            for (j = i + 1; j <= Alg->Dim; j++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                W = round(Alg->Scale * w);
                if (W > std::numeric_limits<int>::max() / 2 / Alg->Precision)
                    W = std::numeric_limits<int>::max() / 2 / Alg->Precision;
                if (Alg->Penalty && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                Alg->NodeSet[j].C[i] = W;
                if (Alg->Asymmetric) {
                    Alg->NodeSet[i].C[j] = W;
                    if (j != i && W > Alg->M)
                        Alg->M = W;
                }
            }
        }
        break;
    case LKH::LOWER_ROW:
        for (i = 2; i <= Alg->Dim; i++) {
            for (j = 1; j < i; j++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                W = round(Alg->Scale * w);
                if (Alg->Penalty && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                Alg->NodeSet[i].C[j] = W;
                if (Alg->Asymmetric) {
                    Alg->NodeSet[j].C[i] = W;
                    if (j != i && W > Alg->M)
                        Alg->M = W;
                }
            }
        }
        break;
	case LKH::UPPER_DIAG_ROW:
        for (i = 1; i <= Alg->Dim; i++) {
            for (j = i; j <= Alg->Dim; j++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                W = round(Alg->Scale * w);
                if (W > std::numeric_limits<int>::max() / 2 / Alg->Precision)
                    W = std::numeric_limits<int>::max() / 2 / Alg->Precision;
                if (Alg->Penalty && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                Alg->NodeSet[j].C[i] = W;
                if (Alg->Asymmetric) {
                    Alg->NodeSet[j].C[i] = W;
                    if (W > Alg->M)
                        Alg->M = W;
                }
            }
        }
        break;
    case LKH::LOWER_DIAG_ROW:
        for (i = 1; i <= Alg->Dim; i++) {
            for (j = 1; j <= i; j++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                W = round(Alg->Scale * w);
                if (W > std::numeric_limits<int>::max() / 2 / Alg->Precision)
                    W = std::numeric_limits<int>::max() / 2 / Alg->Precision;
                if (Alg->Penalty && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                if (j != i)
                    Alg->NodeSet[i].C[j] = W;
                if (Alg->Asymmetric) {
                    Alg->NodeSet[j].C[i] = W;
                    if (W > Alg->M)
                        Alg->M = W;
                }
            }
        }
        break;
    case LKH::UPPER_COL:
        for (j = 2; j <= Alg->Dim; j++) {
            for (i = 1; i < j; i++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                W = round(Alg->Scale * w);
                if (W > std::numeric_limits<int>::max() / 2 / Alg->Precision)
                    W = std::numeric_limits<int>::max() / 2 / Alg->Precision;
                if (Alg->Penalty && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                Alg->NodeSet[j].C[i] = W;
                if (Alg->Asymmetric) {
                    Alg->NodeSet[i].C[j] = W;
                    if (j != i && W > Alg->M)
                        Alg->M = W;
                }
            }
        }
        break;
    case LKH::LOWER_COL:
        for (j = 1; j < Alg->Dim; j++) {
            for (i = j + 1; i <= Alg->Dim; i++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                W = round(Alg->Scale * w);
                if (W > std::numeric_limits<int>::max() / 2 / Alg->Precision)
                    W = std::numeric_limits<int>::max() / 2 / Alg->Precision;
                if (Alg->Penalty && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                Alg->NodeSet[i].C[j] = W;
                if (Alg->Asymmetric) {
                    Alg->NodeSet[j].C[i] = W;
                    if (j != i && W > Alg->M)
                        Alg->M = W;
                }
            }
        }
        break;
    case LKH::UPPER_DIAG_COL:
        for (j = 1; j <= Alg->Dim; j++) {
            for (i = 1; i <= j; i++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                W = round(Alg->Scale * w);
                if (W > std::numeric_limits<int>::max() / 2 / Alg->Precision)
                    W = std::numeric_limits<int>::max() / 2 / Alg->Precision;
                if (Alg->Penalty && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                Alg->NodeSet[j].C[i] = W;
                if (Alg->Asymmetric) {
                    Alg->NodeSet[i].C[j] = W;
                    if (W > Alg->M)
                        Alg->M = W;
                }
            }
        }
        break;
    case LKH::LOWER_DIAG_COL:
        for (j = 1; j <= Alg->Dim; j++) {
            for (i = j; i <= Alg->Dim; i++) {
                if (!fscanf(Alg->ProblemFile, "%lf", &w))
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Missing weight");
                if (j == i)
                    continue;
                W = round(Alg->Scale * w);
                if (W > std::numeric_limits<int>::max() / 2 / Alg->Precision)
                    W = std::numeric_limits<int>::max() / 2 / Alg->Precision;
                if (Alg->Penalty && W < 0)
                    Alg->eprintf("EDGE_WEIGHT_SECTION: Negative weight");
                Alg->NodeSet[i].C[j] = W;
                if (Alg->Asymmetric) {
                    Alg->NodeSet[j].C[i] = W;
                    if (W > Alg->M)
                        Alg->M = W;
                }
            }
        }
        break;
    }
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
    if (Alg->Asymmetric) {
        for (i = 1; i <= Alg->DimensionSaved; i++)
            FixEdge(&Alg->NodeSet[i], &Alg->NodeSet[i + Alg->DimensionSaved]);
        if (Alg->ProblemType == LKH::SOP || Alg->ProblemType == LKH::M1_PDTSP)
            Alg->NodeSet[n].C[1] = 0;
        Alg->Distance = &LKHAlg::Distance_ATSP;
        Alg->WeightType = -1;
    }
}

static void Read_EDGE_WEIGHT_TYPE(LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->EdgeWeightType = Copy(Alg->lkh_strtok(0, Delimiters, &Alg->m_filestream_p))))
        Alg->eprintf("EDGE_WEIGHT_TYPE: string expected");
    for (i = 0; i < strlen(Alg->EdgeWeightType); i++)
		Alg->EdgeWeightType[i] = (char) toupper(Alg->EdgeWeightType[i]);
    if (!strcmp(Alg->EdgeWeightType, "ATT")) {
        Alg->WeightType = LKH::ATT;
		Alg->Distance = &LKHAlg::Distance_ATT;
		Alg->c = &LKHAlg::c_ATT;
		Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "CEIL_2D")) {
        Alg->WeightType = LKH::CEIL_2D;
		Alg->Distance = &LKHAlg::Distance_CEIL_2D;
		Alg->c = &LKHAlg::c_CEIL_2D;
		Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "CEIL_3D")) {
        Alg->WeightType = LKH::CEIL_3D;
		Alg->Distance = &LKHAlg::Distance_CEIL_3D;
		Alg->c = &LKHAlg::c_CEIL_3D;
		Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "EUC_2D") ||
               !strcmp(Alg->EdgeWeightType, "EXACT_2D")) {
        Alg->WeightType = LKH::EUC_2D;
		Alg->Distance = &LKHAlg::Distance_EUC_2D;
		Alg->c = &LKHAlg::c_EUC_2D;
		Alg->CoordType = LKH::TWOD_COORDS;
        if (Alg->Scale == -1 && !strcmp(Alg->EdgeWeightType, "EXACT_2D"))
            Alg->Scale = 1000;
    } else if (!strcmp(Alg->EdgeWeightType, "EUC_3D") ||
               !strcmp(Alg->EdgeWeightType, "EXACT_3D")) {
        Alg->WeightType = LKH::EUC_3D;
		Alg->Distance = &LKHAlg::Distance_EUC_3D;
		Alg->c = &LKHAlg::c_EUC_3D;
		Alg->CoordType = LKH::THREED_COORDS;
        if (Alg->Scale == -1 && !strcmp(Alg->EdgeWeightType, "EXACT_3D"))
            Alg->Scale = 1000;
    } else if (!strcmp(Alg->EdgeWeightType, "EXPLICIT")) {
        Alg->WeightType = LKH::EXPLICIT;
		Alg->Distance = &LKHAlg::Distance_EXPLICIT;
        if (Alg->Scale < 1)
            Alg->Scale = 1;
    } else if (!strcmp(Alg->EdgeWeightType, "FLOOR_2D")) {
        Alg->WeightType = LKH::FLOOR_2D;
		Alg->Distance = &LKHAlg::Distance_FLOOR_2D;
		Alg->c = &LKHAlg::c_FLOOR_2D;
		Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "FLOOR_3D")) {
        Alg->WeightType = LKH::FLOOR_3D;
        Alg->Distance = &LKHAlg::Distance_FLOOR_3D;
        Alg->c = &LKHAlg::c_FLOOR_3D;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "MAN_2D")) {
        Alg->WeightType = LKH::MAN_2D;
        Alg->Distance = &LKHAlg::Distance_MAN_2D;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "MAN_3D")) {
        Alg->WeightType = LKH::MAN_3D;
        Alg->Distance = &LKHAlg::Distance_MAN_3D;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "MAX_2D")) {
        Alg->WeightType = LKH::MAX_2D;
        Alg->Distance = &LKHAlg::Distance_MAX_2D;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "MAX_3D")) {
        Alg->WeightType = LKH::MAX_3D;
        Alg->Distance = &LKHAlg::Distance_MAX_3D;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "GEO")) {
        Alg->WeightType = LKH::GEO;
        Alg->Distance = &LKHAlg::Distance_GEO;
        Alg->c = &LKHAlg::c_GEO;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "GEOM")) {
        Alg->WeightType = LKH::GEOM;
        Alg->Distance = &LKHAlg::Distance_GEOM;
        Alg->c = &LKHAlg::c_GEOM;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "GEO_MEEUS")) {
        Alg->WeightType = LKH::GEO_MEEUS;
        Alg->Distance = &LKHAlg::Distance_GEO_MEEUS;
        Alg->c = &LKHAlg::c_GEO_MEEUS;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "GEOM_MEEUS")) {
        Alg->WeightType = LKH::GEOM_MEEUS;
        Alg->Distance = &LKHAlg::Distance_GEOM_MEEUS;
        Alg->c = &LKHAlg::c_GEOM_MEEUS;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "TOR_2D")) {
        Alg->WeightType = LKH::TOR_2D;
        Alg->Distance = &LKHAlg::Distance_TOR_2D;
        Alg->CoordType = LKH::TWOD_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "TOR_3D")) {
        Alg->WeightType = LKH::TOR_3D;
        Alg->Distance = &LKHAlg::Distance_TOR_3D;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "XRAY1")) {
        Alg->WeightType = LKH::XRAY1;
        Alg->Distance = &LKHAlg::Distance_XRAY1;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "XRAY2")) {
        Alg->WeightType = LKH::XRAY2;
        Alg->Distance = &LKHAlg::Distance_XRAY2;
        Alg->CoordType = LKH::THREED_COORDS;
    } else if (!strcmp(Alg->EdgeWeightType, "SPECIAL")) {
        Alg->WeightType = LKH::SPECIAL;
        Alg->Distance = &LKHAlg::Distance_SPECIAL;
    } else
        Alg->eprintf("Unknown EDGE_WEIGHT_TYPE: %s", Alg->EdgeWeightType);
}

static void Read_FIXED_EDGES_SECTION(LKHAlg *Alg)
{
    LKHAlg::Node *Ni, *Nj, *N, *NPrev = 0, *NNext;
    int i, j, Count = 0;

    CheckSpecificationPart(Alg);
    if (!Alg->FirstNode)
        CreateNodes(Alg);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    if (!Alg->fscanint(Alg->ProblemFile, &i))
        i = -1;
    while (i != -1) {
        if (i <= 0 || i > (Alg->Asymmetric ? Alg->Dimension / 2 : Alg->Dimension))
            Alg->eprintf("FIXED_EDGES_SECTION: Node number out of range: %d",
                    i);
        Alg->fscanint(Alg->ProblemFile, &j);
        if (j <= 0 || j > (Alg->Asymmetric ? Alg->Dimension / 2 : Alg->Dimension))
            Alg->eprintf("FIXED_EDGES_SECTION: Node number out of range: %d",
                    j);
        if (i == j)
            Alg->eprintf("FIXED_EDGES_SECTION: Illegal edge: %d to %d", i, j);
        Ni = &Alg->NodeSet[i];
        Nj = &Alg->NodeSet[Alg->Asymmetric ? j + Alg->Dimension / 2 : j];
        if (!FixEdge(Ni, Nj))
            Alg->eprintf("FIXED_EDGES_SECTION: Illegal fix: %d to %d", i, j);
        /* Cycle check */
        N = Ni;
        Count = 0;
        do {
            NNext = N->FixedTo1 != NPrev ? N->FixedTo1 : N->FixedTo2;
            NPrev = N;
            Count++;
        } while ((N = NNext) && N != Ni);
        if (N == Ni && Count != Alg->Dimension)
            Alg->eprintf("FIXED_EDGES_SECTION: Illegal fix: %d to %d", i, j);
        if (!Alg->fscanint(Alg->ProblemFile, &i))
            i = -1;
    }
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
}

static void Read_GRID_SIZE(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p);

    if (!Token || !sscanf(Token, "%lf", &Alg->GridSize))
        Alg->eprintf("GRID_SIZE: Real expected");
    if (Alg->GridSize < 0)
        Alg->eprintf("GRID_SIZE: non-negative Real expected");
}

static void Read_NODE_COORD_SECTION(LKHAlg *Alg)
{
    LKHAlg::Node *N;
    int Id, i;

    CheckSpecificationPart(Alg);
    if (Alg->CoordType != LKH::TWOD_COORDS && Alg->CoordType != LKH::THREED_COORDS)
        Alg->eprintf("NODE_COORD_SECTION conflicts with NODE_COORD_TYPE: %s",
			Alg->NodeCoordType);
    if (!Alg->FirstNode)
        CreateNodes(Alg);
    N = Alg->FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != Alg->FirstNode);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    for (i = 1; i <= Alg->Dim; i++) {
        if (!Alg->fscanint(Alg->ProblemFile, &Id))
            Alg->eprintf("NODE_COORD_SECTION: Missing nodes");
        if (Id <= 0 || Id > Alg->Dimension)
            Alg->eprintf("NODE_COORD_SECTION: Node number out of range: %d",
                    Id);
        N = &Alg->NodeSet[Id];
        if (N->V == 1)
            Alg->eprintf("NODE_COORD_SECTION: Node number occurs twice: %d",
                    N->Id);
        N->V = 1;
        if (!fscanf(Alg->ProblemFile, "%lf", &N->X))
            Alg->eprintf("NODE_COORD_SECTION: Missing X-coordinate");
        if (!fscanf(Alg->ProblemFile, "%lf", &N->Y))
            Alg->eprintf("NODE_COORD_SECTION: Missing Y-coordinate");
        if (Alg->CoordType == LKH::THREED_COORDS
            && !fscanf(Alg->ProblemFile, "%lf", &N->Z))
            Alg->eprintf("NODE_COORD_SECTION: Missing Z-coordinate");
        if (Alg->Name && !strcmp(Alg->Name, "d657")) {
            N->X = (float) N->X;
            N->Y = (float) N->Y;
        }
    }
    N = Alg->FirstNode;
    do
        if (!N->V && N->Id <= Alg->Dim)
            break;
    while ((N = N->Suc) != Alg->FirstNode);
    if (!N->V)
        Alg->eprintf("NODE_COORD_SECTION: No coordinates given for node %d",
                N->Id);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
    if (Alg->Asymmetric)
        Convert2FullMatrix(Alg);
}

static void Read_NODE_COORD_TYPE(LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->NodeCoordType = Copy(Alg->lkh_strtok(0, Delimiters, &Alg->m_filestream_p))))
        Alg->eprintf("NODE_COORD_TYPE: string expected");
    for (i = 0; i < strlen(Alg->NodeCoordType); i++)
		Alg->NodeCoordType[i] = (char) toupper(Alg->NodeCoordType[i]);
    if (!strcmp(Alg->NodeCoordType, "TWOD_COORDS"))
        Alg->CoordType = LKH::TWOD_COORDS;
    else if (!strcmp(Alg->NodeCoordType, "THREED_COORDS"))
        Alg->CoordType = LKH::THREED_COORDS;
    else if (!strcmp(Alg->NodeCoordType, "NO_COORDS"))
        Alg->CoordType = LKH::NO_COORDS;
    else
        Alg->eprintf("Unknown NODE_COORD_TYPE: %s", Alg->NodeCoordType);
}

static void Read_PICKUP_AND_DELIVERY_SECTION(LKHAlg *Alg)
{
    int Id, i;
    LKHAlg::Node *N = Alg->FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != Alg->FirstNode);
    for (i = 1; i <= Alg->Dim; i++) {
        if (!Alg->fscanint(Alg->ProblemFile, &Id))
            Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: Missing nodes");
        if (Id <= 0 || Id > Alg->Dim)
            Alg->eprintf
                ("PICKUP_AND_DELIVERY_SECTION: Node number out of range: %d",
                 Id);
        N = &Alg->NodeSet[Id];
        if (N->V == 1)
            Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                    "Node number occurs twice: %d", N->Id);
        N->V = 1;
        if (!fscanf(Alg->ProblemFile, "%d %lf %lf %lf %d %d",
                    &N->Demand, &N->Earliest, &N->Latest, &N->ServiceTime,
                    &N->Pickup, &N->Delivery))
            Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                    " Missing data for node %d", N->Id);
        if (N->ServiceTime < 0)
            Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                    "Negative Service Time for node %d", N->Id);
        if (N->Earliest > N->Latest)
            Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                    "Earliest > Latest for node %d", N->Id);
    }
    N = Alg->FirstNode;
    do
        if (!N->V && N->Id <= Alg->Dim)
            break;
    while ((N = N->Suc) != Alg->FirstNode);
    if (!N->V)
        Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: No data given for node %d",
                N->Id);
    if (Alg->ProblemType != LKH::VRPPD) {
        do {
            if (N->Delivery) {
                if (Alg->NodeSet[N->Delivery].Pickup != N->Id ||
                    N->Delivery == N->Id)
                    Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                            "Illegal pairing for node %d", N->Id);
                if (N->Demand < 0)
                    Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                            "Negative demand for delivery node %d", N->Id);
            } else if (N->Pickup) {
                if (Alg->NodeSet[N->Pickup].Delivery != N->Id
                    || N->Pickup == N->Id)
                    Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                            "Illegal pairing for node %d", N->Id);
                if (N->Demand > 0)
                    Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                            "Positive demand for pickup node %d", N->Id);
                if (N->Demand + Alg->NodeSet[N->Pickup].Demand)
                    Alg->eprintf("PICKUP_AND_DELIVERY_SECTION: "
                            "Demand for pickup node %d and demand for delivery "
                            "node %d does not sum to zero", N->Id,
                            N->Pickup);
            }
        } while ((N = N->Suc) != Alg->FirstNode);
    }
}

static void Read_TIME_WINDOW_SECTION(LKHAlg *Alg)
{
    int Id, i;
    LKHAlg::Node *N = Alg->FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != Alg->FirstNode);
    for (i = 1; i <= Alg->Dim; i++) {
        if (!Alg->fscanint(Alg->ProblemFile, &Id))
            Alg->eprintf("TIME_WINDOW_SECTION: Missing nodes");
        if (Id <= 0 || Id > Alg->Dim)
            Alg->eprintf("TIME_WINDOW_SECTION: Node number out of range: %d",
                    Id);
        N = &Alg->NodeSet[Id];
        if (N->V == 1)
            Alg->eprintf("TIME_WINDOW_SECTION: Node number occurs twice: %d",
                    N->Id);
        N->V = 1;
        if (!fscanf(Alg->ProblemFile, "%lf", &N->Earliest))
            Alg->eprintf("TIME_WINDOW_SECTION: Missing earliest time");
        if (!fscanf(Alg->ProblemFile, "%lf", &N->Latest))
            Alg->eprintf("TIME_WINDOW_SECTION: Missing latest time");
        if (N->Earliest > N->Latest)
			Alg->printff("%d: %f %f\n", N->Id, N->Earliest, N->Latest);
        if (N->Earliest > N->Latest)
            Alg->eprintf("TIME_WINDOW_SECTION: Earliest > Latest for node %d",
                    N->Id);
    }
    N = Alg->FirstNode;
    do
        if (!N->V && N->Id <= Alg->Dim)
            break;
    while ((N = N->Suc) != Alg->FirstNode);
    if (!N->V)
        Alg->eprintf("TIME_WINDOW_SECTION: No time window given for node %d",
                N->Id);
}

static void Read_TOUR_SECTION(FILE ** File, LKHAlg *Alg)
{
    LKHAlg::Node *First = 0, *Last = 0, *N, *Na;
    int i, k;

    if (Alg->TraceLevel >= 1) {
		Alg->printff("Reading ");
        if (File == &Alg->InitialTourFile)
			Alg->printff("INITIAL_TOUR_FILE: \"%s\" ... ", Alg->InitialTourFileName);
        else if (File == &Alg->InputTourFile)
			Alg->printff("INPUT_TOUR_FILE: \"%s\" ... ", Alg->InputTourFileName);
        else if (File == &Alg->SubproblemTourFile)
			Alg->printff("SUBPROBLEM_TOUR_FILE: \"%s\" ... ",
				Alg->SubproblemTourFileName);
        else
            for (i = 0; i < Alg->MergeTourFiles; i++)
                if (File == &Alg->MergeTourFile[i])
					Alg->printff("MERGE_TOUR_FILE: \"%s\" ... ",
						Alg->MergeTourFileName[i]);
    }
    if (!Alg->FirstNode)
        CreateNodes(Alg);
    N = Alg->FirstNode;
    do
        N->V = 0;
    while ((N = N->Suc) != Alg->FirstNode);
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension--;
    if (Alg->Asymmetric)
        Alg->Dimension = Alg->DimensionSaved;
    int b = 0;
    if (!Alg->fscanint(*File, &i))
        i = -1;
    else if (i == 0) {
        b = 1;
        i++;
    }
    for (k = 0; k <= Alg->Dimension && i != -1; k++) {
        if (i <= 0 || i > Alg->Dimension)
            Alg->eprintf("(TOUR_SECTION) Node number out of range: %d", i);
        N = &Alg->NodeSet[i];
        if (N->V == 1 && k != Alg->Dimension)
            Alg->eprintf("(TOUR_SECTION) Node number occurs twice: %d", N->Id);
        N->V = 1;
        if (k == 0)
            First = Last = N;
        else {
            if (Alg->Asymmetric) {
                Na = N + Alg->Dimension;
                Na->V = 1;
            } else
                Na = 0;
            if (File == &Alg->InitialTourFile) {
                if (!Na)
                    Last->InitialSuc = N;
                else {
                    Last->InitialSuc = Na;
                    Na->InitialSuc = N;
                }
            } else if (File == &Alg->InputTourFile) {
                if (!Na)
                    Last->InputSuc = N;
                else {
                    Last->InputSuc = Na;
                    Na->InputSuc = N;
                }
            } else if (File == &Alg->SubproblemTourFile) {
                if (!Na)
                    (Last->SubproblemSuc = N)->SubproblemPred = Last;
                else {
                    (Last->SubproblemSuc = Na)->SubproblemPred = Last;
                    (Na->SubproblemSuc = N)->SubproblemPred = Na;
                }
            } else {
                for (i = 0; i < Alg->MergeTourFiles; i++) {
                    if (File == &Alg->MergeTourFile[i]) {
                        if (!Na) {
                            Last->MergeSuc[i] = N;
                            if (i == 0)
                                N->MergePred = Last;
                        } else {
                            Last->MergeSuc[i] = Na;
                            Na->MergeSuc[i] = N;
                            if (i == 0) {
                                Na->MergePred = Last;
                                N->MergePred = Na;
                            }
                        }
                    }
                }
            }
            Last = N;
        }
        if (k < Alg->Dimension) {
            Alg->fscanint(*File, &i);
            if (b)
                if (i >= 0)
                    i++;
        }
        if (k == Alg->Dimension - 1)
            i = First->Id;
    }
    N = Alg->FirstNode;
    do {
        if (!N->V)
            Alg->eprintf("TOUR_SECTION: Node is missing: %d", N->Id);
    } while ((N = N->Suc) != Alg->FirstNode);
    if (File == &Alg->SubproblemTourFile) {
        do {
            if (N->FixedTo1 &&
                N->SubproblemPred != N->FixedTo1
                && N->SubproblemSuc != N->FixedTo1)
                Alg->eprintf("Fixed edge (%d, %d) "
                        "does not belong to subproblem tour", N->Id,
                        N->FixedTo1->Id);
            if (N->FixedTo2 && N->SubproblemPred != N->FixedTo2
                && N->SubproblemSuc != N->FixedTo2)
                Alg->eprintf("Fixed edge (%d, %d) "
                        "does not belong to subproblem tour", N->Id,
                        N->FixedTo2->Id);
        } while ((N = N->Suc) != Alg->FirstNode);
    }
    if (Alg->ProblemType == LKH::HPP)
        Alg->Dimension++;
    if (Alg->Asymmetric)
        Alg->Dimension *= 2;
    if (Alg->TraceLevel >= 1)
		Alg->printff("done\n");
}

static void Read_TYPE(LKHAlg *Alg)
{
    unsigned int i;

    if (!(Alg->Type = Copy(Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p))))
        Alg->eprintf("TYPE: string expected");
    for (i = 0; i < strlen(Alg->Type); i++)
		Alg->Type[i] = (char) toupper(Alg->Type[i]);
    if (!strcmp(Alg->Type, "TSP"))
        Alg->ProblemType = LKH::TSP;
    else if (!strcmp(Alg->Type, "ATSP"))
        Alg->ProblemType = LKH::ATSP;
    else if (!strcmp(Alg->Type, "SOP"))
        Alg->ProblemType = LKH::SOP;
    else if (!strcmp(Alg->Type, "HCP"))
        Alg->ProblemType = LKH::HCP;
    else if (!strcmp(Alg->Type, "HPP"))
        Alg->ProblemType = LKH::HPP;
    else if (!strcmp(Alg->Type, "BWTSP"))
        Alg->ProblemType = LKH::BWTSP;
    else if (!strcmp(Alg->Type, "CCVRP"))
        Alg->ProblemType = LKH::CCVRP;
    else if (!strcmp(Alg->Type, "CVRP") || !strcmp(Alg->Type, "DCVRP"))
        Alg->ProblemType = LKH::CVRP;
    else if (!strcmp(Alg->Type, "ACVRP"))
        Alg->ProblemType = LKH::ACVRP;
    else if (!strcmp(Alg->Type, "CVRPTW"))
        Alg->ProblemType = LKH::CVRPTW;
    else if (!strcmp(Alg->Type, "MLP"))
        Alg->ProblemType = LKH::MLP;
    else if (!strcmp(Alg->Type, "OVRP"))
        Alg->ProblemType = LKH::OVRP;
    else if (!strcmp(Alg->Type, "PDPTW"))
        Alg->ProblemType = LKH::PDPTW;
    else if (!strcmp(Alg->Type, "PDTSP"))
        Alg->ProblemType = LKH::PDTSP;
    else if (!strcmp(Alg->Type, "PDTSPF") || !strcmp(Alg->Type, "PDTSPF"))
        Alg->ProblemType = LKH::PDTSPF;
    else if (!strcmp(Alg->Type, "PDTSPL") || !strcmp(Alg->Type, "PDTSPL"))
        Alg->ProblemType = LKH::PDTSPL;
    else if (!strcmp(Alg->Type, "TRP") || !strcmp(Alg->Type, "MTRP") ||
             !strcmp(Alg->Type, "MTRPD"))
        Alg->ProblemType = LKH::TRP;
    else if (!strcmp(Alg->Type, "RCTVRP"))
        Alg->ProblemType = LKH::RCTVRP;
    else if (!strcmp(Alg->Type, "RCTVRPTW"))
        Alg->ProblemType = LKH::RCTVRPTW;
    else if (!strcmp(Alg->Type, "TSPTW"))
        Alg->ProblemType = LKH::TSPTW;
    else if (!strcmp(Alg->Type, "VRPB"))
        Alg->ProblemType = LKH::VRPB;
    else if (!strcmp(Alg->Type, "VRPBTW"))
        Alg->ProblemType = LKH::VRPBTW;
    else if (!strcmp(Alg->Type, "VRPSPD") ||
             !strcmp(Alg->Type, "VRPSPDTW") ||
             !strcmp(Alg->Type, "VRPMPD") ||
             !strcmp(Alg->Type, "VRPMPDTW") || !strcmp(Alg->Type, "MVRPB"))
        Alg->ProblemType = LKH::VRPPD;
    else if (!strcmp(Alg->Type, "1-PDTSP"))
        Alg->ProblemType = LKH::ONE_PDTSP;
    else if (!strcmp(Alg->Type, "M-PDTSP"))
        Alg->ProblemType = LKH::M_PDTSP;
    else if (!strcmp(Alg->Type, "M1-PDTSP"))
        Alg->ProblemType = LKH::M1_PDTSP;
    else if (!strcmp(Alg->Type, "TSPDL"))
        Alg->ProblemType = LKH::TSPDL;
    else if (!strcmp(Alg->Type, "CTSP"))
        Alg->ProblemType = LKH::CTSP;
    else if (!strcmp(Alg->Type, "TOUR")) {
        Alg->ProblemType = LKH::TOUR;
        Alg->eprintf("TYPE: Type not implemented: %s", Alg->Type);
    } else
        Alg->eprintf("Unknown TYPE: %s", Alg->Type);
    Alg->Asymmetric =
        Alg->ProblemType == LKH::ATSP ||
        Alg->ProblemType == LKH::CCVRP ||
        Alg->ProblemType == LKH::ACVRP ||
        Alg->ProblemType == LKH::CVRPTW ||
        Alg->ProblemType == LKH::MLP ||
        Alg->ProblemType == LKH::M_PDTSP ||
        Alg->ProblemType == LKH::M1_PDTSP ||
        Alg->ProblemType == LKH::ONE_PDTSP ||
        Alg->ProblemType == LKH::OVRP ||
        Alg->ProblemType == LKH::PDTSP ||
        Alg->ProblemType == LKH::PDTSPF ||
        Alg->ProblemType == LKH::PDTSPL ||
        Alg->ProblemType == LKH::PDPTW ||
        Alg->ProblemType == LKH::RCTVRP ||
        Alg->ProblemType == LKH::RCTVRPTW ||
        Alg->ProblemType == LKH::SOP ||
        Alg->ProblemType == LKH::TRP ||
        Alg->ProblemType == LKH::TSPDL ||
        Alg->ProblemType == LKH::TSPTW ||
        Alg->ProblemType == LKH::VRPB ||
        Alg->ProblemType == LKH::VRPBTW || Alg->ProblemType == LKH::VRPPD;
}

static void Read_SERVICE_TIME(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p);

    if (!Token || !sscanf(Token, "%lf", &Alg->ServiceTime))
        Alg->eprintf("SERVICE_TIME: Real expected");
    if (Alg->ServiceTime < 0)
        Alg->eprintf("SERVICE_TIME: < 0");
}

static void Read_SERVICE_TIME_SECTION(LKHAlg *Alg)
{
    int Id, i;
	LKHAlg::Node *N;

    for (i = 1; i <= Alg->Dim; i++) {
        Alg->fscanint(Alg->ProblemFile, &Id);
        if (Id <= 0 || Id > Alg->Dim)
            Alg->eprintf("SERVICE_TIME_SECTION: Node number out of range: %d",
                    Id);
        N = &Alg->NodeSet[Id];
        if (!fscanf(Alg->ProblemFile, "%lf", &N->ServiceTime))
            Alg->eprintf("SERVICE_TIME_SECTION: "
                    "Missing service time for node %d", Id);
    }
}

/*
 The ReadTour function reads a tour from a file.
 
 The format is as follows:
 
 OPTIMUM = <Real>
 Known optimal tour length. A run will be terminated as soon as a tour
 length less than or equal to optimum is achieved.
 Default: MINUS_INFINITY.
 
 TOUR_SECTION :
 A tour is specified in this section. The tour is given by a list of integers
 giving the sequence in which the nodes are visited in the tour. The tour is
 terminated by a -1.
 
 EOF
 Terminates the input data. The entry is optional.
 
 Other keywords in TSPLIB format may be included in the file, but they are
 ignored.
 */

void LKHAlg::ReadTour(char *FileName, FILE ** File)
{
    char *Line, *Keyword, *Token;
    unsigned int i;
    int Done = 0;

    if (!(*File = fopen(FileName, "r")))
        eprintf("Cannot open tour file: \"%s\"", FileName);
    while ((Line = ReadLine(*File))) {
        if (!(Keyword = lkh_strtok(Line, Delimiters, &m_filestream_p)))
            continue;
        for (i = 0; i < strlen(Keyword); i++)
            Keyword[i] = (char) toupper(Keyword[i]);
        if (!strcmp(Keyword, "COMMENT") ||
            !strcmp(Keyword, "DEMAND_SECTION") ||
            !strcmp(Keyword, "DEPOT_SECTION") ||
            !strcmp(Keyword, "DISPLAY_DATA_SECTION") ||
            !strcmp(Keyword, "DISPLAY_DATA_TYPE") ||
            !strcmp(Keyword, "EDGE_DATA_FORMAT") ||
            !strcmp(Keyword, "EDGE_DATA_SECTION") ||
            !strcmp(Keyword, "EDGE_WEIGHT_FORMAT") ||
            !strcmp(Keyword, "EDGE_WEIGHT_SECTION") ||
            !strcmp(Keyword, "EDGE_WEIGHT_TYPE") ||
            !strcmp(Keyword, "FIXED_EDGES_SECTION") ||
            !strcmp(Keyword, "NAME") ||
            !strcmp(Keyword, "NODE_COORD_SECTION") ||
            !strcmp(Keyword, "NODE_COORD_TYPE")
            || !strcmp(Keyword, "TYPE"));
        else if (strcmp(Keyword, "OPTIMUM") == 0) {
            if (!(Token = lkh_strtok(0, Delimiters, &m_filestream_p)) ||
                !sscanf(Token, GainInputFormat, &Optimum))
                eprintf("[%s] (OPTIMUM): Integer expected", FileName);
        } else if (strcmp(Keyword, "DIMENSION") == 0) {
            int Dim = 0;
            if (!(Token = lkh_strtok(0, Delimiters, &m_filestream_p)) ||
                !sscanf(Token, "%d", &Dim))
                eprintf("[%s] (DIMENSION): Integer expected", FileName);
            if (Dim != DimensionSaved && Dim != Dimension) {
				printff("Dim = %d, DimensionSaved = %d, Dimension = %d\n",
                        Dim, DimensionSaved, Dimension);
                eprintf
                    ("[%s] (DIMENSION): does not match problem dimension",
                     FileName);
            }
        } else if (!strcmp(Keyword, "PICKUP_AND_DELIVERY_SECTION")) {
            Read_PICKUP_AND_DELIVERY_SECTION(this);
        } else if (!strcmp(Keyword, "SERVICE_TIME_SECTION")) {
            Read_SERVICE_TIME_SECTION(this);
        } else if (!strcmp(Keyword, "TOUR_SECTION")) {
            Read_TOUR_SECTION(File, this);
            Done = 1;
        } else if (!strcmp(Keyword, "EOF"))
            break;
        else
            eprintf("[%s] Unknown Keyword: %s", FileName, Keyword);
    }
    if (!Done)
        eprintf("Missing TOUR_SECTION in tour file: \"%s\"", FileName);
    fclose(*File);
}

static void Read_RISK_THRESHOLD(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p);

    if (!Token || !sscanf(Token, "%d", &Alg->RiskThreshold))
        Alg->eprintf("RISK_THRESHOLD: Integer expected");
}

static void Read_SALESMEN(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p);

    if (!Token || (Alg->Salesmen == 1 && !sscanf(Token, "%d", &Alg->Salesmen)))
        Alg->eprintf("SALESMEN/VEHICLES: Integer expected");
    if (Alg->Salesmen <= 0)
        Alg->eprintf("SALESMEN/VEHICLES: <= 0");
}

static void Read_SCALE(LKHAlg *Alg)
{
    char *Token = Alg->lkh_strtok(0, Delimiters,&Alg->m_filestream_p);

    if (!Token || !sscanf(Token, "%d", &Alg->Scale))
        Alg->eprintf("SCALE: Integer expected");
    if (Alg->Scale < 1)
        Alg->eprintf("SCALE: < 1");
}

static void Convert2FullMatrix(LKHAlg *Alg)
{
    int n = Alg->DimensionSaved, i, j;
	LKHAlg::Node *Ni, *Nj;

    if (Alg->Scale < 1)
        Alg->Scale = 1;
    if (n > Alg->MaxMatrixDimension) {
		Alg->OldDistance = Alg->Distance;
        Alg->Distance = &LKHAlg::Distance_Asymmetric;
        for (i = 1; i <= n; i++) {
            Ni = &Alg->NodeSet[i];
            Nj = &Alg->NodeSet[i + n];
            Nj->X = Ni->X;
            Nj->Y = Ni->Y;
            Nj->Z = Ni->Z;
            FixEdge(Ni, Nj);
        }
        return;
    }
    assert(Alg->CostMatrix = (int *) calloc((size_t) n * n, sizeof(int)));
    for (i = 1; i <= n; i++) {
        Ni = &Alg->NodeSet[i];
        Ni->C = &Alg->CostMatrix[(size_t) (i - 1) * n] - 1;
    }
    for (i = 1; i <= Alg->Dim; i++) {
        Ni = &Alg->NodeSet[i];
        for (j = i + 1; j <= Alg->Dim; j++) {
            Nj = &Alg->NodeSet[j];
            Ni->C[j] = Nj->C[i] = (Alg->*(Alg->Distance))(Ni, Nj);
        }
    }
    for (i = 1; i <= n; i++)
        FixEdge(&Alg->NodeSet[i], &Alg->NodeSet[i + n]);
    Alg->c = 0;
    Alg->Distance = &LKHAlg::Distance_ATSP;
    Alg->WeightType = -1;
}
}