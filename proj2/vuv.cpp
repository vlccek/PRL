#include <array>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <variant>
#include <vector>
#include <algorithm>
#include <cassert>
#include <map>

using Edge = std::pair<int, int>;


int getCountOfRank()
{
    int numProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    return numProcesses;
}

int getRank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}



int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <binary_tree>" << std::endl;
        std::cerr << "Please provide a binary tree as a command line argument." << std::endl;
        std::cerr << "Example: " << argv[0] << "148" << std::endl;
        return 1;
    }

    MPI_Init(&argc, &argv);

    int rank = getRank();
    int rankCount = getCountOfRank();


    auto treeReprAsString = std::string(argv[1]);
    int treeElementCount = treeReprAsString.length();

    // edge cases
    if (treeElementCount <= 1)
    {
        std::cerr << "Error: The binary tree is withotu edges. " << std::endl;
        MPI_Finalize();
        return 1;
    }

    if (treeElementCount <= 2)
    {
        printf("%c:%d,%c:%d", treeReprAsString[0], 0, treeReprAsString[1], 1);
    }


    std::vector<std::pair<int, int>> edgeList;
    std::vector<std::vector<int>> edgeTreeReprezentation;
    /* reprezetens tree as follow:
     * var[0] = {(0, 1), (0, 2)} // means that thers 2 edges from 0 to 1 and 2
     */
// TODO: vymyslet jak tohle počítat paralelně.
    // its is faster to compute this than to gather that from all ranks
    for (int i = 0; i < treeElementCount; ++i)
    {
        int RchildIdx = (2 * i + 2 < treeElementCount) ? 2 * i + 2 : -1;
        int LchildIdx = (2 * i + 1 < treeElementCount) ? 2 * i + 1 : -1;
        int parentIdx = (i > 0) ? (i - 1) / 2 : -1;

        const int numOfChildrens = (LchildIdx != -1) + (RchildIdx != -1);


        std::vector<int> tmp;


        if (parentIdx != -1)
        {
            edgeList.emplace_back(i, parentIdx);
            tmp.emplace_back(edgeList.size() - 1);
        }


        if (LchildIdx != -1)
        {
            edgeList.emplace_back(i, LchildIdx);

            tmp.emplace_back(edgeList.size() - 1);
        }

        if (RchildIdx != -1)
        {
            edgeList.emplace_back(i, RchildIdx);
            tmp.emplace_back(edgeList.size() - 1);
        }


        edgeTreeReprezentation.push_back(tmp);
    }

    if (edgeList.size() != rankCount)
    {
        std::fprintf(stderr,"Edge list size: %d and number of proccess %d. This number should by equal! \n", edgeList.size(), rankCount);
        return 2;
    }


    std::map<std::pair<int, int>, int> edgeToIndex;
    for (int i = 0; i < edgeList.size(); ++i)
    {
        edgeToIndex[edgeList[i]] = i;
    }

    std::vector<int> backEdgePartnerIndices(edgeList.size());

    for (int i = 0; i < edgeList.size(); ++i)
    {
        // Získání aktuální hrany (u, v)
        std::pair<int, int> currentEdge = edgeList[i];
        int u = currentEdge.first;
        int v = currentEdge.second;

        // Definování opačné (zpětné) hrany (v, u)
        std::pair<int, int> reverseEdge = {v, u};

        // Nalezení indexu opačné hrany v mapě
        auto it = edgeToIndex.find(reverseEdge);
        if (it != edgeToIndex.end())
        {
            // Opačná hrana byla nalezena, její index je it->second
            backEdgePartnerIndices[i] = it->second;
        }
    }
    // real parallel section

    /** etour*/
    std::pair<int, int> myEdge = edgeList[rank];
    int myEdgeIndex = edgeToIndex[myEdge];
    std::pair<int, int> reverseEdge = {myEdge.second, myEdge.first};


    auto it_rev = edgeToIndex.find(reverseEdge);
    int reverseEdgeIndex = it_rev->second;
    int original_successor;

    int vertexNumber = reverseEdge.first;

    auto& nextElements = edgeTreeReprezentation[vertexNumber];

    std::array<int, 2> myEtourResult = {myEdgeIndex, 0}; // second element is just place holder here

    int reverseEdgePos = -1;
    for (int i = 0; i < nextElements.size(); ++i)
    {
        if (nextElements[i] == reverseEdgeIndex)
        {
            reverseEdgePos = i; // Našli jsme pozici
            break;
        }
    }

    if (reverseEdgePos != -1)
    {
        int successorPos = (reverseEdgePos + 1) % nextElements.size(); // Cyklický posun
        int successorEdgeIndex = nextElements[successorPos]; // Získání skutečného indexu následnické hrany
        original_successor = nextElements[successorPos];

        myEtourResult[1] = successorEdgeIndex; // Uložení správného indexu následníka
    }
    // std::printf("%d: myEtour = (%d, %d)\n", rank, myEtourResult[0], myEtourResult[1]);
    // std::printf("id: %d: myEtour = ((%d, %d), (%d, %d)\n", rank, edgeList[myEtourResult[0]].first, edgeList[myEtourResult[0]].second, edgeList[myEtourResult[0]].second, edgeList[myEtourResult[1]].second);


    //*end etour*/
    // root computing:
    int N = edgeList.size();

    auto rightChildRoot = edgeList[1];
    std::pair<int, int> rightChildRootReverse = {rightChildRoot.second, rightChildRoot.first};
    int breakEdgeIndex = edgeToIndex[rightChildRootReverse];


    // ordering etour -- List ranking (Wyllie’s algorithm)


    int my_succ; // AktuÃ¡lny nÃ¡slednÃ­k pre List Ranking
    int my_dist; // AktuÃ¡lna vzdialenosÅ¥ KU KONCU

    // SprÃ¡vna InicializÃ¡cia pre vzdialenosÅ¥ KU KONCU
    my_succ = original_successor; // PouÅ¾i predtÃ½m vypoÄÃ­tanÃ©ho nÃ¡slednÃ­ka
    my_dist = (rank == breakEdgeIndex) ? 0 : 1; // Koniec mÃ¡ 0, ostatnÃ­ 1

    // SprÃ¡vna Ãºprava pre koniec zoznamu
    if (rank == breakEdgeIndex)
    {
        my_succ = breakEdgeIndex; // Proces na konci ukazuje sÃ¡m na seba
    }

    int n_steps = std::ceil(std::log2(edgeList.size()));
    std::vector<int> all_dist;
    all_dist.resize(N);
    std::vector<int> all_succ;
    all_succ.resize(N);


    for (int i = 0; i < n_steps; ++i)
    {
        MPI_Allgather(&my_dist, 1, MPI_INT, all_dist.data(), 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&my_succ, 1, MPI_INT, all_succ.data(), 1, MPI_INT, MPI_COMM_WORLD);

        if (my_succ != breakEdgeIndex)
        {
            int current_succ_idx = my_succ;
            if (current_succ_idx >= 0 && current_succ_idx < N)
            {
                // BezpeÄnostnÃ¡ kontrola
                my_dist = my_dist + all_dist[current_succ_idx]; // SÄÃ­tame vzdialenosti ku koncu
                my_succ = all_succ[current_succ_idx]; // Skok
            }
            else
            {
                my_succ = breakEdgeIndex; // ZastavÃ­me, ak je index neplatnÃ½
            }
        }
    }

    int finalPosn = N - 1 - my_dist;

    std::vector<int> allFinalPosn;
    allFinalPosn.resize(N);

    MPI_Allgather(&finalPosn, 1, MPI_INT, allFinalPosn.data(), 1, MPI_INT, MPI_COMM_WORLD);


    // std::printf("Ordering of (rank: %d) edges: %d \n", rank, finalPosn);

    //* start
    int my_initial_weight = 0;
    int u = myEdge.first;
    int v = myEdge.second;
    int root_vertex_index = 0;


    bool is_forward = false;
    if (u == root_vertex_index && v != root_vertex_index)
    {
        is_forward = true;
    }
    else if (v != root_vertex_index && u == (v - 1) / 2)
    {
        is_forward = true;
    }

    if (is_forward)
    {
        my_initial_weight = -1; // Forward hrana (dole)
    }
    else
    {
        my_initial_weight = +1; // Backward hrana (hore)
    }

    // std::printf("W: %d: myedge=(%d,%d) w=%d\n", rank, u, v, my_initial_weight);


    std::vector<int> weights;
    weights.resize(N);

    MPI_Allgather(&my_initial_weight, 1, MPI_INT, weights.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> orderedWeights(N);
    for (int r = 0; r < N; ++r) { // PrechÃ¡dzaj RANKY procesov r
        int pos = allFinalPosn[r]; // PozÃ­cia procesu r
        if (pos >= 0 && pos < N) {
             orderedWeights[pos] = weights[r]; // Na pozÃ­ciu uloÅ¾ vÃ¡hu procesu r
        }
        // Inak ignoruj (nemalo by nastaÅ¥ pri sprÃ¡vnom ranku)
    }

    // 2. VypoÄÃ­taj VÅ ETKY sufixovÃ© sÃºÄty sekvenÄne (lokÃ¡lne na kaÅ¾dom procese)
    std::vector<int> suffixSumsOrdered(N);
    if (N > 0) {
        suffixSumsOrdered[N - 1] = orderedWeights[N - 1]; // PoslednÃ½ prvok
        for (int i = N - 2; i >= 0; --i) { // Iteruj od predposlednÃ©ho k prvÃ©mu
            suffixSumsOrdered[i] = orderedWeights[i] + suffixSumsOrdered[i + 1];
        }
    }

    // 3. Zisti finÃ¡lny sufixovÃ½ sÃºÄet pre TENTO proces 'rank'
    int myFinalPosn = allFinalPosn[rank]; // Zisti pozÃ­ciu tohto procesu
    int my_final_suffix_sum = -999; // Default hodnota pre prÃ­pad chyby
    if (myFinalPosn >= 0 && myFinalPosn < N) {
         my_final_suffix_sum = suffixSumsOrdered[myFinalPosn];
    }

    // 4. VypoÄÃ­taj finÃ¡lnu ÃºroveÅ, ak mÃ¡ tento proces forward hranu
    int my_level = -1;
    // 'is_forward' a 'v = myEdge.second' musia byÅ¥ vypoÄÃ­tanÃ© predtÃ½m sprÃ¡vne
    if (is_forward) {
        my_level = my_final_suffix_sum + 1;
    }

    // --- ZhromaÅ¾denie a vÃ½pis finÃ¡lnych ÃºrovnÃ­ ---
    // (TÃºto ÄasÅ¥ uÅ¾ mÃ¡Å¡ na konci svojho kompletnÃ©ho kÃ³du, pouÅ¾ije my_level a v)
    std::vector<int> node_levels(treeElementCount, -1);
    int send_data[2] = {-1, -1};
    if (is_forward) {
        send_data[0] = v;
        send_data[1] = my_level;
    }
    std::vector<int> all_data(2 * N);
    MPI_Allgather(send_data, 2, MPI_INT, all_data.data(), 2, MPI_INT, MPI_COMM_WORLD);

    node_levels[root_vertex_index] = 0;
    for (int i = 0; i < N; ++i) {
        int node_v = all_data[2 * i];
        int level_v = all_data[2 * i + 1];
        if (node_v != -1) node_levels[node_v] = level_v;
    }

    // VÃ½pis finÃ¡lneho poÄ¾a ÃºrovnÃ­ (len na rank 0)
    if (rank == 0) {
        for(int i=0; i < treeElementCount; ++i) {
            std::printf("%c:%d", treeReprAsString[i], node_levels[i]);
            if (i < treeElementCount - 1) {
                std::printf(",");
            }
        }
        std::printf("\n");
    }


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
