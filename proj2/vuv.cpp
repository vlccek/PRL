/**
 * @file vuv.cpp
 *
 * @author Jakub Šebek
 * @date 2025-04-14
 */

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
        // Print usage instructions
        fprintf(stderr, "Usage: %s <binary_tree>\n", argv[0]);
        fprintf(stderr, "Please provide a binary tree as a command line argument.\n");
        fprintf(stderr, "Example: %s 148\n", argv[0]); // Assuming "148" is a static example part
        return 1; // Exit with an error code.
    }
    // Initialize the MPI environment.
    // prepare some variables and compute base values
    MPI_Init(&argc, &argv);
    int rank = getRank();
    int rankCount = getCountOfRank();
    auto treeReprAsString = std::string(argv[1]);
    int treeElementCount = treeReprAsString.length();

    // --- Handle edge cases ---
    //   If the tree is empty or has only one node or two nodes, its really doesn't make sense to compute this on multiple ranks.

    // If the tree has 0 or 1 node, it has no edges.
    if (treeElementCount < 1)
    {
        if (rank == 0)
        {
            std::fprintf(stderr, "Error: The binary tree is without edges. \n");
        }
        MPI_Finalize();
        return 10;
    }


    if (treeElementCount == 1)
    {
        if (rank == 0)
        {
                std::printf("%c:%d", treeReprAsString[0], 0);
        }
        MPI_Finalize();
        return 1;
    }

    // If the tree has exactly 2 nodes, it's a simple root-child case.
    if (treeElementCount == 2)
    {
        if (rank == 0)
            std::printf("%c:%d,%c:%d", treeReprAsString[0], 0, treeReprAsString[1], 1);
        MPI_Finalize();
        return 1;
    }


    // List of all directed edges (u, v) in the tree. Includes parent->child and child->parent.
    std::vector<Edge> edgeList;
    // Adjacency list representation using edge indices from edgeList.
    // edgeTreeReprezentation[i] stores indices of edges connected to vertex i.
    std::vector<std::vector<int>> edgeTreeReprezentation;
    /* Example representation:
      * edgeTreeReprezentation[0] = {0, 1} // Means vertex 0 connects via edgeList[0] and edgeList[1]
      */

    // its is faster to compute this than to gather that from all ranks
    for (int i = 0; i < treeElementCount; ++i)
    {
        int RchildIdx = (2 * i + 2 < treeElementCount) ? 2 * i + 2 : -1;
        int LchildIdx = (2 * i + 1 < treeElementCount) ? 2 * i + 1 : -1;
        int parentIdx = (i > 0) ? (i - 1) / 2 : -1;

        const int numOfChildrens = (LchildIdx != -1) + (RchildIdx != -1);

        // Temporary vector to store edge indices connected to the current node i.
        std::vector<int> tmp;


        // If node i has a parent, add the edge (i, parentIdx) and its reverse (parentIdx, i).
        if (parentIdx != -1)
        {
            edgeList.emplace_back(i, parentIdx);
            tmp.emplace_back(edgeList.size() - 1);
        }

        // If node i has a left child, add the edge (i, LchildIdx).
        if (LchildIdx != -1)
        {
            edgeList.emplace_back(i, LchildIdx);

            tmp.emplace_back(edgeList.size() - 1);
        }

        // If node i has a right child, add the edge (i, RchildIdx).
        if (RchildIdx != -1)
        {
            edgeList.emplace_back(i, RchildIdx);
            tmp.emplace_back(edgeList.size() - 1);
        }


        edgeTreeReprezentation.push_back(tmp);
    }


    // --- Verification ---
    // Check if the number of generated directed edges matches the number of MPI processes.
    // This algorithm assumes one process per directed edge.
    if (edgeList.size() != rankCount)
    {
        std::fprintf(stderr, "Edge list size: %d and number of proccess %d. This number should by equal! \n",
                     edgeList.size(), rankCount);
        return 2;
    }


    std::map<Edge, int> edgeToIndex;
    for (int i = 0; i < edgeList.size(); ++i)
    {
        edgeToIndex[edgeList[i]] = i;
    }

    // Create a vector to store the index of the reverse edge for each edge.
    // backEdgePartnerIndices[i] will store the index of the edge (v, u) if edgeList[i] is (u, v).
    std::vector<int> backEdgePartnerIndices(edgeList.size());
    for (int i = 0; i < edgeList.size(); ++i)
    {
        // Get the current edge (u, v).
        Edge currentEdge = edgeList[i];
        int u = currentEdge.first;
        int v = currentEdge.second;

        // Define the reverse edge (v, u).
        Edge reverseEdge = {v, u};

        // If found, store the index of the reverse edge.
        auto it = edgeToIndex.find(reverseEdge);
        if (it != edgeToIndex.end())
        {
            // Opačná hrana byla nalezena, její index je it->second
            backEdgePartnerIndices[i] = it->second;
        }
    }

    // ========== Real Parallel Section Starts Here ==========

    /** Euler Tour Successor Calculation */
    // Each process determines the successor edge for its assigned edge in an Euler tour.

    // Get the edge assigned to this rank.
    Edge myEdge = edgeList[rank];

    // Get the index of this edge.
    int myEdgeIndex = edgeToIndex[myEdge];

    // Define the reverse edge.
    Edge reverseEdge = {myEdge.second, myEdge.first};

    // Find the index of the reverse edge.
    auto it_rev = edgeToIndex.find(reverseEdge);
    int reverseEdgeIndex = it_rev->second;

    // Variable to store the index of the original successor edge in the Euler tour.
    int original_successor = -1;

    // The vertex from which the successor edge will depart.
    // For edge (u, v), the successor departs from v.
    int vertexNumber = reverseEdge.first;

    // Get the list of edge indices connected to vertex 'v'.
    auto& nextElements = edgeTreeReprezentation[vertexNumber];

    std::array<int, 2> myEtourResult = {myEdgeIndex, 0}; // second element is just place holder here

    // Find the position of the reverse edge (v, u) in the adjacency list of v.
    int reverseEdgePos = -1;
    for (int i = 0; i < nextElements.size(); ++i)
    {
        if (nextElements[i] == reverseEdgeIndex)
        {
            reverseEdgePos = i;
            break;
        }
    }

    if (reverseEdgePos != -1)
    {
        int successorPos = (reverseEdgePos + 1) % nextElements.size();
        int successorEdgeIndex = nextElements[successorPos];
        original_successor = nextElements[successorPos];

        myEtourResult[1] = successorEdgeIndex;
    }

    //====== end etour computing ======


    /** ===== Root Definition for List Ranking ======*/
    // We don't use any algorithm that works in every scenario.
    // We assume that the has atleast 3 nodes and is a binary tree.
    // then the last edge will be always the edge from node with id 1 to ROOT. -- THe edge case is handled at the start.

    auto rightChildRoot = edgeList[1];
    Edge rightChildRootReverse = {rightChildRoot.second, rightChildRoot.first};
    int breakEdgeIndex = edgeToIndex[rightChildRootReverse];

    /** ==== List Ranking (Wyllie’s Algorithm) ==== */
    // Determine the rank (position) of each edge in the linearized Euler tour.

    int N = edgeList.size(); // Total number of edges (and processes).

    int my_succ; // Stores the index of the *current* successor edge during pointer jumping.
    int my_dist; // Stores the *current* distance TO THE END of the list (breakEdgeIndex).


    my_succ = original_successor; // Use the Euler tour successor calculated earlier.
    my_dist = (rank == breakEdgeIndex) ? 0 : 1; // The break edge is distance 0, others start at 1.

    // Adjust the successor for the break edge process.
    // The process holding the break edge points to itself to terminate the list.
    if (rank == breakEdgeIndex)
    {
        my_succ = breakEdgeIndex;
    }

    // Calculate the number of steps needed for pointer jumping.
    int n_steps = std::ceil(std::log2(edgeList.size()));

    std::vector<int> all_dist(N);
    std::vector<int> all_succ(N);

    // Perform pointer jumping for n_steps iterations.
    for (int i = 0; i < n_steps; ++i)
    {
        MPI_Allgather(&my_dist, 1, MPI_INT, all_dist.data(), 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&my_succ, 1, MPI_INT, all_succ.data(), 1, MPI_INT, MPI_COMM_WORLD);

        // If this process is not already at the end of the list (breakEdgeIndex)...
        if (my_succ != breakEdgeIndex)
        {
            int current_succ_idx = my_succ;
            if (current_succ_idx >= 0 && current_succ_idx < N)
            {
                my_dist = my_dist + all_dist[current_succ_idx];
                my_succ = all_succ[current_succ_idx]; // jump
            }
            else
            {
                // Should not happen
                my_succ = breakEdgeIndex;
            }
        }

    }

    // Calculate the final position (0-based index from the start) of this process's edge.
    // Position = Total Edges - 1 - DistanceToEnd
    int finalPosn = N - 1 - my_dist;

    // Buffer to gather the final calculated position from all processes.
    std::vector<int> allFinalPosn;
    allFinalPosn.resize(N);

    MPI_Allgather(&finalPosn, 1, MPI_INT, allFinalPosn.data(), 1, MPI_INT, MPI_COMM_WORLD);

    /** === Level Calculation using Prefix/Suffix Sums on Euler Tour  === */

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
        my_initial_weight = -1; // Forward edge
    }
    else
    {
        my_initial_weight = +1; // Backward edge
    }


    std::vector<int> weights;
    weights.resize(N);

    MPI_Allgather(&my_initial_weight, 1, MPI_INT, weights.data(), 1, MPI_INT, MPI_COMM_WORLD);

    // Create a vector to store weights in the order determined by list ranking.
    std::vector<int> orderedWeights(N);
    for (int r = 0; r < N; ++r)
    {
        int pos = allFinalPosn[r];
        if (pos >= 0 && pos < N)
        {
            orderedWeights[pos] = weights[r];
        }
    }


    /** ==== Calculate Suffix Sums on the ordered weights ==== */
    // Each process now has the globally ordered weights and can compute the suffix sums locally.
    std::vector<int> suffixSumsOrdered(N);
    if (N > 0)
    {
        suffixSumsOrdered[N - 1] = orderedWeights[N - 1];
        for (int i = N - 2; i >= 0; --i)
        {
            // Iteruj od predposlednÃ©ho k prvÃ©mu
            suffixSumsOrdered[i] = orderedWeights[i] + suffixSumsOrdered[i + 1];
        }
    }

    // --- Find the relevant suffix sum for this process ---
    int myFinalPosn = allFinalPosn[rank];
    int my_final_suffix_sum;
    if (myFinalPosn >= 0 && myFinalPosn < N)
    {
        my_final_suffix_sum = suffixSumsOrdered[myFinalPosn];
    }

    // --- Calculate the final level for forward edges ---
    int my_level = -1;
    if (is_forward)
    {
        my_level = my_final_suffix_sum + 1;
    }

    // --- Result Gathering and Output ---
    // Each process holding a forward edge (u, v) now knows the level of node v.
    // Gather this information to construct the final level map for all nodes.

    std::vector<int> node_levels(treeElementCount, -1);
    int send_data[2] = {-1, -1};
    if (is_forward)
    {
        send_data[0] = v;
        send_data[1] = my_level;
    }


    // Buffer to receive all node/level pairs from all processes.
    std::vector<int> all_data(2 * N);
    MPI_Allgather(send_data, 2, MPI_INT, all_data.data(), 2, MPI_INT, MPI_COMM_WORLD);

    // --- Final Assembly ---
    // Set the root node's level to 0.
    node_levels[root_vertex_index] = 0;
    for (int i = 0; i < N; ++i)
    {
        int node_v = all_data[2 * i];
        int level_v = all_data[2 * i + 1];
        if (node_v != -1) node_levels[node_v] = level_v;
    }

    // --- Output ---
    // Rank 0 prints the final result in the specified format.
    if (rank == 0)
    {
        for (int i = 0; i < treeElementCount; ++i)
        {
            std::printf("%c:%d", treeReprAsString[i], node_levels[i]);
            if (i < treeElementCount - 1)
            {
                std::printf(",");
            }
        }
        std::printf("\n");
    }


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
