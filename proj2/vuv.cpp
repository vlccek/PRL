#include <array>
#include <cmath>
#include <iostream>
#include <mpi.h>
#include <variant>
#include <vector>
#include <algorithm>

using Edge = std::pair<int, int>;

int getNodeLevel(int nodeIndex)
{
    return static_cast<int>(std::floor(std::log2(static_cast<double>(nodeIndex) + 1.0)));
}

/**
 * @brief Spočítá skutečný počet dětí daného uzlu v binárním stromě. Na základě indexu a celkového počtu uzlů.
 *
 * @param parentIndex Index uzlu  (předpokládáme >= 0 a < N).
 * @param N Celkový počet uzlů ve stromě (předpokládáme > 0).
 * @return int Počet dětí (0, 1 nebo 2).
 */
int countChildren(int parentIndex, int N)
{
    int LchildIdx = 2 * parentIndex + 1;
    int RchildIdx = 2 * parentIndex + 2;

    return (LchildIdx < N) + (RchildIdx < N);
}

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

std::vector<std::pair<int, int>> eulerPath(int N, int id = 0)
{
    std::vector<std::pair<int, int>> result;
    int RchildIdx = (2 * id + 2 < N) ? 2 * id + 2 : -1;
    int LchildIdx = (2 * id + 1 < N) ? 2 * id + 1 : -1;

    if (LchildIdx != -1)
    {
        result.emplace_back(id, LchildIdx);

        auto tmpres = eulerPath(N, LchildIdx);
        result.insert(result.end(), tmpres.begin(), tmpres.end());
        result.emplace_back(LchildIdx, id);
    }

    if (RchildIdx != -1)
    {
        result.emplace_back(id, RchildIdx);

        auto tmpres = eulerPath(N, RchildIdx);
        result.insert(result.end(), tmpres.begin(), tmpres.end());
        result.emplace_back(RchildIdx, id);
    }

    return result;
}


int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <binary_tree>" << std::endl;
        std::cerr << "Please provide a binary tree as a command line argument." << std::endl;
        std::cerr << "Example: " << argv[0] << "148" << std::endl;
        return 1;
    }
    const int rank = getRank();


    auto treeString = std::string(argv[1]);
    int elementCount = treeString.length();
    char nodeElement = argv[1][rank];

    const int numProcesses = getCountOfRank();
    const int maxLevel = getNodeLevel(elementCount);
    const int myLevel = getNodeLevel(rank);
    const int idInLevel = rank % (myLevel + 1);

    const int parentId = (rank - 1) / 2;

    // if children exists safe the childer id otherwise 0
    int RchildIdx = (2 * rank + 2 < numProcesses) ? 2 * rank + 2 : -1;
    int LchildIdx = (2 * rank + 1 < numProcesses) ? 2 * rank + 1 : -1;

    const int numOfChildrens = (LchildIdx < numProcesses) + (RchildIdx < numProcesses);


    // etour


    std::vector<std::pair<int, int>> localEulerPath = {
        {parentId, rank},
        {LchildIdx, rank},
        {RchildIdx, rank},
    };
    auto tmp = std::remove_if(localEulerPath.begin(), localEulerPath.end(), [](const auto& edge)
    {
        return edge.first == -1 || edge.second == -1;
    }); // remove invalid edges

    std::vector<std::pair<int, int>> eulerPath;
    eulerPath.reserve(elementCount * 2); // preallocate space for the euler path


    MPI_Allgather(localEulerPath.data(), localEulerPath.size() * sizeof(std::pair<int, int>), MPI_BYTE,
                  localEulerPath.data(), eulerPath.size() * sizeof(std::pair<int, int>), MPI_BYTE, MPI_COMM_WORLD);


    return 1;


    MPI_Finalize();
    return 0;
    // TIP See CLion help at <a href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>. Also, you can try interactive lessons for CLion by selecting 'Help | Learn IDE Features' from the main menu.
}
