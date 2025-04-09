#include <iostream>
#include <mpi.h>
#include <cstdint>
#include <vector>


#define early_break 1

int main(int argc, char* argv[])
{
    // Initialize the application

    MPI_Init(&argc, &argv); // Inicializace MPI

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); // Počet procesů
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // ID procesu

    // only rakn 0 is sending number to other ranks
    if (world_rank == 0)
    {
        // openfile
        auto file = fopen("numbers", "r");
        if (file == nullptr)
        {
            std::cerr << "Error opening file" << std::endl;
            return 1;
        }

        // load numbers and print numbers
        for (int i = 0; i < world_size; i++)
        {
            uint8_t num;
            // read one byte from file
            fread(&num, sizeof(uint8_t), 1, file);
            std::cout << (int)num << (i + 1 >= world_size ? "\n" : " ");
            MPI_Send(&num, 1, MPI_UINT8_T, i, 0, MPI_COMM_WORLD);
        }
    }

    uint8_t ranksNumber = 0;
    MPI_Recv(&ranksNumber, 1, MPI_UINT8_T, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    uint8_t finished = 0;

#if early_break
    bool done = false;
    int donecount = 0;
#endif

    int i = world_size;
    while (i-- > 0)
    {
        // if you want to break when the numbers are sorted
# if early_break
        if (done)
        {
            donecount++;
            done = false;
        }
        else
        {
            donecount = 0;
        }

        if (donecount > 1)
        {
            break;
        }
#endif

        // for not repeating code
        // for each throught mod result
        for (auto mod : {0, 1})
        {
            finished = 1; // default value

            uint8_t tmp_RanksNumber; // buffer for receiving neighbor data

            if (world_rank % 2 == mod && world_rank + 1 < world_size)
            {
                MPI_Send(&ranksNumber, 1, MPI_UINT8_T, world_rank + 1, 0, MPI_COMM_WORLD);
                MPI_Recv(&ranksNumber, 1, MPI_UINT8_T, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                finished = 1;
            }
            else if (world_rank % 2 != mod && world_rank > 0)
            {
                MPI_Recv(&tmp_RanksNumber, 1, MPI_UINT8_T, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if (tmp_RanksNumber > ranksNumber)
                {
                    MPI_Send(&ranksNumber, 1, MPI_UINT8_T, world_rank - 1, 0, MPI_COMM_WORLD);
                    ranksNumber = tmp_RanksNumber;
                    finished = 0; // not finished flag
                }
                else
                {
                    MPI_Send(&tmp_RanksNumber, 1, MPI_UINT8_T, world_rank - 1, 0, MPI_COMM_WORLD);
                    finished = 1; // finished flag
                }
            }

#if early_break
            unsigned global_done_int = 0;
            MPI_Allreduce(&finished, &global_done_int, 1, MPI_UINT8_T, MPI_SUM, MPI_COMM_WORLD);
            // sum all flags (0 or 1) to check if all processes are done


            if (global_done_int == world_size) // if all processes are done
            {
                done = true;
            }
#endif


            MPI_Barrier(MPI_COMM_WORLD);
        }
    }


    // collection buffer
    std::vector<uint8_t> recv_data;

    if (world_rank == 0)
    {
        recv_data.resize(world_size); // Rezervujeme místo pro všechna data
    }

    // Sběr dat na procesu 0
    MPI_Gather(&ranksNumber, 1, MPI_UINT8_T, recv_data.data(), 1, MPI_UINT8_T, 0, MPI_COMM_WORLD);


    for (auto i : recv_data)
    {
        std::cout << (int)i << std::endl;
    }


    MPI_Finalize(); // Ukončení MPI

    return 0;
}
