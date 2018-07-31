#include <mpi.h>
#include <iostream>

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0)
    {
        int val = 17;
        int result = MPI_Send(&val, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        if (result == MPI_SUCCESS)
            std::cout << "rank 0 ok" << std::endl;
    }
    else if (rank == 1)
    {
        int val;
        int result = MPI_Recv(&val, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (result == MPI_SUCCESS && val == 17)
            std::cout << "rank 1 ok" << std::endl;

    }
    MPI_Finalize();
    return 0;
}
