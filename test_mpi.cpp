#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>


int main(int argc, char *argv[])
{
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_File outFile;
	MPI_File_open(MPI_COMM_WORLD, "test.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &outFile);

	std::stringstream ss;
	ss << "Hello" << rank << std::endl;
	std::string str = ss.str();
	int strLen = str.length();
	MPI_File_write_at_all(outFile, rank * strLen, str.c_str(), strLen, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_write_at_all(outFile, (rank +size )* strLen, str.c_str(), strLen, MPI_CHAR, MPI_STATUS_IGNORE);

	// // write eof
	// MPI_File_write_at_all(outFile, size * str.size(), "", 0, MPI_CHAR, MPI_STATUS_IGNORE);

	MPI_File_close(&outFile);

	MPI_Finalize();
	return 0;
}