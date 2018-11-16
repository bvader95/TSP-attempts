#include <iostream>
#include "MatrixTSP.h"

int main() {

	const unsigned int amountOfTests = 20;

	std::fstream outputFile;
	outputFile.open("test_output.txt", std::fstream::out);
	std::string filenames[] = { "matrix_6.txt", "matrix_10.txt", "matrix_12.txt", "matrix_13.txt", "matrix_14.txt" };

	for (unsigned int file = 0; file < 4; ++file) {
		MatrixTSP tsp(filenames[file]);
		TSPSolution sol;
		outputFile << "TEST CASE: " << filenames[file]<<"\n";
		std::cout << "TEST CASE: " << filenames[file]<<"\n";
		outputFile << "BrtFrc\tResult\tBnB\tResult\n";
		for (unsigned int test = 0; test < amountOfTests; ++test) {
			std::cout << "TEST "<<test+1<<"\n";
			auto start = std::chrono::steady_clock::now();
			sol = tsp.bruteForce(false);
			auto stop = std::chrono::steady_clock::now();
			outputFile << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()<<"\t";
			outputFile << sol.value<<"\t";
			start = std::chrono::steady_clock::now();
			sol = tsp.branchAndBound(false);
			stop = std::chrono::steady_clock::now();
			outputFile << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "\t";
			outputFile << sol.value << "\n";
		}
	}
	system("pause");
	return 0;

	//MatrixTSP tsp("matrix_14.txt");
	//auto start = std::chrono::steady_clock::now();
	//TSPSolution sol = tsp.branchAndBound(false);
	//auto stop = std::chrono::steady_clock::now();
	//auto timeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
	//std::cout << sol.value<<"\t";
	//for (unsigned int i = 0; i < sol.path.size(); ++i)std::cout << sol.path[i]<<" ";
	//std::cout << std::endl;
	//std::cout << "Time in ms: " << timeInMilliseconds << std::endl;
	//system("pause");
	//return 0;
}