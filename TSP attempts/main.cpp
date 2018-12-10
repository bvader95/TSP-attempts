#include <iostream>
#include "MatrixTSP.h"

int main() {

	//const unsigned int amountOfTests = 20;

	//std::fstream outputFile;
	//outputFile.open("test_output.txt", std::fstream::out);
	//std::string filenames[] = { "matrix_6.txt", "matrix_10.txt", "matrix_12.txt", "matrix_13.txt", "matrix_14.txt" };

	//for (unsigned int file = 0; file < 4; ++file) {
	//	MatrixTSP tsp(filenames[file]);
	//	TSPSolution sol;
	//	outputFile << "TEST CASE: " << filenames[file]<<"\n";
	//	std::cout << "TEST CASE: " << filenames[file]<<"\n";
	//	outputFile << "BrtFrc\tResult\tBnB\tResult\n";
	//	for (unsigned int test = 0; test < amountOfTests; ++test) {
	//		std::cout << "TEST "<<test+1<<"\n";
	//		auto start = std::chrono::steady_clock::now();
	//		sol = tsp.bruteForce(false);
	//		auto stop = std::chrono::steady_clock::now();
	//		outputFile << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()<<"\t";
	//		outputFile << sol.value<<"\t";
	//		start = std::chrono::steady_clock::now();
	//		sol = tsp.branchAndBound(false);
	//		stop = std::chrono::steady_clock::now();
	//		outputFile << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << "\t";
	//		outputFile << sol.value << "\n";
	//	}
	//}
	//system("pause");
	//return 0;
	
	std::string filename;
	int algorithmNo;
	std::cout << "Filename? ";
	std::cin >> filename;
	std::cout << "Algorithm? (0 - brute force, 1 - BnB, 2 - annealing) ";
	std::cin >> algorithmNo;
	MatrixTSP tsp(filename);
	TSPSolution sol;
	std::chrono::time_point<std::chrono::steady_clock> start, stop;
	long long timeInMilliseconds;
	switch (algorithmNo) {
	case 0:
		start = std::chrono::steady_clock::now();
		sol = tsp.bruteForce(true);
		stop = std::chrono::steady_clock::now();
		break;
	case 1:
		start = std::chrono::steady_clock::now();
		sol = tsp.branchAndBound(false);
		stop = std::chrono::steady_clock::now();
		break;
	case 2:
		start = std::chrono::steady_clock::now();
		sol = tsp.simulatedAnnealing(100000, 1, 0.99);
		stop = std::chrono::steady_clock::now();
		break;
	default:
		std::cout << "oops, wrong number, shutting down..."<<std::endl;
		break;
	}
	timeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
	std::cout << sol.value<<"\t";
	for (unsigned int i = 0; i < sol.path.size(); ++i)std::cout << sol.path[i]<<" ";
	std::cout << std::endl;
	std::cout << "Time in ms: " << timeInMilliseconds << std::endl;
	system("pause");
	return 0;
}