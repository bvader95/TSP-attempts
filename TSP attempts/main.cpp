#include <iostream>
#include "PointTSP.h"

int main() {
	PointTSP tsp("burma14.txt");
	auto start = std::chrono::steady_clock::now();
	TSPSolution sol = tsp.branchAndBound();
	auto stop = std::chrono::steady_clock::now();
	auto timeInMilliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
	std::cout << sol.value<<"\t";
	for (unsigned int i = 0; i < sol.path.size(); ++i)std::cout << sol.path[i]<<" ";
	std::cout << std::endl;
	std::cout << "Time in ms: " << timeInMilliseconds << std::endl;
	system("pause");
	return 0;
}