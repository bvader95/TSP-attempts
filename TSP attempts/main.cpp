#include <iostream>
#include "SymmetricalTSP.h"

int main() {
	SymmetricalTSP tsp("burma14.txt");
	TSPSolution sol = tsp.branchAndBound();
	std::cout << sol.value<<std::endl;
	for (unsigned int i = 0; i < sol.path.size(); ++i)std::cout << sol.path[i]<<" ";
	std::cout << std::endl;
	system("pause");
	return 0;
}