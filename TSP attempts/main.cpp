#include <iostream>
#include "PointTSP.h"

int main() {
	PointTSP tsp("westernsahara29.txt");
	TSPSolution sol = tsp.localSearch();
	std::cout << sol.value<<"\t";
	for (unsigned int i = 0; i < sol.path.size(); ++i)std::cout << sol.path[i]<<" ";
	std::cout << std::endl;
	system("pause");
	return 0;
}