#include <iostream>
#include "SymmetricalTSP.h"

int main() {
	SymmetricalTSP tsp("berlin52.txt");
	std::cout << tsp.branchAndBound().lowerBound<<std::endl;
	system("pause");
	return 0;
}