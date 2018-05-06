#include <iostream>
#include "SymmetricalTSP.h"

int main() {
	SymmetricalTSP tsp("berlin52.txt");
	std::cout << tsp.printAll();
	system("pause");
	return 0;
}