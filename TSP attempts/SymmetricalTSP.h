#pragma once
#include <utility>
#include <string>
#include <fstream>
#include <iostream>

/*
* A class representing an instance of a travelling salesman problem
* where the distance between A and B is the same as between B and A
*/
class SymmetricalTSP {
private:
	long long int n;//the amount of points
	std::pair<double, double> *coords;//coordinates of the points
public:
	//the constructor loading an instance of a problem from a file
	SymmetricalTSP(std::string filename);
	~SymmetricalTSP();
	std::string printAll();
};

SymmetricalTSP::SymmetricalTSP(std::string filename) {
	n = 0;
	std::ifstream input(filename);
	if (!input.is_open()) {
		std::cerr << "File couldn't be opened!" << std::endl;
		return;
	}
	input >> n;
	if (n <= 0) {
		std::cerr << "Malformed data (you can't give me a negative amount of points)" << std::endl;
		return;
	}
	coords = new std::pair<double, double>[n];
	//TODO: check for when there's not enough data and print an error
	for (int point = 0; point < n; ++point) {
		double x, y;
		input >> x >> y;
		coords[point] = std::pair<double, double>(x, y);
	}
	input.close();
}

std::string SymmetricalTSP::printAll() {
	if (n <= 0) {
		return "Class wasn't initialized";
	}
	std::string output="";
	output = output + "Number of points: "+std::to_string(n)+"\n";
	output = output + "Coordinates:\n";
	output = output + "\tX\tY\n";
	for (int point = 0; point < n; ++point) {
		output.append(std::to_string(point+1)+":\t"
			+std::to_string(coords[point].first)+"\t"
			+std::to_string(coords[point].second)+"\n");
	}
	return output;
}

SymmetricalTSP::~SymmetricalTSP(){
	delete[] coords;
}

