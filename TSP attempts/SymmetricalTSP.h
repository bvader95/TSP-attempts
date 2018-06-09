#pragma once
#include <utility> //pairs
#include <string>
#include <fstream> //reading from files
#include <iostream>
#include <cmath>   //square root f'n
#include <cfloat>  //max float values

/*
* A class representing an instance of a travelling salesman problem
* where the distance between A and B is the same as between B and A
*/
class SymmetricalTSP {
private:
	unsigned int n;//the amount of points
	std::pair<double, double> *coords;//coordinates of the points
public:
	//the constructor loading an instance of a problem from a file
	SymmetricalTSP(std::string filename);
	~SymmetricalTSP();
	std::string printAll();
	double branchAndBound();
	double getDistance(unsigned int v1, unsigned int v2);//an auxillary helper function for calculating distances between two vertices
};

SymmetricalTSP::SymmetricalTSP(std::string filename) {
	n = 0;
	std::ifstream input(filename);
	if (!input.is_open()) {
		std::cerr << "File couldn't be opened!" << std::endl;
		return;
	}
	input >> n;
	coords = new std::pair<double, double>[n];
	//TODO: check for when there's not enough data and print an error
	for (unsigned int point = 0; point < n; ++point) {
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
	for (unsigned int point = 0; point < n; ++point) {
		output.append(std::to_string(point+1)+":\t"
			+std::to_string(coords[point].first)+"\t"
			+std::to_string(coords[point].second)+"\n");
	}
	return output;
}

SymmetricalTSP::~SymmetricalTSP(){
	delete[] coords;
}

double SymmetricalTSP::branchAndBound() {
	//Note: cost of a tour T can be presented as a half of a sum of two edges adjacent to each vertex belonging to the tour,
	//which means an initial lower cound can be calculated as a half of a sum of two shortest edges coming out of each vertex
	//Calculating the initial lower bound
	double lowerBound = 0;
	double shortest1, shortest2;//shortest and second shortest edges found
	for (unsigned int v1 = 0; v1 < n; ++v1) {
		shortest1 = DBL_MAX;
		shortest2 = DBL_MAX;
		for (unsigned int v2 = 0; v2 < n; ++v2) {
			if (v1 == v2)continue;
			double distance = getDistance(v1, v2);
			if (distance < shortest1) {//the distance is the shortest edge found coming from a vertex v1
				shortest2 = shortest1;
				shortest1 = distance;
			}
			else if (distance < shortest2) {//the distance isn't the shortest edge found coming from v1, but is the second shortest
				shortest2 = distance;
			}
		}
		lowerBound = lowerBound + (shortest1 + shortest2) / 2;
	}
	return lowerBound;
}

double SymmetricalTSP::getDistance(unsigned int v1, unsigned int v2) {
	std::pair<double, double> *vOne=&coords[v1], *vTwo=&coords[v2];
	double x = (vTwo->first - vOne->first);
	double y = (vTwo->second - vOne->second);
	return sqrt(x*x+y*y);
}