#pragma once
#include <utility> //pairs
#include <string>
#include <fstream> //reading from files
#include <iostream>
#include <cmath>   //square root f'n
#include <cfloat>  //max float values
#include <vector>

/*
* A class representing an instance of a travelling salesman problem
* where the distance between A and B is the same as between B and A
*/

struct BnBSolution {
	std::vector<unsigned int> vertices;//TODO: consider replacing with a priority queue
	std::vector<bool> visited;
	double lowerBound;
	double value;
};

class SymmetricalTSP {
private:
	unsigned int n;//the amount of points
	std::pair<double, double> *coords;//coordinates of the points
public:
	//the constructor loading an instance of a problem from a file
	SymmetricalTSP(std::string filename);
	~SymmetricalTSP();
	std::string printAll();
	BnBSolution branchAndBound();
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

BnBSolution SymmetricalTSP::branchAndBound() {
	//an array containing the length of the shortest and second shortest edge coming out of a vertice, 
	//used when calculating lower bounds
	double *shortest = new double[n], *secondShortest = new double[n];
	for (unsigned int i = 0; i < n; ++i) {
		shortest[i] = DBL_MAX;
		secondShortest[i] = DBL_MAX;
	}
	double distance;
	for (unsigned int v1 = 0; v1 < n; ++v1) {
		for (unsigned int v2 = v1 + 1; v2 < n; ++v2) {	//it's symmetrical TSP, so E(V1, V2)=E(V2, V1), 
			distance = getDistance(v1, v2);
			if (distance < shortest[v1]) {				//distance is shorter than the shortest so far
				secondShortest[v1] = shortest[v1];
				shortest[v1] = distance;
			}
			else if (distance < secondShortest[v1]) {	//distance is longer than the shortest so far,
														//but shorter than the second shortest
				secondShortest[v1] = distance;
			}
														//second verse, same as the first, but for V2
			if (distance < shortest[v2]) {				//distance is shorter than the shortest so far
				secondShortest[v2] = shortest[v2];
				shortest[v2] = distance;
			}
			else if (distance < secondShortest[v2]) {	//distance is longer than the shortest so far, 
														//but shorter than the second shortest
				secondShortest[v2] = distance;
			}
		}
	}
	//Note: cost of a tour T can be presented as a half of a sum of two edges adjacent to each vertex belonging to the tour,
	//which means an initial lower cound can be calculated as a half of a sum of two shortest edges coming out of each vertex
	//Calculating the initial lower bound
	double lowerBound = 0;
	for (unsigned int i = 0; i < n; ++i) {
		lowerBound = lowerBound + shortest[i] + secondShortest[i];
	}
	lowerBound = lowerBound / 2;
	//creating an 'solution' to initialize the queue
	BnBSolution best;
	best.lowerBound = lowerBound;
	best.vertices.reserve(n);//making the "vertices" vector just big enough to store the entire path
	best.visited.resize(n, false);//making the "visited" vector store 
	best.vertices.push_back(0);//adding vertex 0 as the first element of the path, which we can do w/o loss of generality
	best.visited[0] = true;
	best.value=UINT32_MAX;
	delete[] shortest;
	delete[] secondShortest;
	return best;
}

double SymmetricalTSP::getDistance(unsigned int v1, unsigned int v2) {
	std::pair<double, double> *vOne=&coords[v1], *vTwo=&coords[v2];
	double x = (vTwo->first - vOne->first);
	double y = (vTwo->second - vOne->second);
	return sqrt(x*x+y*y);
}