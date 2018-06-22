#pragma once
#include <utility> //pairs
#include <string>
#include <fstream> //reading from files
#include <iostream>
#include <cmath>   //square root f'n
#include <cfloat>  //max float values
#include <vector>
#include <stack>
#include <algorithm>


/**
* A structure representing a node of a tree created by the 
* Branch and Bound algorithm, containing vertices in order of visiting, 
* the length of the path, the information of what vertex was visited,
* and the node's lower bound
*/

struct BnBNode {
	std::vector<unsigned int> path;
	std::vector<bool> visited;
	double lowerBound;
	double value;
};

/**
 * A structure representing a solution for a TSP problem, 
 * containing vertices of graph, in order of visiting,
 * and the length of the path.
 */

struct TSPSolution {
	std::vector<unsigned int> path;
	double value;
};

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
	TSPSolution bruteForce();
	TSPSolution branchAndBound();
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
	std::string output = "";
	output = output + "Number of points: " + std::to_string(n) + "\n";
	output = output + "Coordinates:\n";
	output = output + "\tX\tY\n";
	for (unsigned int point = 0; point < n; ++point) {
		output.append(std::to_string(point + 1) + ":\t"
			+ std::to_string(coords[point].first) + "\t"
			+ std::to_string(coords[point].second) + "\n");
	}
	return output;
}

SymmetricalTSP::~SymmetricalTSP() {
	delete[] coords;
}

TSPSolution SymmetricalTSP::bruteForce() {
	std::vector<unsigned int> permutation;
	permutation.reserve(n);
	for (unsigned int i = 0; i < n; ++i) permutation.push_back(i);
	std::vector<unsigned int> initialState = permutation;
	double current;
	TSPSolution best;
	best.value = DBL_MAX;
	short printCounter = 0;
	do {
		std::next_permutation(permutation.begin() + 1, permutation.end());
		current = 0;
		for (unsigned int i = 0; i < n; ++i) {
			current = current + getDistance(permutation[i], permutation[(i + 1) % n]);
		}
		if (printCounter ==0  || best.value > current) {
			for (unsigned int i = 0; i < n; ++i) std::cout << permutation[i] << " ";
			std::cout << "\t" << current << "\t" << best.value << std::endl;
		}
		printCounter++;
		if (best.value > current) {
			best.path = permutation;
			best.value = current;
		}
	} while (permutation != initialState);
	return best;
}

TSPSolution SymmetricalTSP::branchAndBound() {
	//an array containing the length of the shortest and second shortest edge coming out of a vertex, 
	//used when calculating lower bounds
	std::vector<double> shortest;
	std::vector<double> secondShortest;
	shortest.resize(n, DBL_MAX);
	secondShortest.resize(n, DBL_MAX);
	double distance;
	for (unsigned int v1 = 0; v1 < n; ++v1) {
		for (unsigned int v2 = v1 + 1; v2 < n; ++v2) {	//it's symmetrical TSP, so E(V1, V2)=E(V2, V1), so we can do that 
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
	//which means an initial lower bound can be calculated as a half of a sum of two shortest edges coming out of each vertex
	//Calculating the initial lower bound
	double lowerBound = 0;
	for (unsigned int i = 0; i < n; ++i) {
		lowerBound = lowerBound + shortest[i] + secondShortest[i];
	}
	lowerBound = lowerBound / 2;
	//creating an initial root 'solution' to initialize the stack
	BnBNode best;
	best.lowerBound = lowerBound;
	best.path.reserve(n);//making the "path" vector just big enough to store the entire path
	best.visited.resize(n, false);//making the "visited" vector store n 'false's
	best.path.push_back(0);//adding vertex 0 as the first element of the path, which we can do w/o loss of generality
	best.visited[0] = true;
	best.value = UINT32_MAX;
	std::stack<BnBNode> s;//fun fact: if we used a queue here, the algorithm would be just a slower brute-force
	s.push(best);
	BnBNode current, child;
	short printCounter = 0;
	while (!s.empty()) {
		current = s.top();
		s.pop();
		if (printCounter == 0) {
			std::cout << current.lowerBound << "\t" << best.value << "\t";
			for (unsigned int i = 0; i < current.path.size(); ++i) std::cout << current.path[i] << " ";
			std::cout << std::endl;
		}
		++printCounter;
		if (current.path.size() == n) {//leaf! add the last edge, check if value is shorter
									   //than what we have so far, if so, save
			current.value = current.value + getDistance(current.path.back(), current.path.front());
			if (current.value < best.value) {
				best = current;
			}
		}
		else {//not a leaf! calculate the lower bound of current node's kids,
			  //and if it's not greater than what we've got, put it in the queue
			for (unsigned int v = 1; v < n; ++v) {
				if (current.visited[v] == false) {
					child = current;
					child.visited[v] = true;
					if (child.path.size() == 1) {//the way of calculating bound and value differs for level 2 nodes
						child.lowerBound = child.lowerBound
							- (shortest[child.path.back()] + shortest[v]) / 2
							+ getDistance(child.path.back(), v);
						child.value = getDistance(child.path.back(), v);
					}
					else {
						child.lowerBound = child.lowerBound
							- (secondShortest[child.path.back()] + shortest[v]) / 2
							+ getDistance(child.path.back(), v);
						child.value = child.value + getDistance(child.path.back(), v);
					}
					child.path.push_back(v);
					if (child.lowerBound <= best.value) s.push(child);
				}
			}
		}
	}
	TSPSolution output;
	output.path = best.path;
	output.value = best.value;
	return output;
}

double SymmetricalTSP::getDistance(unsigned int v1, unsigned int v2) {
	std::pair<double, double> *vOne = &coords[v1], *vTwo = &coords[v2];
	double x = (vTwo->first - vOne->first);
	double y = (vTwo->second - vOne->second);
	return sqrt(x*x + y*y);
}