#pragma once
#include "AuxiliaryTSPh.h"	//some structures I wrote
#include <utility>			//pairs
#include <string>			//filenames for fstream, among other things
#include <fstream>			//reading from files
#include <iostream>			//printing on screen to show the progress of the current f'n
#include <cmath>			//square root f'n
#include <cfloat>			//max float values
#include <stack>			//Branch and Bound; travelling through the solution tree
#include <algorithm>		//fiddling with arrays
#include <chrono>			//time stuff
#include <random>			//std::default_random_engine()


/*
* A class representing an instance of a travelling salesman problem
* that takes coordinates of the cities as argument (as opposed to distances between them)
*/

class PointTSP {
private:
	//the main fields of the class, storing the params of the problem
	unsigned int n;//the amount of points
	std::pair<double, double> *coords;//coordinates of the points
	//auxilliary helper functions that have no business being public
	double getDistance(unsigned int v1, unsigned int v2);
	double calculatePathsLength(TSPSolution &sol);
	TSPSolution singleTwoOpt(TSPSolution sol, unsigned int begin, unsigned int end);
	TSPSolution swapTwoVertices(TSPSolution sol, unsigned int v1, unsigned int v2);
	TSPSolution generateRandomSolution();
	TSPSolution generateNearestNeighbourSolution();
	TSPSolution generateNearestNeighbourSolution(unsigned int initialVertex);
public:
	//the constructor loading an instance of a problem from a file
	PointTSP(std::string filename);
	~PointTSP();
	std::string printAll();
	TSPSolution bruteForce(bool showProgress);
	TSPSolution bruteForce();
	TSPSolution branchAndBound(bool showProgress);
	TSPSolution branchAndBound();
	TSPSolution localSearch();
};

PointTSP::PointTSP(std::string filename) {
	n = 0;
	std::ifstream input(filename);
	if (!input.is_open()) {
		std::cerr << "File couldn't be opened!" << std::endl;
		return;
	}
	input >> n;
	coords = new std::pair<double, double>[n];
	unsigned int readVariables = 0;
	double buffer;
	std::pair<double, double> bufferPair;
	while (input>>buffer) {
		if (readVariables >= 2 * n)break;//stopping after having read enough data
		if (readVariables % 2 == 0) {
			bufferPair.first = buffer;
		}
		else {
			bufferPair.second = buffer;
			coords[readVariables/2] = bufferPair;
		}
		readVariables++;
	}
	if (readVariables < 2 * n) {
		delete[] coords;
		n = 0;
		std::cerr << "File didn't contain enough variables, check the validity of the data!"<<std::endl;
	}
	input.close();
}

std::string PointTSP::printAll() {
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

PointTSP::~PointTSP() {
	delete[] coords;
}

TSPSolution PointTSP::bruteForce() {
	return bruteForce(false)
}
TSPSolution PointTSP::bruteForce(bool showProgress) {
	TSPSolution current, best;
	current.path.reserve(n);
	for (unsigned int i = 0; i < n; ++i) current.path.push_back(i);
	std::vector<unsigned int> initialState = current.path;
	best.value = DBL_MAX;
	short printCounter = 0;
	do {
		std::next_permutation(current.path.begin() + 1, current.path.end());
		calculatePathsLength(current);
		if ((printCounter == 0 || best.value > current.value) && showProgress) {
			for (unsigned int i = 0; i < n; ++i) std::cout << current.path[i] << " ";
			std::cout << "\t" << current.value << "\t" << best.value << std::endl;
		}
		printCounter++;
		if (best.value > current.value) {
			best = current;
		}
	} while (current.path != initialState);
	return best;
}

TSPSolution PointTSP::branchAndBound() {
	return branchAndBound(false);
}
TSPSolution PointTSP::branchAndBound(bool showProgress) {
	//an array containing the length of the shortest and second shortest edge coming out of a vertex, 
	//used when calculating lower bounds
	std::vector<double> shortest;
	std::vector<double> secondShortest;
	shortest.resize(n, DBL_MAX);
	secondShortest.resize(n, DBL_MAX);
	double distance;
	for (unsigned int v1 = 0; v1 < n; ++v1) {
		for (unsigned int v2 = 0; v2 < n; ++v2) { 
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

	//Note: cost of a tour T can be presented as a half of a sum 
	//of two edges adjacent to each vertex belonging to the tour,
	//which means an initial lower bound can be calculated as a half
	//of a sum of two shortest edges coming out of each vertex
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
		if (printCounter == 0 && showProgress) {
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

double PointTSP::getDistance(unsigned int v1, unsigned int v2) {
	std::pair<double, double> *vOne = &coords[v1], *vTwo = &coords[v2];
	double x = (vTwo->first - vOne->first);
	double y = (vTwo->second - vOne->second);
	return sqrt(x*x + y*y);
}

double PointTSP::calculatePathsLength(TSPSolution &sol) {
	double length = 0.0;
	for (unsigned int i = 0; i < n; ++i)length = length + getDistance(sol.path[i], sol.path[(i + 1)%n]);
	sol.value = length;
	return length;
}


TSPSolution PointTSP::localSearch() {
	//generating the initial solution
	TSPSolution bestOverall = generateNearestNeighbourSolution();
	//TSPSolution bestOverall = generateRandomSolution();
	TSPSolution bestInCycle = bestOverall, neighbour;
	//finding the best result that the chosen neighbour-finding solution can produce
	while (true) {
		for (unsigned int v1 = 0; v1 < n; ++v1) {
			for (unsigned int v2 = 0; v2 < n; ++v2) {
				if (v1==v2) continue;//one-vertex-long "path"
				neighbour = singleTwoOpt(bestOverall, v1, v2+1);//see note for singleTwoOpt for explanation why there's a +1
				//neighbour = swapTwoVertices(bestOverall, v1, v2);
				if (neighbour.value < bestInCycle.value) { 
					bestInCycle = neighbour;
				}
			}
		}
		if (bestOverall.value == bestInCycle.value) break;
		else {
			bestOverall = bestInCycle;
		}
		std::cout << bestInCycle.value<<"\t";
		for (unsigned int i = 0; i < n; ++i) std::cout << bestInCycle.path[i] << " ";
		std::cout << std::endl;
	}
	return bestOverall;
}


//IMPORTANT NOTE:
//This function does the 2-opt for the range [begin, end),
//because that's how the functions from <algorithm> operate
//and it looks better than adding +1 everywhere
TSPSolution PointTSP::singleTwoOpt(TSPSolution sol, unsigned int begin, unsigned int end) {
	TSPSolution output=sol;
	if (begin > end) {
		std::rotate(output.path.begin(), output.path.begin()+begin, output.path.end());
		std::reverse(output.path.begin(), output.path.begin() + n - begin + end);
	}
	else {
		std::reverse(output.path.begin() + begin, output.path.begin() + end);
	}
	calculatePathsLength(output);
	return output;
}

TSPSolution PointTSP::swapTwoVertices(TSPSolution sol, unsigned int v1, unsigned int v2) {
	TSPSolution output = sol;
	std::swap(output.path[v1], output.path[v2]);
	calculatePathsLength(output);
	return output;
};

inline TSPSolution PointTSP::generateRandomSolution(){
	TSPSolution output;
	output.path.reserve(n);
	for (unsigned int i = 0; i < n; ++i)output.path.push_back(i);
	std::shuffle(output.path.begin(), output.path.end(),
		std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
	calculatePathsLength(output);
	return output;
}

inline TSPSolution PointTSP::generateNearestNeighbourSolution() {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::uniform_int_distribution<unsigned int> dist(0, n - 1);
	return generateNearestNeighbourSolution(dist(std::default_random_engine(seed)));
}

inline TSPSolution PointTSP::generateNearestNeighbourSolution(unsigned int initialVertex){
	TSPSolution output;
	std::vector<bool> visited;//visited[x] is true if vertex x was visited when generating, false otherwise
	output.path.push_back(initialVertex);
	output.value = 0;
	visited.resize(n, false);
	visited[initialVertex] = true;
	unsigned int currentVertex, nextVertex;
	double shortestEdge;
	while (output.path.size()<n){
		//initialization
		currentVertex = output.path.back();
		shortestEdge = DBL_MAX;
		//find the shortest path from the current vertex that wouldn't create a cycle
		for (unsigned int i = 0; i < n; ++i) {
			if (currentVertex == i)continue;
			if (!visited[i] && getDistance(currentVertex, i) < shortestEdge) {
				shortestEdge = getDistance(currentVertex, i);
				nextVertex = i;
			}
		}
		//add that path to the length and mark the relevant vertex as visited
		output.value = output.value + shortestEdge;
		output.path.push_back(nextVertex);
		visited[nextVertex] = true;
	}
	//close the cycle
	output.value = output.value + getDistance(output.path.back(), output.path.front());
	return output;
}
