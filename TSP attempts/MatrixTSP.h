#pragma once
#include "AuxiliaryTSPh.h"	//some structures I wrote
#include <string>			//filenames for fstream, among other things
#include <fstream>			//reading from files
#include <iostream>			//printing on screen to show the progress of the current f'n
//#include <cmath>			//square root f'n
#include <cfloat>			//max float values
#include <stack>			//Branch and Bound; travelling through the solution tree
#include <algorithm>		//fiddling with arrays
#include <chrono>			//time stuff
#include <random>			//std::default_random_engine()

/*
* A class representing an instance of a travelling salesman problem
* that takes coordinates of the cities as argument (as opposed to distances between them)
*/

class MatrixTSP {
private:
	//the main fields of the class, storing the params of the problem
	unsigned int n;//the amount of points
	std::vector<std::vector<double>> distances; //distance matrix
	//auxilliary helper functions that have no business being public
	double getDistance(unsigned int v1, unsigned int v2);
	double calculatePathsLength(TSPSolution &sol);
	TSPSolution singleTwoOpt(TSPSolution sol, unsigned int begin, unsigned int end);
	TSPSolution swapTwoVertices(TSPSolution sol, unsigned int v1, unsigned int v2);
	TSPSolution generateRandomSolution();
	TSPSolution generateNearestNeighbourSolution();
	TSPSolution generateNearestNeighbourSolution(unsigned int initialVertex);
	TSPSolution generateRandomNeighbour(TSPSolution solution);
	TSPSolution orderCrossover(TSPSolution parent1, TSPSolution parent2);
	TSPSolution cycleCrossover(TSPSolution parent1, TSPSolution parent2);
	void cutoffSelection(std::vector<TSPSolution> &population, unsigned int initialPopulation);
	void rouletteWheelSelection(std::vector<TSPSolution> &population, unsigned int initialPopulation);
public:
	//the constructor loading an instance of a problem from a file
	MatrixTSP(std::string filename);
	std::string printAll();
	TSPSolution bruteForce(bool showProgress);
	TSPSolution bruteForce();
	TSPSolution branchAndBound(bool showProgress);
	TSPSolution branchAndBound();
	TSPSolution localSearch();
	TSPSolution simulatedAnnealing(double coolingCoefficient);
	TSPSolution geneticAlgorithm(unsigned int initialPopulation, char crossoverOp, char selectionMethod);
};

MatrixTSP::MatrixTSP(std::string filename) {
	n = 0;
	std::ifstream input(filename);
	if (!input.is_open()) {
		std::cerr << "File couldn't be opened!" << std::endl;
		return;
	}
	input >> n;
	//last-minute changes, I feel like a student again already!
	distances.resize(n);
	for (unsigned int i = 0; i < n; ++i)distances[i].resize(n);
	unsigned int x,y;//variables marking the spot where the read value will be placed in the matrix
	x = 0, y = 0;
	double buffer;
	//reading from the file,
	//this probably could've been a for loop, but y'know...
	while (input>>buffer) {
		distances[x][y] = buffer;
		++y;
		if (y == n) {//row finished, move to the next one
			y = 0;
			x++;
		}
		if (x == n) {//all rows filled
			break;
		}
	}
	if (x!=n) {//not all rows were filled
		n = 0;
		std::cerr << "File didn't contain enough variables, check the validity of the data!" << std::endl;
	}
	input.close();
}

std::string MatrixTSP::printAll() {
	if (n <= 0) {
		return "Class wasn't initialized";
	}
	std::string output = "";
	output = output + "Number of points: " + std::to_string(n) + "\n";
	output = output + "Coordinates:\n";
	output = output + "\tX\tY\n";
	for (unsigned int row = 0; row < n; ++row) {
		for (unsigned int col = 0; col < n; ++col) {
			output = output + std::to_string(distances[row][col]) + "\t";
		}
		output = output + "\n";
	}
	return output;
}

TSPSolution  MatrixTSP::bruteForce() {
	return bruteForce(false);
}

TSPSolution MatrixTSP::bruteForce(bool showProgress) {
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

TSPSolution MatrixTSP::branchAndBound() {
	return branchAndBound(false);
}
TSPSolution MatrixTSP::branchAndBound(bool showProgress) {
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

double MatrixTSP::getDistance(unsigned int v1, unsigned int v2) {
	return distances[v1][v2];
}

double MatrixTSP::calculatePathsLength(TSPSolution &sol) {
	double length = 0.0;
	for (unsigned int i = 0; i < n; ++i)length = length + getDistance(sol.path[i], sol.path[(i + 1) % n]);
	sol.value = length;
	return length;
}


TSPSolution MatrixTSP::localSearch() {
	//generating the initial solution
	TSPSolution bestOverall = generateNearestNeighbourSolution();
	//TSPSolution bestOverall = generateRandomSolution();
	TSPSolution bestInCycle = bestOverall, neighbour;
	//finding the best result that the chosen neighbour-finding solution can produce
	while (true) {
		for (unsigned int v1 = 0; v1 < n; ++v1) {
			for (unsigned int v2 = 0; v2 < n; ++v2) {
				if (v1 == v2) continue;//one-vertex-long "path"
				neighbour = singleTwoOpt(bestOverall, v1, v2 + 1);//see note for singleTwoOpt for explanation why there's a +1
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
		//std::cout << bestInCycle.value << "\t";
		//for (unsigned int i = 0; i < n; ++i) std::cout << bestInCycle.path[i] << " ";
		//std::cout << std::endl;
	}
	std::sort(bestOverall.path.begin(), bestOverall.path.end());
	return bestOverall;
}

inline TSPSolution MatrixTSP::simulatedAnnealing(double coolingCoefficient) {
	//the <random> stuff that will be used to decide if we're accepting a worse solution or not
	std::uniform_real_distribution<double> dist(0, 1);
	TSPSolution best = generateNearestNeighbourSolution(), bestest = best, current;
	double chanceToAccept, diceRoll;
	int endCounter = 0;
	//the initial temperature is an absolute value of a difference between two randomly-generated solutions, divided by log(0.7)
	//so that the initial chance to accept a solution is around 70%
	double temperature = abs((generateNearestNeighbourSolution().value-best.value)/log(0.7));
	while (true) {
		current = generateRandomNeighbour(best);
		if (best.value >= current.value) {//we found a better solution, accept it unconditionally
			best = current;
			if (bestest.value>=current.value) bestest = current;
		}
		else {//the solution isn't better, but we might accept it based on "temperature"
			chanceToAccept=exp((best.value-current.value)/temperature);
			diceRoll = dist(std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));//yeah, it's a weird dice
			if (chanceToAccept < 0.00001) {
				endCounter++;
				if (endCounter>n) break;//ending when the chance of accepting a worse solution falls below a specific threshold a certain amount of times
			}
			if (chanceToAccept > diceRoll)best = current;
		}
		temperature *= coolingCoefficient;
		//std::cout <<"T= "<< temperature << std::endl;
	}
	if (best.value > bestest.value) return bestest;
	return best;
}


//IMPORTANT NOTE:
//This function does the 2-opt for the range [begin, end),
//because that's how the functions from <algorithm> operate
//and it looks better than adding +1 everywhere
TSPSolution MatrixTSP::singleTwoOpt(TSPSolution sol, unsigned int begin, unsigned int end) {
	TSPSolution output = sol;
	if (begin > end) {
		std::rotate(output.path.begin(), output.path.begin() + begin, output.path.end());
		std::reverse(output.path.begin(), output.path.begin() + n - begin + end);
	}
	else {
		std::reverse(output.path.begin() + begin, output.path.begin() + end);
	}
	calculatePathsLength(output);
	return output;
}

TSPSolution MatrixTSP::swapTwoVertices(TSPSolution sol, unsigned int v1, unsigned int v2) {
	TSPSolution output = sol;
	std::swap(output.path[v1], output.path[v2]);
	calculatePathsLength(output);
	return output;
};

TSPSolution MatrixTSP::generateRandomSolution() {
	TSPSolution output;
	output.path.reserve(n);
	for (unsigned int i = 0; i < n; ++i)output.path.push_back(i);
	std::shuffle(output.path.begin(), output.path.end(),
		std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
	calculatePathsLength(output);
	return output;
}

TSPSolution MatrixTSP::generateNearestNeighbourSolution() {
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::uniform_int_distribution<unsigned int> dist(0, n - 1);
	return generateNearestNeighbourSolution(dist(std::default_random_engine(seed)));
}

TSPSolution MatrixTSP::generateNearestNeighbourSolution(unsigned int initialVertex) {
	TSPSolution output;
	std::vector<bool> visited;//visited[x] is true if vertex x was visited when generating, false otherwise
	output.path.push_back(initialVertex);
	output.value = 0;
	visited.resize(n, false);
	visited[initialVertex] = true;
	unsigned int currentVertex, nextVertex;
	double shortestEdge;
	while (output.path.size()<n) {
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


TSPSolution MatrixTSP::generateRandomNeighbour(TSPSolution solution) {
	std::uniform_int_distribution<unsigned int> dist(0, n - 1);
	unsigned int begin, end;
	begin = dist(std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
	end = dist(std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
	return singleTwoOpt(solution, begin, end);
}

TSPSolution MatrixTSP::orderCrossover(TSPSolution parent1, TSPSolution parent2){
	//creating a vector storing data on what vertices are present in the child
	std::vector<bool> presentInChild;
	presentInChild.resize(n, false);
	//picking two different random indices 
	std::uniform_int_distribution<unsigned int> dist(0, n - 1);
	unsigned int begin, end;
	begin = dist(std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
	do {
		end = dist(std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));
	} while (begin == end);
	//swappping them if "begin" comes before "end"
	if (begin > end)std::swap(begin, end);
	//copy the bit between the "begin" and "end" indices in the parent1 path to the child 
	TSPSolution child;
	child.path.reserve(n);
	for (unsigned int i = begin; i <= end; ++i) {
		child.path.push_back(parent1.path[i]);
		presentInChild[parent1.path[i]] = true;
	}
	//copy all the other vertices from parent2, in the order they are in it
	for (unsigned int i = 0; i < n; ++i) {
		if (!presentInChild[parent2.path[i]])child.path.push_back(parent2.path[i]);
	}
	calculatePathsLength(child);
	return child;
}

TSPSolution MatrixTSP::cycleCrossover(TSPSolution parent1, TSPSolution parent2)
{
	//creating a vector storing data on what vertices were used to make a cycle
	std::vector<bool> presentInCycle;
	presentInCycle.resize(n, false);
	//filling up an array with positions of cities in parent1
	//p1pos[x] means that the city X is present at parent1.path[p1pos[x]]
	std::vector<unsigned int> p1pos;
	p1pos.resize(n);
	for (int i = 0; i < n; ++i) {
		p1pos[parent1.path[i]] = i;
	}
	//the vector of vectors, each containing a CX cycle
	std::vector<std::vector<unsigned int>> cycles;
	//creating the cycles
	for (int i = 0; i < n; ++i) {
		if (presentInCycle[i] == true) continue;
		std::vector<unsigned int> newCycle;
		newCycle.push_back(i);
		presentInCycle[i] = true;
		//the next position in the CX cycle is as follows:
		//get the value V from parent2 from the same position as the last position in cycle
		//then find the position of V in parent1
		unsigned int nextInCycle = p1pos[parent2.path[i]];
		while (nextInCycle != newCycle[0]){
			newCycle.push_back(nextInCycle);
			presentInCycle[nextInCycle] = true;
			nextInCycle = p1pos[parent2.path[nextInCycle]];
		}
		cycles.push_back(newCycle);
	}
	TSPSolution child;
	child.path.resize(n, 999999999);
	//copying values from parents, alternating based on cycles
	for (unsigned int cycleNo = 0; cycleNo < cycles.size(); ++cycleNo) {
		std::vector<unsigned int> currentCycle = cycles[cycleNo];
		for (unsigned int i = 0; i < currentCycle.size(); ++i) {
			if (cycleNo % 2 == 0) {
				child.path[currentCycle[i]]=parent1.path[currentCycle[i]];
			}
			else {
				child.path[currentCycle[i]] = parent2.path[currentCycle[i]];
			}
		}
	}
	calculatePathsLength(child);
	return child;
}

void MatrixTSP::cutoffSelection(std::vector<TSPSolution>& population, unsigned int initialPopulation){
	//sort the population table
	std::sort(population.begin(), population.end());
	//keep the N best solutions
	population.resize(initialPopulation);
}

void MatrixTSP::rouletteWheelSelection(std::vector<TSPSolution>& population, unsigned int initialPopulation){ 
	std::sort(population.begin(), population.end());
	double bestSol = population[0].value, worstSol = population[population.size()-1].value;
	//selection time! 
	//the <random> stuff that will be used to decide if we're accepting a solution or not
	std::uniform_real_distribution<double> dist(0, 1);
	for (int i = 0; i < population.size(); ++i) {
		double diceRoll = dist(std::default_random_engine(std::chrono::system_clock::now().time_since_epoch().count()));//yeah, it's a weird dice
		double fitFunctionResult = 2*(population[i].value-bestSol)/(bestSol);
		if (diceRoll < fitFunctionResult) {//the solution was deemed unfit, REMOVE IT! and adjust i to take that into account
			population.erase(population.begin()+i);
			--i;
		}
	}
	//cutting down the population if it's too big;
	if (population.size() > initialPopulation << 8 ) {
		population.resize(initialPopulation<<8);
	}
}

TSPSolution MatrixTSP::geneticAlgorithm(unsigned int initialPopulation, char crossoverOp, char selectionMethod){
	std::vector<TSPSolution> population;//TODO: think whether using a different structure is a good idea
	//generate some initial solutions differing from each other
	population.push_back(generateRandomSolution());
	for (unsigned int i = 0; i < initialPopulation; ++i) {
		TSPSolution random = generateRandomSolution();
		do {
			random = generateRandomSolution();
		} while (random == population[population.size()-1]);//to ensure no repeats
		population.push_back(random);
	}
	
	//main loop
	for (unsigned int loopIteration = 0; loopIteration < 100; ++loopIteration) {
		//generate children
		unsigned int currentPopulation = population.size();
		for (unsigned int i = 0; i < currentPopulation-1; i = i + 2) {
			switch (crossoverOp) {
				case 1:
					population.push_back(cycleCrossover(population[i], population[i + 1]));
					population.push_back(cycleCrossover(population[i + 1], population[i]));
				break;
				default:
					population.push_back(orderCrossover(population[i], population[i + 1]));
					population.push_back(orderCrossover(population[i + 1], population[i]));
				break;
			}
		}
		switch (selectionMethod) {
		case 1: 
			rouletteWheelSelection(population, initialPopulation);
			break;
		default:
			cutoffSelection(population, initialPopulation);
			break;
		}
	}
	return population[0];
}
