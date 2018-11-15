#pragma once
#include <vector>

/*
Structures used in both PointTSP and MatrixTSP
*/

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
