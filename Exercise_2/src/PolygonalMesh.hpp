#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary {

struct PolygonalMesh
{
	unsigned int NumCell0Ds;
	unsigned int NumCell1Ds;
	unsigned int NumCell2Ds;
	
	vector<unsigned int> Cell0DsId; 
	vector<unsigned int> Cell1DsId;
	vector<unsigned int> Cell2DsId;
	
	/* the Coordinates and the Extrema are matrices because the Paraview
	export function by Professor Vicini needs Eigen matrices as arguments.*/
	
	Eigen::MatrixXd Cell0DsCoordinates; 
    Eigen::MatrixXi Cell1DsExtrema; 
	
	/*for vertices and edges, I used a vector of vectors. 
	Each vector contains another vector with all the edges and vertices 
	associated to a specific Polygon */
	
	vector<vector<unsigned int>> Cell2DsVertices;
	vector<vector<unsigned int>> Cell2DsEdges;
	
	/*To store the markers, I used a map, which associates an unsigned int
	(the Marker), to the list of points having that marker */
	
	map<unsigned int, list<unsigned int>> Cell0DMarkers;
	map<unsigned int, list<unsigned int>> Cell1DMarkers;
};

}