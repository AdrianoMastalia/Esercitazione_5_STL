#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Eigen/Eigen"
using namespace Eigen;
using namespace std;
namespace PolygonalLibrary
{

	
bool ImportMesh(PolygonalMesh& mesh)
{
	/*I try to read the file containing the 0d cells, in case of 
	success I write the markers with the associated points, so I can
	compare the output with the image provided in the exercise */
	
	if(!ImportCell0Ds(mesh))
	{
		cout<<"Error in the reading of 0Ds properties";
		return false;
	}
	else
	{
		cout<<"Markers associated to their points: "<<endl;
		for(auto &Marker : mesh.Cell0DMarkers)
		{
			cout<<"Marker value "<<Marker.first<<", Point IDs associated to the Marker: "<<endl;
			for(auto element : Marker.second)
				cout<<element<<" ";
			cout<<endl;
		}
	}
	if(!ImportCell1Ds(mesh))
	{
		cout<<"Error in the reading of 1Ds properties";
		return false;
	}
	else
	{
		cout<<"Markers associated to each segment: "<<endl;
		for(auto &Marker : mesh.Cell1DMarkers)
		{
			cout<<"Marker value "<<Marker.first<<", segment IDs associated to the Marker: "<<endl;
			for(auto element : Marker.second)
				cout<<element<<" ";
			cout<<endl;
		}
	}
	
	/*As I didn't track the points associated to the 0 marker, because in applications 
	we are interested in the points on the boundary only. As all polygon is associated to the 
	marker 0, I didnt't store this information, so I don'have to print it here neither */
	
	if(!ImportCell2Ds(mesh))
	{
		cout<<"Error in the reading of 2Ds properties";
		return false;
	}
	
	/*with the following functions, I control that there are no zero length edges
	as well as no zero area polygons */
	
	edges_length(mesh);
	check_area(mesh);
	return true;
} 


bool ImportCell0Ds(PolygonalMesh& mesh)
{
	/* /!\ IMPORTANT: THE FOLLOWING CODE WORKS IF WHE CONSIDER THE STRUCT MESH WITH
	DIFFERENT PROPERTIES: IN PARTICULAR HERE I CONSIDERED THE COORDINATES AND EXTREMA
	NOT AS MATRICES BUT AS VECTORS OF 2D EIGEN VECTORS. AS SAID IN THE POLYGONALMESH.HPP FILE,
	THE CURRENT CODE USES EIGEN MATRICES IN ORDER TO EXPORT THE IMAGES OF THE MESH, WITHOUT
	THE IMAGE TO BE EXPORTED, THIS WOULD HAVE BEEN ANOTHER WAY TO STORE THE MESH, USING ONLY ONE 
	WHILE LOOP /!\ */
	
    /*ifstream file("Cell0Ds.csv");
    if (!file.is_open())
	{
		cout<<"Errore nell'apertura del file Array";
        return false;
	}

    string tmp;
	getline(file, tmp);
	unsigned int id;
	unsigned int marker;
	Vector2d coordinates;
    while (getline(file, tmp))
    {
		replace(tmp.begin(),tmp.end(),';',' ');
		istringstream datas0D(tmp);
		datas0D >> id >> marker >> mesh.Cell0DsCoordinates(0, id) >> mesh.Cell0DsCoordinates(1, id);
		mesh.Cell0DsId.push_back(id);
		if(marker != 0)
		{
			auto it = mesh.Cell0DMarkers.find(marker);
			if(it != mesh.Cell0DMarkers.end())
				mesh.Cell0DMarkers[marker].push_back(id);
			else
				mesh.Cell0DMarkers.insert({marker, {id}});
		}
    }
	mesh.NumCell0Ds = mesh.Cell0DsId.size();
	file.close();
	return true;*/
	
	ifstream file("Cell0Ds.csv");

    if(file.fail())
        return false;
	
	/*I put all the lines of the file in this list of strings */
	
    list<string> listLines;

    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    /* I remove the header of the file from the list */
	
    listLines.pop_front();
	/*the remaining number of lines equals the number of 0d Cells */
	
    mesh.NumCell0Ds = listLines.size();
    if (mesh.NumCell0Ds == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }
	/*Here, I allocate memory for the vector of Ids and I crete the coordinates as a 
	3 x number_of_cells matrix. The third line of the matrix is for an hypothetical z coordinate, which 
	will be equal to zero in our case, as we are dealing with a two-dimensional mesh */
	
    mesh.Cell0DsId.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (string& line : listLines)
    {
		/*Here, I replace the space with ; and convert the string to a stringstream, so I can take 
		the data more easily out of it */
		
		replace(line.begin(),line.end(),';',' ');
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2d coord;
		
		/*Now, in the first row of the Cell0DsCoordinates matrix we'll have the x-values of the points 
		of the mesh, and in the second row the y-coordinates. The id of the point is used to denote 
		the column in wich the poin (x,y) can be found */
		
        converter >>  id >> marker >> mesh.Cell0DsCoordinates(0, id) >> mesh.Cell0DsCoordinates(1, id);
        mesh.Cell0DsId.push_back(id);

	/*If the current marker is different from zero (so the point belongs to the boundary of the mesh), 
	I have to memorize it in the map */
	
        if(marker != 0)
		{
			/*I search the marker in the ones already memorized as keys in the map,
			and if it hasn't been memorized yet, I add the couple (marker, {id}), where 
			the second term is a list composed by only one element, the id of the current point */
			
			auto it = mesh.Cell0DMarkers.find(marker);
			if(it != mesh.Cell0DMarkers.end())
				mesh.Cell0DMarkers[marker].push_back(id);
			else
				mesh.Cell0DMarkers.insert({marker, {id}});
		}

    }

    return true;
}

/*In the followin function, I implement the same stuff as in the function for 0d Cells, the
first commented part is the other way to store data of the mesh using vectors and not eigen matrices */

bool ImportCell1Ds(PolygonalMesh& mesh)
{
	/* ifstream file("Cell1Ds.csv");
    if (!file.is_open())
	{
		cout<<"Errore nell'apertura del file Array";
        return false;
	}

    string tmp;
	getline(file, tmp);
	unsigned int id;
	unsigned int marker;
	Vector2i Origin_End;
    while (getline(file, tmp))
    {
		replace(tmp.begin(),tmp.end(),';',' ');
		istringstream datas1D(tmp);
		datas1D >> id >> marker >> mesh.Cell1DsExtrema(0, id) >>  mesh.Cell1DsExtrema(1, id);
		mesh.Cell1DsId.push_back(id);
		
		if(marker != 0)
		{
			auto it = mesh.Cell1DMarkers.find(marker);
			if(it != mesh.Cell1DMarkers.end())
				mesh.Cell1DMarkers[marker].push_back(id);
			else
				mesh.Cell1DMarkers.insert({marker, {id}});
		}
    }
	mesh.NumCell1Ds = mesh.Cell1DsId.size();
	file.close();
	return true; */
	
	
	ifstream file("Cell1Ds.csv");

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();

    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();
    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsId.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

    for (string& line : listLines)
    {
		replace(line.begin(),line.end(),';',' ');
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2i vertices;
		
		/*Here, the matrix Cell1DsExtrema will have two rows, the first with the origins of the edges, and
		the second with the ends of the edges. As in the coordinates case, the columns indicate the id of
		the edge (Origin, End)*/
		
        converter >>  id >> marker >>  mesh.Cell1DsExtrema(0, id) >>  mesh.Cell1DsExtrema(1, id);
        mesh.Cell1DsId.push_back(id);
		if(marker != 0)
		{
			auto it = mesh.Cell1DMarkers.find(marker);
			if(it != mesh.Cell1DMarkers.end())
				mesh.Cell1DMarkers[marker].push_back(id);
			else
				mesh.Cell1DMarkers.insert({marker, {id}});
		}
        
    }
    return true;
}


bool ImportCell2Ds(PolygonalMesh& mesh)
{
	ifstream file("Cell2Ds.csv");
    if (!file.is_open())
	{
		cout<<"Errore nell'apertura del file Array";
        return false;
	}

    string tmp;
	getline(file, tmp);
	unsigned int id;
	unsigned int marker;
	unsigned int num_vertices;
	unsigned int num_edges;
	
	/*In this part, there is a difference with the other two functions, as I don't have to specify
	the dimensions of any matrix, I can use a single while loop to store the lines of the file*/
	
	 while (getline(file, tmp))
    {
		replace(tmp.begin(),tmp.end(),';',' ');
		istringstream datas2D(tmp);
		
		/*As I don't know the number of vertices of the polygon, 
		I have to memorize the line up to this point*/
		
		datas2D >> id >> marker >> num_vertices;
		mesh.Cell2DsId.push_back(id);
		
		/*I create the vector of vertices with the correct size, and then I memorize all
		the vertices associated to a cell in a vector. This vector is then appended in the bigger
		vector containing all the vertices of the polygons */
		
		vector<unsigned int> Vertices(num_vertices);
		for(int i = 0; i<num_vertices; i++)
			datas2D >> Vertices[i];
		mesh.Cell2DsVertices.push_back(Vertices);
		
		/*then I repeat the procedure with the edges of each polygon */
		
		datas2D>>num_edges;
		vector<unsigned int> Edges(num_edges);
		for(int i = 0; i<num_edges; i++)
			datas2D >> Edges[i];
		mesh.Cell2DsEdges.push_back(Edges);
    }
	mesh.NumCell2Ds = mesh.Cell2DsId.size();
	file.close();
    return true;
}

/*the following functions control that there are no zero length edges in the mesh */

bool edges_length(PolygonalMesh& mesh)
{
	
	/*With the outer loop, we are considering every different polygon of the mesh, the inner loop
	distinguishes the different edges associated to the i-th polygon */
	
	for(unsigned int i = 0; i<mesh.NumCell2Ds; i++)
	{
		for(unsigned int j = 0; j<mesh.Cell2DsEdges[i].size();j++)
		{	
			/* the vector edges contains all the edges of the i-th polygon.
			at edges[j] I can find the id of the j-th edge of the i-th polygon*/
			
			vector<unsigned int>& edges = mesh.Cell2DsEdges[i];
			
			/*In the Cell1DsExtrema matrix, at column edges[j] there is the couple (O,E) of Origin and 
			End IDs of the j-th edge in the i-th polygon (recall that each edge is formed by two points). 
			I store in the variables Origin_index and End_Index the IDs of Origin and End of the edge */
			
			int& Origin_index = mesh.Cell1DsExtrema(0,edges[j]);
			int& End_index = mesh.Cell1DsExtrema(1,edges[j]);
			
			/*Having the ID of the origin and end of the edge, I can take their coordinates, remembering that in 
			Cell0DsCoordinates, in the column denoted by Origin_index and End_index, there are the (x,y) coordinates
			of the Origin and End points of the current edge. So with Cell0DsCoordinates(0,Origin_index) I store the x-coordinate 
			of the origin, with Cell0DsCoordinates(1,Origin_index) I store the y-coordinate of the origin and so on*/
			
			double& X_Origin = mesh.Cell0DsCoordinates(0,Origin_index);
			double& Y_Origin = mesh.Cell0DsCoordinates(1,Origin_index);
			double& X_End = mesh.Cell0DsCoordinates(0,End_index);
			double& Y_End = mesh.Cell0DsCoordinates(1,End_index);	
			
			/*the distance between points is the euclidean Distance, given by Pithagora's Theorem */
			
			double distance = sqrt(pow(X_Origin-X_End,2)+pow(Y_Origin-Y_End,2));
			
			/*Here, I set a tolerance of 1e-16 on the length of the edge, If there are smaller edges in the mesh, 
			they will be counted as zero length edges. This is due to the fact that numerical cancellation can be dangerous
			if the points are too close one another */
			
			if(distance < 1e-16)
			{
				cout<<"There is an error at polygon with ID "<< i <<", edge with ID "<<edges[j]<<" has 0 length"<<endl;
				return false;
			}
		}
	}
	cout<<"There are no length zero edges "<<endl;
	return true;
}

/* /!\ IMPORTANT: IN THE FOLLOWING FUNCTION, I CALCULATE THE AREA OF EVERY POLYGON AND MAKE SURE THAT IT'S DIFFERENT
FROM ZERO, BUT I USED THE FORMULA GIVEN BY PROFESSOR VICINI, WHICH DOESN'T CONSIDER POLYGONS WITH OVERLAPPING EDGES. 
SO THIS PART OF THE CODE WORKS ONLY WITH THE ASSUMPTION THAT THERE ARE NO OVERLAPPING EDGES OF THE POLYGON /!\ */

bool check_area(PolygonalMesh& mesh)
{
	/*The outer loop concerns the i-th polygon of the mesh, while the inner loop takes the 
	j-th vertex of the i-th polygon */
	
	for(unsigned int i = 0; i<mesh.NumCell2Ds; i++)
	{
		double area = 0.0;
		unsigned int n = mesh.Cell2DsVertices[i].size();
		for(unsigned int j = 0; j < n; j++)
		{
			/* Here, I take the IDs of the points which form the i-th polygon. Note that, in order to 
			apply the formula seen in the course theory, the points have to be ordered with the last one going
			back to the first point, one way to do this is by considering the modulo n operation on the point P2
			(recall that n is the number of vertices of the i-th polygon, so everything is well-defined) */
			
			unsigned int& P1_id = mesh.Cell2DsVertices[i][j];
			unsigned int& P2_id = mesh.Cell2DsVertices[i][(j+1)%n];
			
			/*As done in the function above, inside Cell0DsCoordinates, at column given by P1_id or P2_id, 
			there are the (x,y) coordinates of the vertices P1 and P2, so I can access these data considering the first (0)
			or second (1) row of the matrix, at column P1_id or P2_id. */
			
			double& X_P1 = mesh.Cell0DsCoordinates(0,P1_id);
			double& Y_P1 = mesh.Cell0DsCoordinates(1,P1_id);
			double& X_P2 = mesh.Cell0DsCoordinates(0,P2_id);
			double& Y_P2 = mesh.Cell0DsCoordinates(1,P2_id);
			
			area += X_P1*Y_P2-X_P2*Y_P1;
		}
		
		area = 0.5*abs(area);
		
		/*I also put a tolerance on the area of the polygon, the reason is the same as the one described in the previous function */
		
		if(area < 1e-12)
		{
			cout<<"There is an error at polygon with ID "<<i<<": it has zero area"<<endl;
			return false;
		}
	}
	cout<<"There are no zero area polygons"<<endl;
	return true;
}

}
