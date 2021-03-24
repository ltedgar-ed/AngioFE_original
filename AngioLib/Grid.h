///////////////////////////////////////////////////////////////////////
// Grid.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// The GRID class stores the regular quadrilateral grid that is fit to 
// the simulation domain.  The nodes of the grid store field 
// information which influences microvessel growth.  Interpolation within
// the grid is handled by trilinear shape functions.
///////////////////////////////////////////////////////////////////////



#pragma once

#include "Segment.h"
#include "mat3.h"
#include <vector>

class Data;
class Input;
class Elem;
class Node;

using namespace std;

typedef vector<vector<vector<double> > > Stack;

// Define the number of nodes in the x, y, and z direction
const int xnodes = 76;          // Node distribution according to total domain size of 3822 x 2548 x 200
const int ynodes = 51;
const int znodes = 20;


class Grid  
{
public:
	///// GRID: Member Functions /////
	Grid(Input &input);                       // GRID Constructor - Requires loading condition indicator and collagen density input
	virtual ~Grid();

    virtual void create_grid(Input &input);

	// GRID.IsOutsideBox - Inline function that determines if a newly created vessel Segment is outside of the domain
	bool IsOutsideBox(const Segment& seg);

    int findnnumber(Node node);
    
    // GRID.ijk - Given a current postion, relays the I, J, and K nodes that describe the element that contains that position
	void ijk(double xpt, double ypt, double zpt, int &I, int &J, int &K, Data &data);
		
	// GRID.natcoordinates - Accepts a position in global coordinates and determines the position in natural coorindates for the particular grid element
    void natcoordinates(double &xix, double &xiy, double &xiz, double xpt, double ypt, double zpt, int I, int J, int K);
    void natcoord(double &xix, double &xiy, double &xiz, double xpt, double ypt, double zpt, int elem_num);
    
    int findelem(double xpt, double ypt, double zpt);
    int enhanced_findelem(double xpt, double ypt, double zpt);

    // GRID.shapefunctions - Determine the shape function values for a given position in natural coorindates
    void shapefunctions(double (&shapeF)[8], double xix, double xiy, double xiz);
    
    // GRID.shapefun_d1 - Determine the first order derivative of the shape functions for a particular node
    void shapefun_d1(double (&dshapeF)[3], const double xix, const double xiy, const double xiz, int node);
    
    void nattoglobal(double &xpt, double &ypt, double &zpt, double xix, double xiy, double xiz, int elem_num);

	int elem_find_neighbor(int elem_num,int neighbor_id);
	
	double find_density_scale(double coll_den);
    
    //void find_density(Segment &seg, Grid &grid, Data &data);
    //void smooth_density(double alpha[3][3][3], double density[xnodes][ynodes][znodes], int num_nodes[3]);
	//void Gradient(Grid &grid, Data &data);
	//void clear_density(Grid &grid);
    
public:
	int load_cond;
	
	int xnodes;
	int ynodes;
	int znodes;
	int Nn;                                                     // GRID.Nn - Total number of nodes in the grid
	int Ne;                                                     // GRID.Ne - Total number of elements in the grid
	
	double xrange[2];                                           // GRID.xrange - Array containing the minimum and maximum values in the x direction
	double yrange[2];                                           // GRID.yrange - Array containing the minimum and maximum values in the y direction
	double zrange[2];                                           // GRID.zrange - Array containing the minimum and maximum values in the z direction
	
	int num_nodes[3];                                           // GRID.num_nodes - Array containing the number of nodes in the x, y, and z direction
		
	vector<double> x;                                           // GRID.x - Array containing the nodal positions in the x direction
	vector<double> y;                                           // GRID.y - Array containing the nodal positions in the y direction
	vector<double> z;                                           // GRID.z - Array containing the nodal positions in the z direction

    vector<Elem> ebin;
                   
	Stack theta;                                                // GRID.theta - Matrix containing the angle of collagen fiber orientation in the x-y plane for each of the nodes 
	Stack eta;                                                  // GRID.eta - Matrix containing the angle of collagen fiber orientation off of the x-y plane for each of the nodes
	
	Stack coll_fib_x;
	Stack coll_fib_y;
	Stack coll_fib_z;

	vector<Node> nodes;
		
	double coll_den;
	double den_scale;                                           // GRID.den_scale - Scaling factor based off of the collagen density within the domain
	double a, b, c;                                             // GRID.a, GRID.b, GRID.c - Parameters that describe the function that determines den_scale
	                                                            //                          given the command line input of collagen density.
	char x_bctype; 
	char y_bctype;
	char z_bctype;

	char frontbc;
	char rightbc;
	char backbc;
	char leftbc;
	char bottombc;
	char topbc;
};




   
    
    
    
///////////////////////////////////////////////////////////////////////
// Inline Member Functions
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////
// IsOutsideBox
///////////////////////////////////////////////////////////////////////////////////

// GRID.IsOutsideBox - Determines if a vessel segment has grown outside of the domain
//      Input:       const Segment& seg - Newly created vessel Segment
//
//      Output:      Boolean operator 'true' if segment is outside the domain

inline bool Grid::IsOutsideBox(const Segment& seg) //determine if segment is outside box
{
	if (seg.x[0] < xrange[0] || seg.x[0] > xrange[1] ||
		seg.x[1] < xrange[0] || seg.x[1] > xrange[1] ||
		seg.y[0] < yrange[0] || seg.y[0] > yrange[1] ||
		seg.y[1] < yrange[0] || seg.y[1] > yrange[1] ||
		seg.z[0] < zrange[0] || seg.z[0] > zrange[1] ||
		seg.z[1] < zrange[0] || seg.z[1] > zrange[1])
		return true;
	else
		return false;
}

