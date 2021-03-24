///////////////////////////////////////////////////////////////////////
// Elem.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
//  The ELEM class contains the basic hexahedral element used to 
//  construct the grid.  The ELEM class also contains the NODE class
//  and the FACE class. 
///////////////////////////////////////////////////////////////////////



#pragma once

#include <list>
#include "vect3.h"
#include <vector>
class Segment;
using namespace std;

class Node
{
  public:  
    double x;
    double y;
    double z;
    
	double x0;
	double y0;
	double z0;

    double theta;
    double eta;
    
	double theta0;
	double eta0;

	double ecm_den;
    double ecm_den0;
	
	int id;
        
	bool updated;

	vect3 ecm_den_grad;
	vect3 u;
	vect3 collfib;
	vect3 collfib0;

	vector<double> ecm_den_store;
	vector<vect3> ecm_fibril_store;

 public:
    Node() : x(0.), y(0.), z(0.), x0(0.), y0(0.), z0(0.), theta(0.), eta(0.), theta0(0.), eta0(0.), ecm_den(0.), ecm_den0(0.), id(0), updated(false) {}

	bool operator == (const Node& n2)
	{
		if ((x == n2.x) && (y == n2.y) && (z == n2.z))
			return true;
		else 
			return false;
	}
};



class Face
{
  public:
    bool BC;
    char bc_type;
    
  public:
    Face() : BC(false), bc_type('n') {}
};



class Elem
{
  ///// ELEM: Member Functions /////
  public:
    Elem() : elem_num(-1), volume(0.), volume0(0.), alpha(0.) {}
	  
	  // Find the dimensions of the bounding box for the element
    double bb_xmin();
    double bb_xmax();
    double bb_ymin();
    double bb_ymax();
	double bb_zmin();
    double bb_zmax();
    
    // Find the dimensions of the inner box for the element
    double ib_xmin();
    double ib_xmax();
    double ib_ymin();
    double ib_ymax();
	double ib_zmin();
    double ib_zmax();
    
  ///// ELEM: Member Fields /////
  public:
	int elem_num;       // Element Identifier
    
	double volume;
	double volume0;

	double alpha;
	vect3 fiber_orient;

	//list<list<Segment>::iterator > resident_segs;

    Node * n1;            // Bottom, lower, left node
    Node * n2;            // Bottom, lower, right node
    Node * n3;            // Bottom, upper, left node
    Node * n4;            // Bottom, upper, right node
    Node * n5;            // Top, lower, left node 
    Node * n6;            // Top, lower, right node    
    Node * n7;            // Top, upper, left node
    Node * n8;            // Top, upper, right node
    
    Face f1;            // Front face -y
    Face f2;            // Right face +x
    Face f3;            // Back face +y
    Face f4;            // Left face -x
    Face f5;            // Top face +z
    Face f6;            // Bottom face -z

};



