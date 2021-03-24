///////////////////////////////////////////////////////////////////////
// Culture.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// The CULTURE class contains all the functions that describe how the 
// SEGMENT class is to grow and orient itself. These functions are the 
// rule-set that arrange the line segments into vascular networks, 
// mimicking angiogenesis. 
///////////////////////////////////////////////////////////////////////



#pragma once

#include <list>
#include "BC.h"

class Grid;
class Data;
class Segment;

using namespace std;


class Culture  
{
  ///// CULTURE: Member Functions ///// 
  public:
	Culture();                                                  // CULTURE Constructor
	virtual ~Culture();

    // CULTURE.createInitFrag - Seed the initial fragments within the domain
	Segment createInitFrag(Data &data, Grid &grid, int i,list<Segment> &frag);
	
	// CULTURE.createNewSeg - Create a new segment at the tip of an existing segment
	Segment createNewSeg(list<Segment>::iterator it, Grid &grid, Data &data,int k, list<Segment> &frag);
	
	// CULTURE.findLength - Determine the length of a newly created segment based on the growth function g(t)
	//double findLength(double xpt, double ypt, double zpt, Grid &grid, Data &data);
	
	// CULTURE.findAngle - Determine the orientation angle of a newly created segment based on the information stored in GRID
	vect3 findAngle(list<Segment>::iterator, double xpt, double ypt, double zpt, Grid &grid, Data &data);
	
	// CULTURE.findCollAngle - Obtain the component of new vessel orientation determined by local collagen fiber orientation 
	vect3 findCollAngle(double xpt, double ypt, double zpt, Grid &grid, Data &data);
	
	double findDenScale(double xpt, double ypt, double zpt, Grid &grid);
	
	//double findDenAngle(list<Segment>::iterator it, double xpt, double ypt, double zpt, Grid &grid, Data &data, int ang_type);

	// CULTURE.connectSegment - Create a new segment connecting two existing segments that are fusing through anastomosis
	Segment connectSegment(list<Segment>::iterator it, list<Segment>::iterator it2, int k, int kk, Grid &grid, Data &data, list<Segment> &frag);
	
	// CULTURE.CheckForIntersection - Check a newly created segment to see if it physically intersections with any existing segments
	void CheckForIntersection(Segment &seg,list<Segment> &frag, Data &data, list<Segment>::iterator it);
	
	// CULTURE.findIntersect - Find the coordinates at which two segments intersect
	double findIntersect(double a[3], double b[3], double c[3], double d[3], double intersectpt[3]);
	
	// CULTURE.intersectPlane - Determine if a segment encounters one of the boundary planes, find the coordinates of the intersection point
	bool intersectPlane(Grid& grid, Segment &Seg, int n, double intersectpt[3]);
	
	// CULTURE.PeriodicBC - If a segment encounters one of the boundary planes, enforce the periodic boundary conditions
	Segment PeriodicBC(Segment &seg,Grid &grid,list<Segment> &frag,Data &data);

  ///// CULTURE: Member Fields /////	
  public:
    BC bc;

		double W[4];                                                                        // W[1] = Weight for vessel density
};


