///////////////////////////////////////////////////////////////////////
// Segment.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// Microvessels are represent by a collection of line segments. 
// Growth is represented by the addition of new segments onto the 
// active tips of exisiting segments. Within the simulation, these 
// line segments are found in the SEGMENT class.
///////////////////////////////////////////////////////////////////////



#pragma once

#include "vect3.h"

class Segment  
{
  ///// SEGMENT: Member Functions /////
  public:
    Segment();                                                  // SEGMENT Constructor
	virtual ~Segment();
    
    void findlength();
    void findunit();
	//void findphi();

  ///// SEGMENT: Member Fields /////
  public:
	double x[2];                                                // SEGMENT.x - Array containing the x-coordinates of the segment tips
	double y[2];                                                // SEGMENT.y - Array containing the y-coordinates of the segment tips
	double z[2];                                                // SEGMENT.z - Array containing the z-coordinates of the segment tips
	
    double x0[2];
	double y0[2];
	double z0[2];
	
    vect3 uvect;
	
	int tip[2];                                                 // SEGMENT.tip - Array indicating that the segment's nodes are active
    double length;                                              // SEGMENT.length - Length of the segment (in um), Euclidean distance between the nodes
	
	//double phi1;                                                // SEGMENT.phi1 - Orientation of the segment within the x-y plane, phi1 is 
	                                                            //                the angle between the segment and the x-axis (in radians)
	//double phi2;                                                // SEGMENT.phi2 - Orientation of the segment off of the x-y plane, phi2 is
	                                                            //                the angle between the segment and the x-y plane (in radians)
	
	int label;                                                  // SEGMENT.label - Label that indicates which initial fragment the segment orginated from
	int vessel;                                                 // SEGMENT.vessel - Label that indicates which vessel the segment belongs to
	int seg_num;
	
	int BCdead;                                                 // SEGMENT.BCdead - Boolean flag that indicates that the segment has encountered a boundary condition
	double TofBirth;                                            // SEGMENT.TofBirth - Time point at which the segment was created
	double Recent_branch;                                       // SEGMENT.Recent_branch - Indicates that the segment was recently involved in the formation of a branch
	bool init_branch;                                           // SEGMENT.init_branch - Boolean flag that indicates whether or not the initial fragment should form a branch
	                                                            //                       at t = 0
	int sprout;                                                 // SEGMENT.sprout - Marker that indicates which end of the parent vessel the segment sprouted from.  1 for the +1 end, 
	                                                            //                  -1 for the -1 end, 9 for an initially seeded fragment,
	                                                            //                  0 for a frament not touched by the growth routine for some reason.
    int anast;
    
	int tip_elem[2];                                               // SEGMENT.seg_elem - Indicates which element the segment is currently occupying

	bool elem_tagged;

	int bdyf_id[2];

	bool mark_of_death;
	int death_label;

	int tip_BC[2];

	int seg_conn[2][2];

	//int line_num;
};


