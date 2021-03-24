///////////////////////////////////////////////////////////////////////
// Data.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// The DATA class stores all the important parameters about the 
// simulation.
///////////////////////////////////////////////////////////////////////



#pragma once

class Grid;
class Input;

class Data  
{
  ///// DATA: Member Functions ///// 
  public:
    Data(Input &input, Grid &grid);                                 // DATA Constructor - Requires boundary condition indicator and GRID
	virtual ~Data();

  ///// DATA: Member Fields /////
  public:
    int NFRAGS;                                                 // DATA.NFRAGS - Number of initial microvessel fragments
	double branch_chance;                                       // DATA.branch_chance - Probability of forming a new branch
	double maxt;                                                // DATA.maxt - Time at which the simulation ends (in days)
	double t;                                                   // DATA.t - Current time (in days)
	double dt;                                                  // DATA.dt - Variable time step (in days)
	
	double dt_a, dt_b, dt_c;
    double n;    
    
    double a1,a2,a3;                                            // DATA.a1, DATA.a2, DATA.a3 - Parameters for the branching curve

	                                                            // Microvessel growth rate is modeled using the Boltzmann sigmoid curve 
	                                                            //          g(t) = y0 + a/(1+e^-(t-x0)/b)
	double y0;                                                  // DATA.y0 - Bottom of sigmoid curve
	double a;                                                   // DATA.a - Distance between top and bottom of the curve (a + y0 = top of curve)
	double x0;                                                  // DATA.x0 - Time point at which t is halfway between top & bottom
	double b;                                                   // DATA.b - Steepness of the curve
	double d;                                                   // DATA.d - Initial value of the curve (t = 0)

	int num_branches;                                           // DATA.num_branches - Counter indicating the number of branches formed during the simulation
	int num_anastom;                                            // DATA.num_anastom - Counter indicating the number of anastomoses formed during the simulation
	int num_zdead;
	int nsegs;                                                  // DATA.nsegs - Counter that stores in current number of Segments within the simulation domain
	                                        
	double total_length;                                        // DATA.total_length - Total vascular length within the domain (sum of the length of all Segments) (in um)
    double vessel_width;                                        // DATA.vessel_width - Width of Segments (in um) 
    
    double dx;                                                  // DATA.dx - Spatial discretization step size in the x-direction (in um)
	double dy;                                                  // DATA.dy - Spatial discretization step size in the y-direction (in um)
	double dz;                                                  // DATA.dz - Spatial discretization step size in the z-direction (in um)
	
	bool branch;                                                // DATA.branch - Boolean flag that indicates to the model that the Segment being created is the result of a new branch
    int num_vessel;                                             // DATA.num_vessel - Counter that indicates the next vessel ID number of new Segments
    double length_adjust;                                       // DATA.length_adjust - Length adjuster, scales the length of new segments, allows fine tuning of total vascular length produced
                                                                //                      by the simulation.
    double anast_dist;                                          // DATA.anast_dist - If the shortest distance vector between a segment and any other segment is less than this value,
                                                                //                   then the microvessels fuse through anastomosis.
    
    double vess_length;                                                       
    double old_length;
    
    time_t ranseed;

	int num_lines;                                                                
};


