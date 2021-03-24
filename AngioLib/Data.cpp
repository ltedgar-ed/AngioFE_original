///////////////////////////////////////////////////////////////////////
// Data.cpp
///////////////////////////////////////////////////////////////////////



// Include:
#include "stdafx.h"
#include <math.h>

#include "Grid.h"
#include "Data.h"
#include "Input.h"
#include "angio3d.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Data::Data(Input &input, Grid &grid)                                // Constructor for DATA object
{                                                               // Input: Branching chance defined in command line, GRID data object
	
	///// BEGIN USER INPUT /////
	
	NFRAGS = input.infrag;                                      // NFRAGS - Define the number of intial fragments
	maxt = static_cast<double>(input.imax_time);                // Total number of days (Default: 7 days)
	dt = static_cast<double>(input.idt);                        // Initial time step (Default: 0.25 days)
    length_adjust = static_cast<double>(input.ilngth_adj);      // length_adjust - Define factor used to scale the length of new segments (no scale = 1.0)  
    anast_dist = static_cast<double>(input.ianst_dst);          // anast_dist - Distance threshold for anastomosis (in um)    
    branch_chance = input.ibrnch_ch;                            // branch_chance - Define from command line argument 'bc'
    
    ///// END USER INPUT /////
    
    
    t = 0.0;                                                    // Initial time (Default: 0 days)
    dt_a = 0.0637;
    dt_b = 9.0957;
    dt_c = 2.6073;
    n = 1.0;

    total_length = 0.;
    
    a1 = -1.2653;                                               // a1, a2, a3 - Parameters for branching curve
	a2 = 1.535;
	a3 = 0.1108;

    a = 1900.0;                                              // a, b, x0, y0 - Parameters for growth curve
    b = 1.4549;
    x0 = 4.9474;
    y0 = -19.1278;

    d = y0 + a/(1+pow(E,x0/b));                                 // d - Initial value of growth curve (t = 0)

	vessel_width = 7;                                           // vessel_width - Diameter of microvessels (Default: 7 um)
	num_anastom = 0;                                            // num_anastom - Initialize anastomoses counter
	num_branches = 0;                                           // num_branches - Initialize branching counter
	num_zdead = 0;
	
	branch = false;                                             // branch - Set branching indicator to 'false'

	dx = (grid.xrange[1] - grid.xrange[0])/(grid.num_nodes[0]-1);   // dx - Spatial step size in x direction, range of x divided by the number of nodes
	dy = (grid.yrange[1] - grid.yrange[0])/(grid.num_nodes[1]-1);   // dy - Spatial step size in y direction, range of y divided by the number of nodes
	dz = (grid.zrange[1] - grid.zrange[0])/(grid.num_nodes[2]-1);   // dz - Spatial step size in z direction, range of z divided by the number of nodes

	nsegs = 0;                                                  // nsegs - Initialize segment counter
    num_vessel = NFRAGS - 1;                                    // num_vessel - Initialize vessel counter
        
	vess_length = d;
	old_length = d;

	num_lines = -1;
}


Data::~Data()                                                   // Destructor for DATA object
{

}
