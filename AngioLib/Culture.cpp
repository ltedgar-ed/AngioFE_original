///////////////////////////////////////////////////////////////////////
// Culture.cpp
///////////////////////////////////////////////////////////////////////



// Include:
#include "stdafx.h"
#include <iostream>

#include "Culture.h"
#include "Data.h"
#include "Grid.h"
#include "Elem.h"

using namespace std;


///////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////////

Culture::Culture()                                              // Constructor for CULTURE object
{
	W[0] = 0;
	W[1] = 0;
	W[2] = 0;
	W[3] = 0;
}

Culture::~Culture()                                             // Destructor for CULTURE object
{

}



///////////////////////////////////////////////////////////////////////
// Member Functions
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// createInitFrag
///////////////////////////////////////////////////////////////////////

// CULTURE.createInitFrag - Seed an initial fragment
//      Input:  - DATA object
//              - GRID object
//              - Index that indicates which initial fragment this segment is (i)
//              - Segment container (frag)
//
//      Output: - Initial fragment (seg)

Segment Culture::createInitFrag(Data &data, Grid &grid, int i, list<Segment> &frag)
{
    Segment seg;                                                // Declare SEGMENT seg
	seg.length = data.d;                                        // Set seg length to value of growth function at t = 0
    
    int elem_num = 0;
    double xix, xiy, xiz = {0.};
    double xpt, ypt, zpt = {0.};
    
	//// Seed normally
	elem_num = int((float(rand())/RAND_MAX)*grid.Ne);
	xix = 2*((float(rand())/RAND_MAX) - 0.5);
	xiy = 2*((float(rand())/RAND_MAX) - 0.5);
	xiz = 2*((float(rand())/RAND_MAX) - 0.5);
/*	elem_num = 4609;//247//665;//int((float(rand())/RAND_MAX)*grid.Ne);
	xix = 0;//2*((float(rand())/RAND_MAX) - 0.5);
	xiy = 0;//2*((float(rand())/RAND_MAX) - 0.5);
	xiz = 0;//2*((float(rand())/RAND_MAX) - 0.5);
*/    grid.nattoglobal(xpt, ypt, zpt, xix, xiy, xiz, elem_num);
	
	//// Seed outside of the plug only
	//bool out_plug = false;
	//double plug_size = 1250.;
	//double r = 0.;
	//double x_cent = (grid.xrange[1] - grid.xrange[0])/2;
	//double y_cent = (grid.yrange[1] - grid.yrange[0])/2;
	//double z_cent = (grid.zrange[1] - grid.zrange[0])/2;
	//
	//while (out_plug == false){
	//	elem_num = int((float(rand())/RAND_MAX)*grid.Ne);
	//	xix = 2*((float(rand())/RAND_MAX) - 0.5);
	//	xiy = 2*((float(rand())/RAND_MAX) - 0.5);
	//	xiz = 2*((float(rand())/RAND_MAX) - 0.5);
 //   
	//	grid.nattoglobal(xpt, ypt, zpt, xix, xiy, xiz, elem_num);
	//	r = sqrt(pow(xpt - x_cent,2) + pow(ypt - y_cent,2) + pow(zpt - z_cent,2));

	//	if (r >= plug_size)
	//		out_plug = true;}
    
	//// Seed inside microchannels only
	//bool out_channel = false;
	//
	//while (out_channel == false){
	//	elem_num = int((float(rand())/RAND_MAX)*grid.Ne);
	//	xix = 2*((float(rand())/RAND_MAX) - 0.5);
	//	xiy = 2*((float(rand())/RAND_MAX) - 0.5);
	//	xiz = 2*((float(rand())/RAND_MAX) - 0.5);
	//	grid.nattoglobal(xpt, ypt, zpt, xix, xiy, xiz, elem_num);

	//	if (((ypt >= 774) && (ypt <= 1774)) || ((ypt >= 3322) && (ypt <= 4322)))
	//		out_channel = true;}
	

	// Determine vessel orientation based off of collagen fiber orientation
    seg.x[0] = xpt;                                      // Define xpt as x-coordinate of segment tip
    seg.y[0] = ypt;                                      // Define ypt as y-coordinate of segment tip
    seg.z[0] = zpt;                                      // Define zpt as z-coordinate of segment tip
    seg.tip_elem[0] = elem_num;
    
    //seg.phi1 = findCollAngle(xpt, ypt, zpt, grid, data, 1);     // Determine phi1 for the segment using the collagen fiber orientation stored in GRID
    //seg.phi2 = findCollAngle(xpt, ypt, zpt, grid, data, 2);     // Determine phi2 for the segment using the collagen fiber orientation stored in GRID
	//seg.phi2 = 0.75*seg.phi2;                                    // Scale phi2 to make initial fragments more planar (x-y plane)

	seg.uvect = findCollAngle(xpt, ypt, zpt, grid, data);
	//seg.uvect.z = 0.2*seg.uvect.z;

	//// End of new segment is origin plus length component in each direction	
	//seg.x[1] = seg.x[0]+seg.length*cos(seg.phi2)*cos(seg.phi1); // Determine the x-coordinate of the end point using the length vector and orientation angles                 
	//seg.y[1] = seg.y[0]+seg.length*cos(seg.phi2)*sin(seg.phi1); // Determine the y-coorindate of the end point using the length vector and orientation angles
	//seg.z[1] = seg.z[0]+seg.length*sin(seg.phi2);               // Determine the z-coordinate of the end point using the length vector and orientation angles

	// End of new segment is origin plus length component in each direction	
	seg.x[1] = seg.x[0] + seg.length*seg.uvect.x;				  // Determine the x-coordinate of the end point using the length vector and orientation angles                 
	seg.y[1] = seg.y[0] + seg.length*seg.uvect.y;				  // Determine the y-coorindate of the end point using the length vector and orientation angles
	seg.z[1] = seg.z[0] + seg.length*seg.uvect.z;                 // Determine the z-coordinate of the end point using the length vector and orientation angles
	seg.findlength();

	seg.tip[0] = -1;                                            // Set the tip at the start point of the segment as -1 
	seg.tip[1] = 1;                                             // Set the tip at the end point of the segment as +1
	
	data.total_length = data.total_length + fabs(seg.length);    // Update total length
	
	seg.label = i;                                              // Give the segment the appropriate label
	seg.vessel = seg.label;                                     // Set the segment vessel as the segment label
	seg.TofBirth = data.t;                                      // Store the segment's time of birth
	
	seg.sprout = 9;                                             // Set sprout as 9, indicating that this segment is an initial fragment
	
	elem_num = grid.findelem(seg.x[1], seg.y[1], seg.z[1]);
		
	if (elem_num == -1)
		bc.checkBC(seg, 1, grid, data, frag);
	else
		seg.tip_elem[1] = elem_num;


	seg.x0[0] = seg.x[0];
	seg.x0[1] = seg.x[1];
	seg.y0[0] = seg.y[0];
	seg.y0[1] = seg.y[1];
	seg.z0[0] = seg.z[0];
	seg.z0[1] = seg.z[1];	
	
	//data.num_lines++;
	//seg.line_num = data.num_lines;
	
	return seg;                                                 // Return the new segment 
}



///////////////////////////////////////////////////////////////////////
// createNewSeg
///////////////////////////////////////////////////////////////////////

// CULTURE.createNewSeg - Create a new segment at the tip of an existing segment
//      Input:  - Iterator for Segment container 'frag', points to the parent segment (it)
//              - GRID object
//              - DATA object
//              - Index indicating which tip of the parent segment is forming the new segment (k)
//              - Segment container (frag)
//
//      Output: - Newly created segment (seg)

Segment Culture::createNewSeg(list<Segment>::iterator it,Grid &grid, Data &data, int k, list<Segment> &frag)
{
	Segment seg;                                                // Declare SEGMENT seg
	seg.Recent_branch = it->Recent_branch/2;                    // Halve the Recent_branch indicator
    
	++data.nsegs;                                               // Iterate the total segment counter +1 
	seg.seg_num = data.nsegs;				
	
	//seg.line_num = it->line_num;

    int elem_num = 0;
    
	if (it->tip[k] == 1)                                          // If the parent vessel active tip is set as +1...
	{
		//seg.length = findLength(it->x[1],it->y[1],it->z[1],grid,data);      // Determine length of new segment                   
		seg.length = data.vess_length;

	
		double den_scale = 1.0;
		den_scale = findDenScale(it->x[1], it->y[1], it->z[1], grid);
			
		//seg.phi1 = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data,1);    // Determine the angle between the new segment and the x-axis
		//seg.phi2 = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data,2);    // Determine the angle between the new segment and the x-y plane
		seg.uvect = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data); 
		
		seg.length = den_scale*data.vess_length;

		//if (data.branch)                                                   // If new segment is a branch...
		//{
		//	if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
		//			seg.phi1 = ang_dom(it->phi1+pi/2);
		//        else
		//	        seg.phi1 = ang_dom(it->phi1-pi/2);
		//	
		//	//if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
		//	//		seg.phi2 = ang_dom(it->phi2+pi/2);
		// //       else
		//	//        seg.phi2 = ang_dom(it->phi2-pi/2);
		//	
		//	//seg.phi2 = 0.1*seg.phi2;                                           // force the branch to be mostly planar (x-y plane)
		//}
		//
		
		if (data.branch)                                                   // If new segment is a branch...
		{
			//double x = seg.uvect.x;
			//double y = seg.uvect.y;
			//double H = sqrt(x*x + y*y);
			//double theta = acos(x/H);

			//if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
			//		theta = theta + pi/2;
		 //       else
			//        theta = theta - pi/2;
			//
			//seg.uvect.x = H*cos(theta);
			//seg.uvect.y = H*sin(theta);
			////seg.uvect.z = 0;

			//seg.uvect = seg.uvect/seg.uvect.norm();
		
			//vect3 coll_fib = findCollAngle(it->x[1], it->y[1], it->z[1], grid, data);

			//if (acos(coll_fib*seg.uvect) > pi/2)
			//	coll_fib = -coll_fib;

			//seg.uvect = (seg.uvect + coll_fib)/2;
			
			vect3 coll_fib = findCollAngle(it->x[1], it->y[1], it->z[1], grid, data);
			vect3 newseg;
			newseg = coll_fib - seg.uvect*(seg.uvect*coll_fib)*0.5;
			newseg = newseg/newseg.norm();
			seg.uvect = newseg;

			data.num_lines++;
			//seg.line_num = data.num_lines;
		}

		seg.label = it->label;                                  // Transfer the label of the parent segment to the new segment
		seg.vessel = it->vessel;                                // Transfer the vessel number of the parent segment to the new segment	
		
		seg.x[0] = it->x[1];                                    // Set the origin of new segment as the active tip of the previous segment
		seg.y[0] = it->y[1];
		seg.z[0] = it->z[1];
		seg.tip_elem[0] = it->tip_elem[1];
		seg.seg_conn[0][0] = it->seg_num;

		//seg.x[1] = seg.x[0]+seg.length*cos(seg.phi2)*cos(seg.phi1);     // Determine the x-coordinate of the end point using the length vector and orientation angles
		//seg.y[1] = seg.y[0]+seg.length*cos(seg.phi2)*sin(seg.phi1);     // Determine the y-coordinate of the end point using the length vector and orientation angles
		//seg.z[1] = seg.z[0]+seg.length*sin(seg.phi2);                   // Determine the z-coordinate of the end point using the length vector and orientation angles
    
		seg.x[1] = seg.x[0] + seg.length*seg.uvect.x;					  // Determine the x-coordinate of the end point using the length vector and orientation angles
		seg.y[1] = seg.y[0] + seg.length*seg.uvect.y;					  // Determine the y-coordinate of the end point using the length vector and orientation angles
		seg.z[1] = seg.z[0] + seg.length*seg.uvect.z;	                  // Determine the z-coordinate of the end point using the length vector and orientation angles
		seg.findlength();

        seg.tip[1] = 1;                                         // Turn on end tip of new segment
		seg.tip[0] = 0;                                         // Turn off origin tip of new segment
		
		seg.TofBirth = data.t;                                  // Stamp segment with time of birth
		it->tip[k] = 0;                                         // Turn off previous segment tip
		seg.bdyf_id[k] = it->bdyf_id[k];
        
        seg.sprout = 1;                                         // Set sprout for the new segment as +1, indicating this segment originated from a +1 tip
					    
		if (it->seg_conn[1][0] == 0)
			it->seg_conn[1][0] = seg.seg_num;
		else
			it->seg_conn[1][1] = seg.seg_num;
				
		elem_num = grid.findelem(seg.x[1], seg.y[1], seg.z[1]);
		
		if (elem_num == -1)
			bc.checkBC(seg, 1, grid, data, frag);
		else
			seg.tip_elem[1] = elem_num;

	}
	
	
	else if (it->tip[k] == -1)                                    // If the parent vessel active tip is set as +1...
	{
	    //seg.length = -1*findLength(it->x[0],it->y[0],it->z[0],grid,data);  // Determine length of new segment 
		seg.length = -data.vess_length;
		
		double den_scale = 1.0;
		den_scale = findDenScale(it->x[0], it->y[0], it->z[0], grid);
			
		//seg.phi1 = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data,1);    // Determine the angle between the new segment and the x-axis
		//seg.phi2 = findAngle(it,it->x[1],it->y[1],it->z[1],grid,data,2);    // Determine the angle between the new segment and the x-y plane
		
		seg.uvect = findAngle(it,it->x[0],it->y[0],it->z[0],grid,data);

		seg.length = -den_scale*data.vess_length;
				
				
		//if (data.branch)                                                   // If new segment is a branch...
		//{
		//	if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
		//			seg.phi1 = ang_dom(it->phi1+pi/2);
		//        else
		//	        seg.phi1 = ang_dom(it->phi1-pi/2);
		//	
		//	//if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
		//	//		seg.phi2 = ang_dom(it->phi2+pi/2);
		// //       else
		//	//        seg.phi2 = ang_dom(it->phi2-pi/2);
		//	
		//	//seg.phi2 = 0.1*seg.phi2;                                           // force the branch to be mostly planar (x-y plane)
		//}

		if (data.branch)                                                   // If new segment is a branch...
		{
			//double x = seg.uvect.x;
			//double y = seg.uvect.y;
			//double H = sqrt(x*x + y*y);
			//double theta = acos(x/H);

			//if (float(rand())/RAND_MAX < 0.5)									// force the phi1 to be perpindicular with the parent segment                       
			//		theta = theta + pi/2;
		 //       else
			//        theta = theta - pi/2;
			//
			//seg.uvect.x = H*cos(theta);
			//seg.uvect.y = H*sin(theta);
			////seg.uvect.z = 0;

			//seg.uvect = seg.uvect/seg.uvect.norm();
		
			//vect3 coll_fib = findCollAngle(it->x[0], it->y[0], it->z[0], grid, data);

			//if (acos(coll_fib*seg.uvect) > pi/2)
			//	coll_fib = -coll_fib;

			//seg.uvect = (seg.uvect + coll_fib)/2;

			vect3 coll_fib = findCollAngle(it->x[1], it->y[1], it->z[1], grid, data);
			vect3 newseg;
			newseg = coll_fib - seg.uvect*(seg.uvect*coll_fib)*0.5;
			newseg = newseg/newseg.norm();

			seg.uvect = newseg;
			
			data.num_lines++;
			//seg.line_num = data.num_lines;
		}

		seg.label = it->label;                                  // Transfer the label of the parent segment to the new segment
		seg.vessel = it->vessel;                                // Transfer the vessel number of the parent segment to the new segment
		
		seg.x[1] = it->x[0];                                    // Set the origin of new segment as the active tip of the previous segment
		seg.y[1] = it->y[0];
		seg.z[1] = it->z[0];
		seg.tip_elem[1] = it->tip_elem[0];
		seg.seg_conn[1][0] = it->seg_num;
		
		//seg.x[0] = seg.x[1]+seg.length*cos(seg.phi2)*cos(seg.phi1);     // Determine the x-coordinate of the end point using the length vector and orientation angles
		//seg.y[0] = seg.y[1]+seg.length*cos(seg.phi2)*sin(seg.phi1);     // Determine the y-coordinate of the end point using the length vector and orientation angles
		//seg.z[0] = seg.z[1]+seg.length*sin(seg.phi2);                   // Determine the z-coordinate of the end point using the length vector and orientation angles
		
		seg.x[0] = seg.x[1] + seg.length*seg.uvect.x;					  // Determine the x-coordinate of the end point using the length vector and orientation angles
		seg.y[0] = seg.y[1] + seg.length*seg.uvect.y;					  // Determine the y-coordinate of the end point using the length vector and orientation angles
		seg.z[0] = seg.z[1] + seg.length*seg.uvect.z;	                  // Determine the z-coordinate of the end point using the length vector and orientation angles
		seg.findlength();

		seg.tip[0] = -1;                                        // Turn on end tip of new segment
		seg.tip[1] = 0;                                         // Turn off origin tip of new segment
		
		seg.TofBirth = data.t;                                  // Stamp segment with time of birth
		it->tip[k] = 0;                                         // Turn off previous segment tip
		seg.bdyf_id[k] = it->bdyf_id[k];

        seg.sprout = -1;                                        // Set sprout for the new segment as -1, indicating this segment originated from a -1 tip

		if (it->seg_conn[0][0] == 0)
			it->seg_conn[0][0] = seg.seg_num;
		else
			it->seg_conn[0][1] = seg.seg_num;
		
		
		elem_num = grid.findelem(seg.x[0], seg.y[0], seg.z[0]);
		
		if (elem_num == -1)
			bc.checkBC(seg, 0, grid, data, frag);
		else
			seg.tip_elem[0] = elem_num;
	
		//if (it->sprout == 9){
		//	data.num_lines++;
		//	seg.line_num = data.num_lines;}
	}
	
	seg.x0[0] = seg.x[0];
	seg.x0[1] = seg.x[1];
	seg.y0[0] = seg.y[0];
	seg.y0[1] = seg.y[1];
	seg.z0[0] = seg.z[0];
	seg.z0[1] = seg.z[1];

	//if (it->Recent_branch == 1){
	//	data.num_lines++;
	//	seg.line_num = data.num_lines;}
		
	return seg;                                                 // Return the new segment 
}



///////////////////////////////////////////////////////////////////////
// findLength
///////////////////////////////////////////////////////////////////////
// CULTURE.findLength - Determine the length of a newly created segment based on the growth function g(t)
//      Input:  - Coordinates of the active tip that is sprouting the new segment (xpt, ypt, zpt)
//              - GRID object
//              - DATA object
//
//      Output: - Magnitude of length vector (in um)

//double Culture::findLength(double xpt, double ypt, double zpt, Grid &grid, Data &data)
//{
//	double lc;                                                  // lc - Length calculation obtained from growth function g(t)
//	double nt;                                                  // nt - Number of active tips
//    double length;
//    
//	lc = data.a/(1.+pow(E,-(data.t-data.x0)/data.b));
//	lc -= data.a/(1.+pow(E,-(data.t-data.dt-data.x0)/data.b));
//	
//	nt = double(2*data.NFRAGS + data.num_branches - data.num_anastom);
//	//nt = double(2*data.NFRAGS + data.num_branches - data.num_anastom - data.num_zdead);
//	
//	if (nt <= 0)
//		//nt = double(2*data.NFRAGS + data.num_branches - data.num_anastom - data.num_zdead);
//		nt = 4*data.NFRAGS;
//	
//	length = (data.NFRAGS*lc/nt)*grid.den_scale*data.length_adjust;
//	
//	return ;
//}



///////////////////////////////////////////////////////////////////////
// findAngle
///////////////////////////////////////////////////////////////////////

// CULTURE.findAngle - Determine the orientation angle of a newly created segment
//      Input:  - Iterator for segment container which points to the parent segment (it)
//              - Coordinates of the active tip that is sprouting the new segment (xpt, ypt, zpt)
//              - GRID object
//              - DATA object
//              - Integer describing the angle type: 1 for phi1, 2 for phi2 (ang_type)
//
//      Output: - Segment orientation angle phi1 or phi1 (in radians)

vect3 Culture::findAngle(list<Segment>::iterator it, double xpt, double ypt, double zpt, Grid &grid,Data &data)
{
	vect3 angle;

	double den_scale = 1.0;
	
	den_scale = findDenScale(xpt, ypt, zpt, grid);

	//double W[4] = {10*(1/grid.den_scale), 0, 0, 100};
	//double W[4] = {10, 0, 0, 100};                       // W[0] = Weight for collagen orientation
	//double W[4] = {10, 0, 0, 50};                                                                        // W[1] = Weight for vessel density
	                                                                        // W[2] = Weight for random component
	                                                                        // W[3] = Weight for previous vessel direction
	                                        
    
    vect3 coll_angle;      // Component of new vessel orientation resulting from collagen fiber orientation
    vect3 den_angle;       // Component of new vessel orientation resulting from vessel density
    vect3 ran_angle;       // Component of new vessel orientation resulting from random walk
    vect3 per_angle;       // Component of new vessel orientation resulting from previous vessel direction        
        
    
    // Find the component of the new vessel direction determined by collagen fiber orientation    
    coll_angle = findCollAngle(xpt, ypt, zpt, grid, data);

	per_angle = it->uvect;

	angle = (coll_angle*W[0] + per_angle*W[3])/(W[0]+W[3]);

	angle = angle/angle.norm();

	return angle;
	
}



///////////////////////////////////////////////////////////////////////
// findCollAngle
///////////////////////////////////////////////////////////////////////


vect3 Culture::findCollAngle(double xpt, double ypt, double zpt, Grid &grid, Data &data)
{       
    //double coll_angle = 0.0;
    vect3 coll_angle;

    double xix, xiy, xiz;
    double shapeF[8];
    
    int elem_num = grid.findelem(xpt, ypt, zpt);
    
    if (elem_num < 0)
		return coll_angle;
		
	Elem elem;
    elem = grid.ebin[elem_num];
        
    // Convert to natural coordinates -1 <= Xi <= +1
    grid.natcoord(xix, xiy, xiz, xpt, ypt, zpt, elem_num);        
    
    // Obtain shape function weights
    grid.shapefunctions(shapeF, xix, xiy, xiz);
    
    // Determine component of new vessel direction due to nodal collagen fiber orientation
	coll_angle = ((*elem.n1).collfib)*shapeF[0] + ((*elem.n2).collfib)*shapeF[1] + ((*elem.n3).collfib)*shapeF[2] + ((*elem.n4).collfib)*shapeF[3] + ((*elem.n5).collfib)*shapeF[4] + ((*elem.n6).collfib)*shapeF[5] + ((*elem.n7).collfib)*shapeF[6] + ((*elem.n8).collfib)*shapeF[7];
	
	coll_angle = coll_angle/coll_angle.norm();

    return coll_angle;
}



///////////////////////////////////////////////////////////////////////
// findDenScale
///////////////////////////////////////////////////////////////////////

double Culture::findDenScale(double xpt, double ypt, double zpt, Grid &grid)
{
	double coll_den = 0.0;
    double den_scale = 1.0;

    double xix, xiy, xiz;
    double shapeF[8];

	int elem_num = grid.findelem(xpt, ypt, zpt);
    
    if (elem_num < 0)
		return den_scale;
		
	Elem elem;
    elem = grid.ebin[elem_num];
        
    // Convert to natural coordinates -1 <= Xi <= +1
    grid.natcoord(xix, xiy, xiz, xpt, ypt, zpt, elem_num);        
    
    // Obtain shape function weights
    grid.shapefunctions(shapeF, xix, xiy, xiz);
    
	coll_den = shapeF[0]*(*elem.n1).ecm_den + shapeF[1]*(*elem.n2).ecm_den + shapeF[2]*(*elem.n3).ecm_den + shapeF[3]*(*elem.n4).ecm_den + shapeF[4]*(*elem.n5).ecm_den + shapeF[5]*(*elem.n6).ecm_den + shapeF[6]*(*elem.n7).ecm_den + shapeF[7]*(*elem.n8).ecm_den;
    
	den_scale = grid.find_density_scale(coll_den);

    if (den_scale < 0.)
		den_scale = 0.;
	
	return den_scale;
}



///////////////////////////////////////////////////////////////////////
// findDenAngle
///////////////////////////////////////////////////////////////////////

//double Culture::findDenAngle(double xpt, double ypt, double zpt, Grid &grid, Data &data, int ang_type)
//{
//    double den_angle=0;
//    double angle_det, grad_vec[3], x_vec[3], z_vec[3], dens_ang, rand_ang;
//    
//    int I, J, K;
//	grid.ijk(xpt, ypt, zpt, I, J, K, data);                   // Determine which cell the tip is in
//    
//    double xix, xiy, xiz;
//    double shapeF[8];
//      
//    // Convert to natural coordinates -1 <= Xi <= +1
//    grid.natcoordinates(xix, xiy, xiz, xpt, ypt, zpt, I, J, K);        
//        
//    // Obtain shape function weights
//    grid.shapefunctions(shapeF, xix, xiy, xiz);
//
//	double dens_dot;
//	
//	grad_vec[0] = 0.;
//	grad_vec[1] = 0.;
//	grad_vec[2] = 0.;
//	
//	x_vec[0] = 1.;
//	x_vec[1] = 0.;
//	x_vec[2] = 0.;
//	
//	z_vec[0] = 0.;
//	z_vec[1] = 0.;
//	z_vec[2] = 1.;
//	
//	
//	if (ang_type == 1)
//	{
//		grad_vec[0] = -(shapeF[0]*grid.drho_dx[I][J][K] + shapeF[1]*grid.drho_dx[I][J+1][K] + shapeF[2]*grid.drho_dx[I+1][J][K] + shapeF[3]*grid.drho_dx[I+1][J+1][K]);
//		grad_vec[0] += -(shapeF[4]*grid.drho_dx[I][J][K+1] + shapeF[5]*grid.drho_dx[I][J+1][K+1] + shapeF[6]*grid.drho_dx[I+1][J][K+1] + shapeF[7]*grid.drho_dx[I+1][J+1][K+1]);
//
//		grad_vec[1] = -(shapeF[0]*grid.drho_dy[I][J][K] + shapeF[1]*grid.drho_dy[I][J+1][K] + shapeF[2]*grid.drho_dy[I+1][J][K] + shapeF[3]*grid.drho_dy[I+1][J+1][K]);
//		grad_vec[1] += -(shapeF[4]*grid.drho_dy[I][J][K+1] + shapeF[5]*grid.drho_dy[I][J+1][K+1] + shapeF[6]*grid.drho_dy[I+1][J][K+1] + shapeF[7]*grid.drho_dy[I+1][J+1][K+1]);		
//		
//		grad_vec[2] = 0;
//		
//		grad_vec[0] = -grad_vec[0];
//		grad_vec[1] = -grad_vec[1];
//		
//		//determine angle of density gradient
//		if (vec_norm(grad_vec) != 0.0)
//		{
//			dens_dot = vec_dot(grad_vec,x_vec)/vec_norm(grad_vec);
//			den_angle = ang_dom(acos(dens_dot));
//		    //cout << endl << "dens_dot " << dens_dot << "   den_angle " << den_angle << endl;
//		}
//		else
//		    den_angle = 0;
//		
//		if (grad_vec[1] < 0)
//		    den_angle = -den_angle;
//    }
//	
//	if (ang_type == 2)
//	{
//		grad_vec[0] = -(shapeF[0]*grid.drho_dx[I][J][K] + shapeF[1]*grid.drho_dx[I][J+1][K] + shapeF[2]*grid.drho_dx[I+1][J][K] + shapeF[3]*grid.drho_dx[I+1][J+1][K]);
//		grad_vec[0] += -(shapeF[4]*grid.drho_dx[I][J][K+1] + shapeF[5]*grid.drho_dx[I][J+1][K+1] + shapeF[6]*grid.drho_dx[I+1][J][K+1] + shapeF[7]*grid.drho_dx[I+1][J+1][K+1]);
//
//		grad_vec[1] = -(shapeF[0]*grid.drho_dy[I][J][K] + shapeF[1]*grid.drho_dy[I][J+1][K] + shapeF[2]*grid.drho_dy[I+1][J][K] + shapeF[3]*grid.drho_dy[I+1][J+1][K]);
//		grad_vec[1] += -(shapeF[4]*grid.drho_dy[I][J][K+1] + shapeF[5]*grid.drho_dy[I][J+1][K+1] + shapeF[6]*grid.drho_dy[I+1][J][K+1] + shapeF[7]*grid.drho_dy[I+1][J+1][K+1]);		
//				
//		grad_vec[2] = -(shapeF[0]*grid.drho_dz[I][J][K] + shapeF[1]*grid.drho_dz[I][J+1][K] + shapeF[2]*grid.drho_dz[I+1][J][K] + shapeF[3]*grid.drho_dz[I+1][J+1][K]);
//		grad_vec[2] += -(shapeF[4]*grid.drho_dz[I][J][K+1] + shapeF[5]*grid.drho_dz[I][J+1][K+1] + shapeF[6]*grid.drho_dz[I+1][J][K+1] + shapeF[7]*grid.drho_dz[I+1][J+1][K+1]);		
//		
//		grad_vec[1] = -grad_vec[1];
//		grad_vec[2] = -grad_vec[2];
//					
//		//determine angle of density gradient
//		if (vec_norm(grad_vec) != 0.0)
//		{
//			dens_dot = vec_dot(grad_vec,z_vec)/vec_norm(grad_vec);
//			den_angle = ang_dom(acos(dens_dot));        // Angle between gradient vector and normal to x-y place (z-axis)
//			//cout << endl << "dens_dot " << dens_dot << "   den_angle " << den_angle << endl;
//	        den_angle = (0.5*pi) - 0.5*den_angle;       // Angle between gradient vector and x-y plane	
//		}
//		else
//			den_angle = 0;
//		
//		if (grad_vec[2] < 0)
//		    den_angle = -den_angle;
//	}
//
//    return den_angle;
//}


///////////////////////////////////////////////////////////////////////
// connectSegment
///////////////////////////////////////////////////////////////////////

// creates a new segment to connect close segments

Segment Culture::connectSegment(list<Segment>::iterator it,list<Segment>::iterator it2, int k, int kk, Grid &grid, Data &data, list<Segment> &frag)
 {
 	Segment seg;
 	double oppadj;
 	
 	seg.length = sqrt((it->x[k]-it2->x[kk])*(it->x[k]-it2->x[kk])+
 		(it->y[k]-it2->y[kk])*(it->y[k]-it2->y[kk])+(it->z[k]-it2->z[kk])*(it->z[k]-it2->z[kk]));
 
 	/*if (it->x[k] != it2->x[kk])
 	{
 		oppadj = (it2->y[kk]-it->y[k])/(it2->x[kk]-it->x[k]);
 		seg.phi1 = atan(oppadj);
		oppadj = (it2->z[kk]-it->z[k])/(it2->x[kk]-it->x[k]);
		seg.phi2 = atan(oppadj);
 	}
 	else
	{
 		seg.phi1 = pi/2;
		seg.phi2 = pi/2;
	}*/
		
 
 	seg.x[0] = it->x[k];
 	seg.y[0] = it->y[k];
	seg.z[0] = it->z[k];
	seg.tip_elem[0] = it->tip_elem[k];
	seg.x[1] = it2->x[kk];
 	seg.y[1] = it2->y[kk];
	seg.z[1] = it2->z[kk];
	seg.tip_elem[1] = it2->tip_elem[kk];
 	
	seg.x0[0] = seg.x[0];
	seg.x0[1] = seg.x[1];
	seg.y0[0] = seg.y[0];
	seg.y0[1] = seg.y[1];
	seg.z0[0] = seg.z[0];
	seg.z0[1] = seg.z[1];

 	seg.TofBirth = data.t;
 	seg.label = it->label;
 	seg.vessel = it->vessel;
 	
 	seg.tip[0] = 0;
 	seg.tip[1] = 0;
 	seg.anast = 1;
 	it->anast = 1;
 	it2->anast = 1;
	seg.findlength();

 	++data.num_anastom;
 	
	++data.nsegs;
	seg.seg_num = data.nsegs;

	seg.seg_conn[0][0] = it->seg_num;
	
	if (k == 0)
		it->seg_conn[0][0] = seg.seg_num;
	else
		it->seg_conn[1][0] = seg.seg_num;

	seg.seg_conn[1][0] = it2->seg_num;

	if (kk == 0)
		it2->seg_conn[0][1] = seg.seg_num;
	else
		it2->seg_conn[1][1] = seg.seg_num;
	
	return seg;
 }
 
///////////////////////////////////////////////////////////////////////
// CheckForIntersection
///////////////////////////////////////////////////////////////////////
// Description: checks for intersection between a passed segment and all other existing segments that are not members
// of vessel containing the segment

void Culture::CheckForIntersection(Segment &seg,list<Segment> &frag, Data &data, list<Segment>::iterator it)
{
	double p1[3], p2[3], pp1[3], pp2[3]; //tip points for segment to check intersection
	list<Segment>::iterator it2; //loop over segments in list
	double intersectpt[3]; //the intersection pt of vectors
	
	intersectpt[0] = 0;
	intersectpt[1] = 0;
	
	double lambda;
	
	if (seg.tip[1] == 1)
	{
		p1[0] = seg.x[0];
		p2[0] = seg.x[1];
		p1[1] = seg.y[0];
		p2[1] = seg.y[1];
		p1[2] = seg.z[0];
		p2[2] = seg.z[1];
	}
	else if (seg.tip[0] == -1)
	{
		p1[0] = seg.x[1];
		p2[0] = seg.x[0];
		p1[1] = seg.y[1];
		p2[1] = seg.y[0];
		p1[2] = seg.z[1];
		p2[2] = seg.z[0];
	}

	for (it2 = frag.begin(); it2 != frag.end(); ++it2)
	{
		if (it->label!=it2->label)
		{
			/*
			if (it->length > 0)
			{
				pp1[0] = it2->x[0];
				pp1[1] = it2->y[0];
				pp2[0] = it2->x[1];
				pp2[1] = it2->y[1];
			}
			else
			{
				pp1[0] = it2->x[1];
				pp1[1] = it2->y[1];
				pp2[0] = it2->x[0];
				pp2[1] = it2->y[0];
			}*/
			pp1[0] = it2->x[0];
			pp1[1] = it2->y[0];
			pp1[2] = it2->z[0];
			pp2[0] = it2->x[1];
			pp2[1] = it2->y[1];
			pp2[2] = it2->z[1];

			lambda = findIntersect(p1,p2,pp1,pp2,intersectpt);
			if (lambda >=0 && lambda <=1)
			{
				if (seg.tip[1] == 1)
				{
					seg.x[1] = intersectpt[0];
					seg.y[1] = intersectpt[1];
					seg.z[1] = intersectpt[2];
					seg.tip[1] = 0;
				}
				else if (seg.tip[0] == -1)
				{
					seg.x[0] = intersectpt[0];
					seg.y[0] = intersectpt[1];
					seg.z[0] = intersectpt[2];
					seg.tip[0] = 0;
				}
				cout << "3D intersection" << endl;
				++data.num_anastom;
			}
				
		}
	}


	return;
}


///////////////////////////////////////////////////////////////////////
// findIntersect
///////////////////////////////////////////////////////////////////////

// variables:
// a : start x,y point on segment1
// b : end x,y point on segment1
// c : start x,y point on segment2
// d : end x,y point on segment2
// t1 : displacement vector for segment1 (b-a)
// t2 : displacement vector for segment2 (d-c)
// lambda : distance to move along t1 from its start point
// mu : distance to move along t2 from its start point

// General case for the minimum distance between two lines
// d^2 = [BCD + B^2E + C^2F + A(D^2-4EF)]/[C^2 - 4AE]
// where
// A = t1.t1
// B = 2(a.t1 - t1.c)
// C = 2(t1.t2)
// D = 2(t2.c - t2.a)
// E = t2.t2
// F = a.a + c.c
// This will have a unique solution if
// | 2A -C |
// | -C	2E | is non-zero
// closest point on t1 is : a + lambda*t1
// where
// lambda = (C*mu - B)/(2A)
// and
// mu = (2*A*D + B*C)/(C^2 - 4*A*E) 


double Culture::findIntersect(double a[3], double b[3], double c[3], double d[3], double intersectpt[3])
{
	double t1[3], t2[3];
	double A, B, C, D, E, F, min_dist;
	double lambda, mu;
	
	t1[0] = b[0] - a[0];
	t1[1] = b[1] - a[1];
	t1[2] = b[2] - a[2];

	t2[0] = d[0] - c[0];
	t2[1] = d[1] - c[1];
	t2[2] = d[2] - c[2];

	A = vec_dot(t1,t1);
	B = 2*(vec_dot(a,t1) - vec_dot(t1,c));
	C = 2*(vec_dot(t1,t2));
	D = 2*(vec_dot(t2,c) - vec_dot(t2,a));
	E = vec_dot(t2,t2);
	F = vec_dot(a,a)+ vec_dot(c,c);

	if (4*A*E-C*C != 0)
	{
		min_dist = (B*C*D + B*B*E + C*C*F + A*(D*D-4*E*F))/(C*C - 4*A*E);
		if (min_dist <= 7)
		{
			mu = (2*A*D + B*C)/(C*C - 4*A*E);
			lambda = (C*mu - B)/(2*A);
			intersectpt[0] = a[0] + lambda*t1[0];
			intersectpt[1] = a[1] + lambda*t1[1];
			intersectpt[2] = a[2] + lambda*t1[2];
		}
		else
		{
			return 1000;
		}
	}
	else
	{
		return 1000;
	}

	return lambda;
}


///////////////////////////////////////////////////////////////////////
// intersectPlane
///////////////////////////////////////////////////////////////////////

// line equation P = LP[3] + u*V[3]
// Plane equation N[3].(P[3]-P*[3]) = 0
// solving for u, u={N[3].(P*[3]-LP[3])}/{N[3].V[3]}

// box face numbering
// Front = 0
// Right = 1
// Back = 2+
// Left = 3
// Top = 4
// Bottom = 5


double N0[3], N1[3], N2[3], N3[3], N4[3], N5[3]; //face normals
double P0[3], P1[3], P2[3], P3[3], P4[3], P5[3]; //point on faces


bool Culture::intersectPlane(Grid& grid, Segment &seg, int n, double intersectpt[3])
{

	//front
	N0[0] = 0;
	N0[1] = -1;
	N0[2] = 0;
	P0[0] = 0;
	P0[1] = 0;
	P0[2] = 0;

	//right
	N1[0] = 1;
	N1[1] = 0;
	N1[2] = 0;
	P1[0] = grid.xrange[1];
	P1[1] = 0;
	P1[2] = 0;

	//back
	N2[0] = 0;
	N2[1] = 1;
	N2[2] = 0;
	P2[0] = 0;
	P2[1] = grid.yrange[1];
	P2[2] = 0;

	//left
	N3[0] = -1;
	N3[1] = 0;
	N3[2] = 0;
	P3[0] = 0;
	P3[1] = 0;
	P3[2] = 0;

	//top
	N4[0] = 0;
	N4[1] = 0;
	N4[2] = 1;
	P4[0] = 0;
	P4[1] = 0;
	P4[2] = grid.zrange[1];

	//bottom
	N5[0] = 0;
	N5[1] = 0;
	N5[2] = -1;
	P5[0] = 0;
	P5[1] = 0;
	P5[2] = 0;
	double V[3]; //segment displacement vector
	double V2[3]; //P*[3]-LP[3]
	double LP[3]; //origin of segment
	double u; //scalar weight to move along segment displacement vector


	if (seg.tip[1] == 1)
	{
		V[0] = seg.x[1] - seg.x[0];
		V[1] = seg.y[1] - seg.y[0];
		V[2] = seg.z[1] - seg.z[0];
		LP[0] = seg.x[0];
		LP[1] = seg.y[0];
		LP[2] = seg.z[0];
	}
	else
	{
		V[0] = seg.x[0] - seg.x[1];
		V[1] = seg.y[0] - seg.y[1];
		V[2] = seg.z[0] - seg.z[1];
		LP[0] = seg.x[1];
		LP[1] = seg.y[1];
		LP[2] = seg.z[1];
	}

	switch (n)
	{
		case 0: //front
			V2[0] = P0[0] - LP[0];
			V2[1] = P0[1] - LP[1];
			V2[2] = P0[2] - LP[2];
			u = vec_dot(N0,V2)/vec_dot(N0,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[0] > grid.xrange[0] && intersectpt[0] < grid.xrange[1] &&
					intersectpt[2] > grid.zrange[0] && intersectpt[2] < grid.zrange[1])
					return true;
			}
			break;
		case 1: //right
			V2[0] = P1[0] - LP[0];
			V2[1] = P1[1] - LP[1];
			V2[2] = P1[2] - LP[2];
			u = vec_dot(N1,V2)/vec_dot(N1,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[1] > grid.yrange[0] && intersectpt[1] < grid.yrange[1] &&
					intersectpt[2] > grid.zrange[0] && intersectpt[2] < grid.zrange[1])
					return true;
			}
			break;
		case 2: //back
			V2[0] = P2[0] - LP[0];
			V2[1] = P2[1] - LP[1];
			V2[2] = P2[2] - LP[2];
			u = vec_dot(N2,V2)/vec_dot(N2,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[0] > grid.xrange[0] && intersectpt[0] < grid.xrange[1] &&
					intersectpt[2] > grid.zrange[0] && intersectpt[2] < grid.zrange[1])
					return true;
			}
			break;
		case 3: //left
			V2[0] = P3[0] - LP[0];
			V2[1] = P3[1] - LP[1];
			V2[2] = P3[2] - LP[2];
			u = vec_dot(N3,V2)/vec_dot(N3,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[1] > grid.yrange[0] && intersectpt[1] < grid.yrange[1] &&
					intersectpt[2] > grid.zrange[0] && intersectpt[2] < grid.zrange[1])
					return true;
			}
			break;
		case 4: //top
			V2[0] = P4[0] - LP[0];
			V2[1] = P4[1] - LP[1];
			V2[2] = P4[2] - LP[2];
			u = vec_dot(N4,V2)/vec_dot(N4,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[0] > grid.xrange[0] && intersectpt[0] < grid.xrange[1] &&
					intersectpt[1] > grid.yrange[0] && intersectpt[1] < grid.yrange[1])
					return true;
			}
			break;
		case 5: //bottom
			V2[0] = P5[0] - LP[0];
			V2[1] = P5[1] - LP[1];
			V2[2] = P5[2] - LP[2];
			u = vec_dot(N5,V2)/vec_dot(N5,V);
			if (u>0 && u<1)
			{
				intersectpt[0] = LP[0] + u*V[0];
				intersectpt[1] = LP[1] + u*V[1];
				intersectpt[2] = LP[2] + u*V[2];
				if (intersectpt[0] > grid.xrange[0] && intersectpt[0] < grid.xrange[1] &&
					intersectpt[1] > grid.yrange[0] && intersectpt[1] < grid.yrange[1])
					return true;
			}
			break;
	}
	return false;
}








///////////////////////////////////////////////////////////////////////
// PeriodicBC
///////////////////////////////////////////////////////////////////////

//table of faces opposite to 0,1,..6
double oppface[6];

Segment Culture::PeriodicBC(Segment &seg,Grid &grid,list<Segment> &frag,Data &data)
{

	oppface[0] = grid.yrange[1];
	oppface[1] = grid.xrange[0];
	oppface[2] = grid.yrange[0];
	oppface[3] = grid.xrange[1];
	oppface[4] = grid.zrange[0];
	oppface[5] = grid.zrange[1];
	double unit_vec[3] = {0};
	double length = 0.0;
	double rem_length = 0.0;
	int n = 0;
	double intersectpt[3] = {0};
	
	if (seg.tip[1] == 1)
	{
		unit_vec[0] = (seg.x[1]-seg.x[0]);
		unit_vec[1] = (seg.y[1]-seg.y[0]);
		unit_vec[2] = (seg.z[1]-seg.z[0]);
		length = vec_norm(unit_vec);
		unit_vec[0] /= length;
		unit_vec[1] /= length;
		unit_vec[2] /= length;
		
		for (n=0;n<6;++n)
		{
			if (intersectPlane(grid, seg,n,intersectpt))
			{
				//cout << endl << "Vessel " << seg.label << " Crossed Face " << n << endl;
				seg.x[1] = intersectpt[0];
				seg.y[1] = intersectpt[1];
				seg.z[1] = intersectpt[2];
				if (seg.length < 0)
				{
					seg.length = -sqrt((seg.x[1]-seg.x[0])*(seg.x[1]-seg.x[0]) +
					(seg.y[1]-seg.y[0])*(seg.y[1]-seg.y[0]) + (seg.z[1]-seg.z[0])*(seg.z[1]-seg.z[0]));
				}
				else
				{
					seg.length = sqrt((seg.x[1]-seg.x[0])*(seg.x[1]-seg.x[0]) +
					(seg.y[1]-seg.y[0])*(seg.y[1]-seg.y[0]) + (seg.z[1]-seg.z[0])*(seg.z[1]-seg.z[0]));
				}

				seg.tip[1] = 0;
				seg.BCdead = 1;
				
				/*if (n > 3){
				    data.num_zdead = data.num_zdead+1;
				    return seg;}*/
				
				if (n > 3)
				    data.num_zdead = data.num_zdead+1;
				return seg;
				
				frag.push_front (seg);
				Segment seg2;
				data.num_vessel = data.num_vessel + 1;
				
				switch (n)
				{
				case 0:
					seg2.x[0] = seg.x[1];
					seg2.y[0] = oppface[n];
					seg2.z[0] = seg.z[1];
					break;
					
				case 1:
					seg2.x[0] = oppface[n];
					seg2.y[0] = seg.y[1];
					seg2.z[0] = seg.z[1];
					break;
				case 2:
					seg2.x[0] = seg.x[1];
					seg2.y[0] = oppface[n];
					seg2.z[0] = seg.z[1];
					break;
				case 3:
					seg2.x[0] = oppface[n];
					seg2.y[0] = seg.y[1];
					seg2.z[0] = seg.z[1];
					break;
				case 4:
					seg2.x[0] = seg.x[1];
					seg2.y[0] = seg.y[1];
					seg2.z[0] = oppface[n];
					break;
				case 5:
					seg2.x[0] = seg.x[1];
					seg2.y[0] = seg.y[1];
					seg2.z[0] = oppface[n];
					break;
				}

				if (seg.length > 0) 
				{
					rem_length = length - sqrt((seg.x[1]-seg.x[0])*(seg.x[1]-seg.x[0])
					+ (seg.y[1]-seg.y[0])*(seg.y[1]-seg.y[0])
					+ (seg.z[1]-seg.z[0])*(seg.z[1]-seg.z[0]));
					seg2.x[1] = seg2.x[0] + rem_length*unit_vec[0];
					seg2.y[1] = seg2.y[0] + rem_length*unit_vec[1];
					seg2.z[1] = seg2.z[0] + rem_length*unit_vec[2];
				}
				else 
				{
					rem_length = -length + sqrt((seg.x[1]-seg.x[0])*(seg.x[1]-seg.x[0])
					+ (seg.y[1]-seg.y[0])*(seg.y[1]-seg.y[0])
					+ (seg.z[1]-seg.z[0])*(seg.z[1]-seg.z[0]));
					seg2.x[1] = seg2.x[0] - rem_length*unit_vec[0];
					seg2.y[1] = seg2.y[0] - rem_length*unit_vec[1];
					seg2.z[1] = seg2.z[0] - rem_length*unit_vec[2];
				}
				
				seg2.tip[1] = 1;
				seg2.tip[0] = 0;
				seg2.label = seg.label;
				seg2.vessel = data.num_vessel;
				seg2.sprout = seg.sprout;
				seg2.length = rem_length;
//				seg2.phi1 = seg.phi1;
//				seg2.phi2 = seg.phi2;
				seg2.Recent_branch = seg.Recent_branch;
				seg2.TofBirth = seg.TofBirth;
				if (grid.IsOutsideBox(seg2))
					seg = PeriodicBC(seg2,grid,frag,data);
				else
					return seg2;
			}
			
		}
	}
	
	else
	{
		unit_vec[0] = (seg.x[0]-seg.x[1]);
		unit_vec[1] = (seg.y[0]-seg.y[1]);
		unit_vec[2] = (seg.z[0]-seg.z[1]);
		length = vec_norm(unit_vec);
		unit_vec[0] /= length;
		unit_vec[1] /= length;
		unit_vec[2] /= length;
		for (n=0;n<6;++n)
		{
			if (intersectPlane(grid,seg,n,intersectpt))
			{
				//cout << endl << "Vessel " << seg.label << " Crossed Face " << n << endl;
				seg.x[0] = intersectpt[0];
				seg.y[0] = intersectpt[1];
				seg.z[0] = intersectpt[2];
				if (seg.length < 0)
				{
					seg.length = -sqrt((seg.x[1]-seg.x[0])*(seg.x[1]-seg.x[0]) +
					(seg.y[1]-seg.y[0])*(seg.y[1]-seg.y[0]) + (seg.z[1]-seg.z[0])*(seg.z[1]-seg.z[0]));
				}
				else
				{
					seg.length = sqrt((seg.x[1]-seg.x[0])*(seg.x[1]-seg.x[0]) +
					(seg.y[1]-seg.y[0])*(seg.y[1]-seg.y[0]) + (seg.z[1]-seg.z[0])*(seg.z[1]-seg.z[0]));
				}
				
				
				seg.tip[0] = 0;
                
                /*if (n > 3){
				    data.num_zdead = data.num_zdead+1;
				    return seg;}*/
				
				if (n > 3)
				    data.num_zdead = data.num_zdead+1;
				return seg;
				
				frag.push_front (seg);
				Segment seg2;
				data.num_vessel = data.num_vessel + 1;
				
				switch (n)
				{
				case 0:
					seg2.x[1] = seg.x[0];
					seg2.y[1] = oppface[n];
					seg2.z[1] = seg.z[0];
					break;
					
				case 1:
					seg2.x[1] = oppface[n];
					seg2.y[1] = seg.y[0];
					seg2.z[1] = seg.z[0];
					break;
				case 2:
					seg2.x[1] = seg.x[0];
					seg2.y[1] = oppface[n];
					seg2.z[1] = seg.z[0];
					break;
				case 3:
					seg2.x[1] = oppface[n];
					seg2.y[1] = seg.y[0];
					seg2.z[1] = seg.z[0];
					break;
				case 4:
					seg2.x[1] = seg.x[0];
					seg2.y[1] = seg.y[0];
					seg2.z[1] = oppface[n];
					break;
				case 5:
					seg2.x[1] = seg.x[0];
					seg2.y[1] = seg.y[0];
					seg2.z[1] = oppface[n];
					break;
				}

				if (seg.length < 0) 
				{
					rem_length = -length + sqrt((seg.x[1]-seg.x[0])*(seg.x[1]-seg.x[0])
					+ (seg.y[1]-seg.y[0])*(seg.y[1]-seg.y[0])
					+ (seg.z[1]-seg.z[0])*(seg.z[1]-seg.z[0]));
					seg2.x[0] = seg2.x[1] - rem_length*unit_vec[0];
					seg2.y[0] = seg2.y[1] - rem_length*unit_vec[1];
					seg2.z[0] = seg2.z[1] - rem_length*unit_vec[2];
					
				}
				else 
				{
					rem_length = length - sqrt((seg.x[1]-seg.x[0])*(seg.x[1]-seg.x[0])
					+ (seg.y[1]-seg.y[0])*(seg.y[1]-seg.y[0])
					+ (seg.z[1]-seg.z[0])*(seg.z[1]-seg.z[0]));
					seg2.x[0] = seg2.x[1] + rem_length*unit_vec[0];
					seg2.y[0] = seg2.y[1] + rem_length*unit_vec[1];
					seg2.z[0] = seg2.z[1] + rem_length*unit_vec[2];
				}
				
				seg2.tip[1] = 0;
				seg2.tip[0] = -1;
				seg2.label = seg.label;
				seg2.vessel = data.num_vessel;
				seg2.sprout = seg.sprout;
				seg2.length = rem_length;
//				seg2.phi1 = seg.phi1;
//				seg2.phi2 = seg.phi2;
				seg2.Recent_branch = seg.Recent_branch;
				seg2.TofBirth = seg.TofBirth;
				
				if (grid.IsOutsideBox(seg2))
					seg = PeriodicBC(seg2,grid,frag,data);
				else
					return seg2;
			}
		}
	}
	return seg;
}


