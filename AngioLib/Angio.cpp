///////////////////////////////////////////////////////////////////////
// Angio.cpp
///////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include <iostream>

#include "Angio.h"
#include "angio3d.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Angio::Angio(Input &input) : grid(input), data(input, grid)
{	

	
	first_adj = true;                                           // Declare and define first_adj as 'true', indicates t = 0
    half_cell_l = 0.5*data.dx;                                  // Half the length of one grid element in the x direction, use in determining variable time step
	p = 0.;
	q = 0.;

	kill_off = false;

    cout << endl << "Angiogenesis Growth Model:" << endl << endl;

	killed_segs_stream = fopen("out_dead_segs.ang","wt");

	cult.W[0] = input.iweight1;
	cult.W[3] = input.iweight4;

}

Angio::~Angio()
{
	fclose(killed_segs_stream);
}

///////////////////////////////////////////////////////////////////////
// Member Functions:
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// seedFrag
///////////////////////////////////////////////////////////////////////

void Angio::seedFrag()
{
    int i;
    
    for (i=0; i < data.NFRAGS; ++i)                             // Iterate from 1 to TOTAL NUMBER OF INITIAL FRAGMENTS (i)
	{
		Segment seg;                                            // Declare SEGMENT object 'seg'
		seg = cult.createInitFrag(data,grid,i,frag);            // CULTURE.createInitFrag(data, grid, vessel to be created, iterator, list container)
		                                                        // Creates the initial fragements for t = 0
		
		data.nsegs = data.nsegs+1;                              // Iterate the total segment counter +1
		seg.seg_num = data.nsegs;
		frag.push_back (seg);                                   // Append new segment at the bottom of the list 'frag'
	
		//seg_nodes.push_back(vect3(seg.x[0],seg.y[0],seg.z[0]));
		//seg_nodes.push_back(vect3(seg.x[1],seg.y[1],seg.z[1]));	
	}
	
	//update_elem_tags(grid);

	return;
}



///////////////////////////////////////////////////////////////////////
// initBranch
///////////////////////////////////////////////////////////////////////

void Angio::initBranch()
{
	list<Segment>::iterator it;
	
	for (it = frag.begin(); it != frag.end(); ++it)
	{
	     if (float(rand())/RAND_MAX < (data.a1+data.a2))        // Generate random number between 0 and 1
			it->init_branch = true;                             // If random number less than a1 + a2, then initial branching is true
		else
			it->init_branch = false;                            // If random number is not less than a1 + a2, then initial branching is false
	}
	
	
	return;
}



///////////////////////////////////////////////////////////////////////
// updateTime
///////////////////////////////////////////////////////////////////////

void Angio::updateTime()
{
    ///// Variable time step: Determine time step which causes vessel to grow half a grid cell /////
		
	//q = pow(E,-(data.t-data.x0)/data.b);                                                        // Define 'q' using sigmoidal growth rate function
	//p = -(half_cell_l + half_cell_l*q + data.a*q)/(half_cell_l + half_cell_l*q - data.a);       // Define 'p' based on half the length of the size of a grid element

	//if (p > 0)                                                  // For first time step, p < 0 so dt = dt
	//{
	//	if (first_adj)                                          // If culture just initially seeded...
	//	{
	//		data.dt = (data.t-data.x0+log(p)*data.b);             // Only applied for first time step that's greater than 0 (time step 2)
	//		first_adj = false;
	//	}
	//	else
	//		data.dt = (data.t-data.x0+log(p)*data.b);           // Applied to all preceding time steps
	//}
		
	data.dt = data.dt_a*pow(E,data.dt_b/(data.n + data.dt_c));    
    data.n = data.n + 1.0;

    if (data.dt > (data.maxt-data.t) && (data.maxt - data.t) > 0)       // If dt is bigger than the amount of time left...
		data.dt = data.maxt - data.t;                                       // then just set dt equal to the amount of time left (maxt-t)

	
    data.t = data.t+data.dt;                                    // Update time

    return;
}



///////////////////////////////////////////////////////////////////////
// updateTotalLength
///////////////////////////////////////////////////////////////////////

void Angio::updateTotalLength()
{
    data.total_length = 0.;                                  // Initialize total_length to 0
    list<Segment>::iterator it;
		
		
    for (it = frag.begin(); it != frag.end(); ++it)                 // Iterate through all segments stored in 'a_frag' (it)
    {
        data.total_length = data.total_length + fabs(it->length);        // Calculate total_length        
    }

    return;
}



///////////////////////////////////////////////////////////////////////
// Growth
///////////////////////////////////////////////////////////////////////

void Angio::Growth()
{
    int k;
    list<Segment>::iterator it;                                 // Declare iterator through list container 'it'
    int test = 0;

	//// Elongate the active vessel tips
    for (it = frag.begin(); it != frag.end(); ++it)             // Iterate through each vessel fragment (it)
	{
		//it->findphi();
		test = 1;

		for (k=0; k<=1; ++k)                                        // Iterate through each tip of the vessel fragment (k)
		{
	        if (it->tip[k] != 0)                                        // If tip is active (not 0) then create a new segment
			{
		        Segment seg;                                                // Declare SEGMENT object 'seg'
				seg = cult.createNewSeg(it,grid,data,k,frag);               // CULTURE.createNewSeg(list iterator, grid, data, tip, frag container)
				                                                            // Create new vessel segment on existing segment to 
				                                                            // 'elongate' vessel over time step 
				
				//cult.CheckForIntersection(seg,frag,data,it);                // CULTURE.CheckForIntersection(new segment, fragment container, data, list iterator)
				                                                            // Checks for intersection between new segment and any existing fragment
				                                                               // If true, '3D Insertion' notice will be printed
				
				
				
				frag.push_front (seg);                                      // Append new segment at the top of the list 'frag'
				
				/*if (seg.tip[k] == 1)
					seg_nodes.push_back(vect3(seg.x[1],seg.y[1],seg.z[1]));
				else if (seg.tip[k] == -1)
					seg_nodes.push_back(vect3(seg.x[0],seg.y[0],seg.z[0]));*/
			
			}
            
	    }
		    
		Branch(it);												// Determine which segments form new branches

    }

	Fuse();														// Determine which segments form new anastomoses
	
	//update_elem_tags(grid);
	   
	return;
}



///////////////////////////////////////////////////////////////////////
// updateLength
///////////////////////////////////////////////////////////////////////

void Angio::updateLength()
{
	double lc;                                                  // lc - Length calculation obtained from growth function g(t)
	//double nt;                                                  // nt - Number of active tips
    
	lc = data.a/(1.+pow(E,-(data.t-data.x0)/data.b));
	lc -= data.a/(1.+pow(E,-(data.t-data.dt-data.x0)/data.b));
	
	//nt = double(2*data.NFRAGS + data.num_branches - data.num_anastom);
	////nt = double(2*data.NFRAGS + data.num_branches - data.num_anastom - data.num_zdead);
	//
	//if (nt <= 0)
	//	nt = double(2*data.NFRAGS + data.num_branches - data.num_anastom - data.num_zdead);
	//	//nt = 4*data.NFRAGS;
	
	data.vess_length = lc*data.length_adjust;
	//data.vess_length = lc*grid.den_scale*data.length_adjust;

	return;
}


///////////////////////////////////////////////////////////////////////
// Branch
///////////////////////////////////////////////////////////////////////

void Angio::Branch(list<Segment>::iterator it)
{
    int k;
    //list<Segment>::iterator it;                                         // Declare iterator through list container 'it'
    
    // Generate a random number between 0 and 1. If that number is less than the branching probability at
    // the current time, or initial branching is activated (data.ini_branch = true), then branching occurs
	
    double den_scale = 1.0;

	double xpt = (it->x[1] + it->x[0])/2;
	double ypt = (it->y[1] + it->y[0])/2;
	double zpt = (it->z[1] + it->z[0])/2;
	
	den_scale = cult.findDenScale(xpt, ypt, zpt, grid);
	
	//if ( float(rand())/RAND_MAX < grid.den_scale*data.dt*data.branch_chance/data.t || (it->init_branch == true) )
	if ( float(rand())/RAND_MAX < den_scale*data.dt*data.branch_chance/data.t || (it->init_branch == true) )
    {
	    if ((it->BCdead == 0) && (it->anast == 0))                  // Segments that have encountered a boundary condition or formed an anastomoses may not
		{                                                           // form a branch.
	        
			Segment seg;                                            // Declare SEGMENT object 'seg'
			data.num_branches = data.num_branches + 1;              // Iterate the total number of branches +1
			data.branch = true;                                     // Branching flag set to 'true.' This tells the program that the
			                                                        // new vessel segment being created is arising from a branch      
    		it->init_branch = false;
    				
			if (float(rand())/RAND_MAX < 0.5)                       // Randomly determine which node of the parent segment forms the branch
			    k = 0;
		    else
			    k = 1;
    							
			it->tip[k] = sign(0.5f - float(rand())/RAND_MAX);        // Randomly assign the branch to grow as +1 or -1
    				
			seg = cult.createNewSeg(it,grid,data,k,frag);           // CULTURE.createNewSeg(list iterator, grid, data, tip, frag container)
			                                                            // Create new vessel segment on existing segment to 
			                                                            // create the branching 

			//if (seg.tip[k] == 1)
			//		seg_nodes.push_back(vect3(seg.x[1],seg.y[1],seg.z[1]));
			//	else if (seg.tip[k] == -1)
			//		seg_nodes.push_back(vect3(seg.x[0],seg.y[0],seg.z[0]));
			
            data.num_vessel = data.num_vessel + 1;
		    seg.vessel = data.num_vessel;
    				                    
			it->Recent_branch = 1;
			seg.Recent_branch = 1;                                  // Indicate that this segment just branched, prevents it from 
				                                                        // branching again too soon
    		
			//++data.nsegs;                                           // Iterate the total segment counter +1
			//seg.seg_num = data.nsegs;

			frag.push_front (seg);                                  // Append new segment at the top of the list 'frag'
			
			data.branch = false;                                    // Turn off branching flag once branching algorithm is complete
        }
	}                                                               // End Branching
	

	return;
}



///////////////////////////////////////////////////////////////////////
// Fuse
///////////////////////////////////////////////////////////////////////

void Angio::Fuse()
{
    list<Segment>::iterator it;
           
    for (it = frag.begin(); it != frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		check4anast(it, 0);
        check4anast(it, 1);
	} 
			
    //update_elem_tags(grid);

	return;
}



///////////////////////////////////////////////////////////////////////
// check4anast
///////////////////////////////////////////////////////////////////////

void Angio::check4anast(list<Segment>::iterator it, int k)
{
    if (it->tip[k] == 0)
        return;
	
	if (it->anast == 1)
		return;

	int kk = 0;
    double dist0 = 0.;
    double dist1 = 0.;
    list<Segment>::iterator it2;
	
    for (it2 = frag.begin(); it2 != frag.end(); ++it2)          // Iterate through all segments in frag list container again (it2)
	{                                                           
	    dist0 = (it->x[k] - it2->x[0])*(it->x[k] - it2->x[0]) + (it->y[k] - it2->y[0])*(it->y[k] - it2->y[0]) + (it->z[k] - it2->z[0])*(it->z[k] - it2->z[0]);
        dist1 = (it->x[k] - it2->x[1])*(it->x[k] - it2->x[1]) + (it->y[k] - it2->y[1])*(it->y[k] - it2->y[1]) + (it->z[k] - it2->z[1])*(it->z[k] - it2->z[1]);

        anastomose(dist0, dist1, k, it, it2);
    } 
 	
	//int elem_num = grid.findelem(it->x[k], it->y[k], it->z[k]);
	//
	//scan_4_segs(it, k, elem_num, data, grid);							// Check current element for nearby segs
	//search_nearby_neighbors(it, k, elem_num, data, grid);

	return;           
}



///////////////////////////////////////////////////////////////////////
// search_nearby_neighbors
///////////////////////////////////////////////////////////////////////

//void Angio::search_nearby_neighbors(list<Segment>::iterator it, int k, int elem_num)
//{
//	// Check all 26 neighbors of elem_num
//	int fronty = grid.elem_find_neighbor(elem_num,1);					// Find fronty
//	int righty = grid.elem_find_neighbor(elem_num,2);					// Find righty
//	int backy = grid.elem_find_neighbor(elem_num,3);					// Find backy
//	int lefty = grid.elem_find_neighbor(elem_num,4);					// Find lefty
//	int toppy = grid.elem_find_neighbor(elem_num,5);					// Find toppy
//	int bottomy = grid.elem_find_neighbor(elem_num,6);					// Find bottomy
//		
//	scan_4_segs(it, k, fronty);								// Check fronty for nearby segs (front, center, middle)
//	scan_4_segs(it, k, righty);								// Check righty for nearby segs (core, right, middle)
//	scan_4_segs(it, k, backy);								// Check backy for nearby segs (back, center, middle)
//	scan_4_segs(it, k, lefty);								// Check lefty for nearby segs (core, left, middle)
//	scan_4_segs(it, k, toppy);								// Check toppy for nearby segs (core, center, top)
//	scan_4_segs(it, k, bottomy);							// Check bottomy for nearby segs (core, center, bottom)
//
//	//scan_4_segs(it, k, grid.elem_find_neighbor(righty,5), data, grid);	// Check righty's top (core, right, top)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(righty,6), data, grid);	// Check righty's bottom (core, right, bottom)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(lefty,5), data, grid);	// Check lefty's top (core, left, top)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(lefty,6), data, grid);	// Check lefty's bottom (core, left, bottom)
//	//
//	//scan_4_segs(it, k, grid.elem_find_neighbor(fronty,5), data, grid);	// Check fronty's top (front, center, top)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(grid.elem_find_neighbor(fronty,5),2), data, grid);	// Check fronty's toppy's right (front, right, top)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(grid.elem_find_neighbor(fronty,5),4), data, grid);	// Check fronty's toppy's left (front, left, top)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(fronty,2), data, grid);	// Check fronty's right (front, right, middle)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(fronty,4), data, grid);	// Check fronty's left (front, left, middle)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(fronty,6), data, grid);	// Check fronty's bottom (front, center, bottom)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(grid.elem_find_neighbor(fronty,6),2), data, grid);	// Check fronty's bottomy's right (front, right, bottom)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(grid.elem_find_neighbor(fronty,6),4), data, grid);	// Check fronty's bottomy's left (front, left, bottom)
//
//	//scan_4_segs(it, k, grid.elem_find_neighbor(backy,5), data, grid);	// Check backy's top (back, center, top)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(grid.elem_find_neighbor(backy,5),2), data, grid);	// Check backy's toppy's right (back, right, top)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(grid.elem_find_neighbor(backy,5),4), data, grid);	// Check backy's toppy's left (back, left, top)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(backy,2), data, grid);	// Check backy's right (back, right, middle)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(backy,4), data, grid);	// Check backy's left (back, left, middle)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(backy,6), data, grid);	// Check backy's bottom (back, center, bottom)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(grid.elem_find_neighbor(backy,6),2), data, grid);	// Check backy's bottomy's right (back, right, bottom)
//	//scan_4_segs(it, k, grid.elem_find_neighbor(grid.elem_find_neighbor(backy,6),4), data, grid);	// Check backy's bottomy's left (back, left, bottom)
//
//	return;
//}



///////////////////////////////////////////////////////////////////////
// scan_4_segs
///////////////////////////////////////////////////////////////////////

//void Angio::scan_4_segs(list<Segment>::iterator it, int k, int elem_num)
//{
//	if (elem_num == -1)
//		return;
//
//	list<list<Segment>::iterator >::iterator seg_it;
//	list<Segment>::iterator it2;
//	
//	double dist0 = 0.;
//    double dist1 = 0.;
//
//	for (seg_it = grid.ebin[elem_num].resident_segs.begin(); seg_it != grid.ebin[elem_num].resident_segs.end(); ++seg_it)
//	{
//		it2 = (*seg_it);
//		dist0 = sqrt((it->x[k] - it2->x[0])*(it->x[k] - it2->x[0]) + (it->y[k] - it2->y[0])*(it->y[k] - it2->y[0]) + (it->z[k] - it2->z[0])*(it->z[k] - it2->z[0]));
//        dist1 = sqrt((it->x[k] - it2->x[1])*(it->x[k] - it2->x[1]) + (it->y[k] - it2->y[1])*(it->y[k] - it2->y[1]) + (it->z[k] - it2->z[1])*(it->z[k] - it2->z[1]));
//		anastomose(dist0, dist1, k, it, it2);
//	}
//	
//	return;
//}



///////////////////////////////////////////////////////////////////////
// anastomose
///////////////////////////////////////////////////////////////////////

void Angio::anastomose(double dist0, double dist1, int k, list<Segment>::iterator it, list<Segment>::iterator it2)
{
    //double anast_dist = data.vess_length;
	double anast_dist = data.anast_dist;

	if ((it->anast == 1) || (it2->anast == 1))
        return;
    
    if (it->label == it2->label)
        return;

    int kk = 9;
    
    if (dist0 <= anast_dist)
        kk = 0;
    else if (dist1 <= anast_dist)
        kk = 1;
  
    if (kk == 9)
        return;
                                                
	//if (it2->tip[kk] == 0)                                      // Tip-to-tip anastomose only
    //  return;
     
    Segment seg;                                                // Declare SEGMENT object 'seg'
								
	seg = cult.connectSegment(it,it2,k,kk,grid,data,frag);      // CULTURE.connectSegment(segment 1, segment 2, tip 1, tip 2, grid, data, frag list container)
					                                            // This function will create a segment between to two segments to complete the anastomosis
	frag.push_front (seg);                                      // Append new segment at the top of the list 'frag'                  
	it->tip[k]=0;                                               // Deactivate tip of segment 1 after anastomosis
	it2->tip[kk]=0;                                             // Deactivate tip of segment 2 after anastomosis (tip-tip anastomosis only)

    
	//update_elem_tags(grid);
	
	return;
}



///////////////////////////////////////////////////////////////////////
// createFragList
///////////////////////////////////////////////////////////////////////

list<Segment> Angio::returnFragList()
{
    return frag;
}


///////////////////////////////////////////////////////////////////////
// removeErrors
///////////////////////////////////////////////////////////////////////

void Angio::removeErrors()
{
	if (kill_off == false){
		double length_limit = data.d;
    
		list<Segment>::iterator it;
    
		for (it = frag.begin(); it != frag.end(); ++it){
			if (fabs(it->length) >= 2.0*length_limit){
				it->death_label = -7;
				fprintf(killed_segs_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i\n",it->TofBirth,it->x[0],it->y[0],it->z[0],it->x[1],it->y[1],it->z[1],it->length,it->death_label);
				it = frag.erase(it);}
			if (it == frag.end())
				return;
		}
                
		return;
	}
}



///////////////////////////////////////////////////////////////////////
// displace_vessels
///////////////////////////////////////////////////////////////////////

void Angio::displace_vessels()
{
    int k = 0; int elem_num = 0;
    list<Segment>::iterator it;
    double xpt = 0.; double ypt = 0.; double zpt = 0.;
    double xix = 0.; double xiy = 0.; double xiz = 0.;
    double shapeF[8] = {0.};
    vect3 disp;
    Elem elem;
    
    for (it = frag.begin(); it != frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		for (k=0; k<=1;++k)                                         // Iterate through both segment tips (k)
		{
		    xpt = it->x[k];
		    ypt = it->y[k];
		    zpt = it->z[k];
		    
			int BC_face = 0;
			elem_num = it->tip_elem[k];
		    
			if (elem_num != -1)
			{
				elem = grid.ebin[elem_num];
		    
				grid.natcoord(xix, xiy, xiz, xpt, ypt, zpt, elem_num);
		    
				/*if (xix >= 1.0)
					xix = 0.99;
				if (xix <= -1.0)
					xix = -0.99;
				if (xiy >= 1.0)
					xiy = 0.99;
				if (xiy <= -1.0)
					xiy = -0.99;
				if (xiz >= 1.0)
					xiz = 0.99;
				if (xiz <= -1.0)
					xiz = -0.99;		*/		
				
				grid.shapefunctions(shapeF, xix, xiy, xiz);
		    
				disp.x = shapeF[0]*(*elem.n1).u.x + shapeF[1]*(*elem.n2).u.x + shapeF[2]*(*elem.n3).u.x + shapeF[3]*(*elem.n4).u.x + shapeF[4]*(*elem.n5).u.x + shapeF[5]*(*elem.n6).u.x + shapeF[6]*(*elem.n7).u.x + shapeF[7]*(*elem.n8).u.x;
				disp.y = shapeF[0]*(*elem.n1).u.y + shapeF[1]*(*elem.n2).u.y + shapeF[2]*(*elem.n3).u.y + shapeF[3]*(*elem.n4).u.y + shapeF[4]*(*elem.n5).u.y + shapeF[5]*(*elem.n6).u.y + shapeF[6]*(*elem.n7).u.y + shapeF[7]*(*elem.n8).u.y;
				disp.z = shapeF[0]*(*elem.n1).u.z + shapeF[1]*(*elem.n2).u.z + shapeF[2]*(*elem.n3).u.z + shapeF[3]*(*elem.n4).u.z + shapeF[4]*(*elem.n5).u.z + shapeF[5]*(*elem.n6).u.z + shapeF[6]*(*elem.n7).u.z + shapeF[7]*(*elem.n8).u.z;
		    		    
				it->x[k] = it->x[k] + disp.x;
				it->y[k] = it->y[k] + disp.y;
				it->z[k] = it->z[k] + disp.z;
			}
			else
			{
				//it->x[0] = 0.;
				//it->x[1] = 0.;
				//it->y[0] = 0.;
				//it->y[1] = 0.;
				//it->z[0] = 0.;
				//it->z[1] = 0.;
				//it->tip[0] = 0;
				//it->tip[1] = 0;
				////it->length = 0.;
				it->mark_of_death = true;
				it->death_label = 9;
			}			
	    }
        
		//it->findlength();
		
		//if ((it->tip[0] == -1) || (it->tip[1] == 1))
			//it->findphi();
    }
    
    return;
}   
    


///////////////////////////////////////////////////////////////////////
// find_active_tips()
///////////////////////////////////////////////////////////////////////

void Angio::find_active_tips()
{
	list<Segment>::iterator it;
           
	active_tips.clear();
	
	for (it = frag.begin(); it != frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		if (it->tip[0] != 0){
			active_tips.push_back(it);}
		else if (it->tip[1] != 0){
			active_tips.push_back(it);}
	} 
			
    return;
}



///////////////////////////////////////////////////////////////////////
// find_newborns()
///////////////////////////////////////////////////////////////////////

void Angio::find_newborns()
{
	list<Segment>::iterator it;
           
	newborns.clear();
	
	for (it = frag.begin(); it != frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	{
		if (it->TofBirth == data.t)
			newborns.push_back(it);
	} 
			
    return;
}



///////////////////////////////////////////////////////////////////////
// update_elem_tags
///////////////////////////////////////////////////////////////////////

void Angio::update_elem_tags()
{
	//list<Segment>::iterator it;
	//int elem0 = 0; int elem1 = 0;
	//
	//for (it = frag.begin(); it != frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
	//{
	//	if (it->elem_tagged == false)
	//	{
	//		elem0 = grid.findelem(it->x[0],it->y[0],it->z[0]);
	//		elem1 = grid.findelem(it->x[1],it->y[1],it->z[1]);

	//		if (elem0 == elem1){
	//			grid.ebin[elem0].resident_segs.push_back(it);}
	//	    else{
	//			grid.ebin[elem0].resident_segs.push_back(it);
	//			grid.ebin[elem1].resident_segs.push_back(it);}

	//		it->elem_tagged = true;
	//	}
	//}

	return;
}



///////////////////////////////////////////////////////////////////////
// kill_dead_segs
///////////////////////////////////////////////////////////////////////

void Angio::kill_dead_segs()
{  
	if (kill_off == false){
		list<Segment>::iterator it;
    
		for (it = frag.begin(); it != frag.end(); ++it){
			if (it->mark_of_death == true){
				fprintf(killed_segs_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i\n",it->TofBirth,it->x[0],it->y[0],it->z[0],it->x[1],it->y[1],it->z[1],it->length,it->death_label);
				it = frag.erase(it);}
			if (it == frag.end())
				return;
		}
	}
                
    return;
}



///////////////////////////////////////////////////////////////////////
// initialize_seg_positions
///////////////////////////////////////////////////////////////////////

void Angio::initialize_seg_positions()
{
	list<Segment>::iterator frag_it;

	for (frag_it = frag.begin(); frag_it != frag.end(); ++frag_it)
	{
		frag_it->x0[0] = frag_it->x[0];
		frag_it->x0[1] = frag_it->x[1];
		frag_it->y0[0] = frag_it->y[0];
		frag_it->y0[1] = frag_it->y[1];
		frag_it->z0[0] = frag_it->z[0];
		frag_it->z0[1] = frag_it->z[1];
	}

	return;
}



///////////////////////////////////////////////////////////////////////
// enforce_global_BC
///////////////////////////////////////////////////////////////////////

void Angio::enforce_global_BC()
{
	list<Segment>::iterator it;
    
	double eps = 0.001;
    
	double xmin = grid.xrange[0]*(1 - eps);
	double xmax = grid.xrange[1]*(1 + eps);
	if (xmin == 0.0)
		xmin = -1.0*eps*xmax;
		
	double ymin = grid.yrange[0]*(1 - eps);
	double ymax = grid.yrange[1]*(1 + eps);
	if (ymin == 0.0)
		ymin = -1.0*eps*ymax;
	
	double zmin = grid.zrange[0]*(1 - eps);
	double zmax = grid.zrange[1]*(1 + eps);
	if (zmin == 0.0)
		zmin = -1.0*eps*zmax;
	
	for (it = frag.begin(); it != frag.end(); ++it){
		for (int k=0; k<=1;++k){
			if (it->x[k] < xmin){
				it->x[k] = grid.xrange[0];
				it->tip[k] = 0;}
			
			if (it->x[k] > xmax){
				it->x[k] = grid.xrange[1];
				it->tip[k] = 0;}

			if (it->y[k] < ymin){
				it->y[k] = grid.yrange[0];
				it->tip[k] = 0;}

			if (it->y[k] > ymax){
				it->y[k] = grid.yrange[1];
				it->tip[k] = 0;}

			if (it->z[k] < zmin){
				it->z[k] = grid.zrange[0];
				it->tip[k] = 0;}

			if (it->z[k] > zmax){
				it->z[k] = grid.zrange[1];
				it->tip[k] = 0;}
		}
	}
	
	return;
}
