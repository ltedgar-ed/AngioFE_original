///////////////////////////////////////////////////////////////////////
// Segment.cpp
///////////////////////////////////////////////////////////////////////



// Include:
#include "stdafx.h"

#include "Segment.h"
#include "vect3.h"
#include "angio3d.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Segment::Segment()                                              // Constructor for SEGMENT object   
{
    x[0]= 0;                                                    // x, y, z - Initialize the coordinates of the tips to 0                                                                       
	x[1]= 0;
	y[0]= 0;
	y[1]= 0;
	z[0]= 0;
	z[1]= 0;
	
	x0[0]= 0;                                                    // x, y, z - Initialize the coordinates of the tips to 0                                                                       
	x0[1]= 0;
	y0[0]= 0;
	y0[1]= 0;
	z0[0]= 0;
	z0[1]= 0;

	tip[0]= 0;                                                  // tip - Initialize the activity of the tips to 0
	tip[1]= 0;
	length = 0;                                                 // length - Initialize the length of the segment to 0
	
//	phi1 = 0;                                                   // phi1, phi2 - Initialize the orientation angles of the segment to 0
//	phi2 = 0;
	
	label = 0;                                                  // label - Initialize label to 0
    vessel = 0;                                                 // vessel - Initialize vessel to 0
    
    BCdead = 0;                                                 // BCdead - Set boundary condition indicator to 'false'
	TofBirth = 0;                                               // TofBirth - Initialize time of birth to 0
	Recent_branch = 0;                                          // Recent_branch - Initialize branching indicator to 0
	init_branch = false;                                        // init)branch - Set initial branching flag to 'false'
	sprout = 0;                                                 // sprout - Initialize sprout to 0                                                                                                                                            

    anast = 0;

	elem_tagged = false;

	bdyf_id[0] = -1;
	bdyf_id[1] = -1;
	
	mark_of_death = false;
	death_label = 0;
	
	tip_BC[0] = 0;
	tip_BC[1] = 0;

	tip_elem[0] = -1;
	tip_elem[1] = -1;

	seg_num = 0;

	seg_conn[0][0] = 0;
	seg_conn[0][1] = 0;
	seg_conn[1][0] = 0;
	seg_conn[1][1] = 0;

	//line_num = 0;
}

Segment::~Segment()                                             // Destructor for SEGMENT object
{

}

void Segment::findlength()
{	
	double new_length = sqrt(pow((x[1] - x[0]),2) + pow((y[1] - y[0]),2) + pow((z[1] - z[0]),2));
    
    if (length < 0.0)
        length = -new_length;
	else
		length = new_length;
        
    uvect.x = x[1] - x[0];
	uvect.y = y[1] - y[0];
	uvect.z = z[1] - z[0];

	if (length != 0)
		uvect = uvect/uvect.norm();
	
	return;
}


void Segment::findunit()
{	
	uvect.x = x[1] - x[0];
	uvect.y = y[1] - y[0];
	uvect.z = z[1] - z[0];

	if (uvect.norm() != 0.0)
		uvect = uvect/uvect.norm();
	
	return;
}

//void Segment::findphi()
//{
//    vect3 vvect, xvvect, zvvect;
//           
//	if (length == 0.0)
//		return;
//	
//	vvect.x = x[1] - x[0];
//	vvect.y = y[1] - y[0];    
//	vvect.z = z[1] - z[0];
//	vvect = vvect/vvect.norm();
//
//    xvvect = vect3(vvect.x, vvect.y, 0);
//    
////	double phi10 = phi1;
////	double phi20 = phi2;
//
//	// Calculate phi2
//	phi2 = asin(vvect.z);
//
//	double alt_phi2 = 0.;
//
//	if (vvect.x < 0)
//		if (phi2 > 0)
//			alt_phi2 = pi - phi2;
//		else if (phi2 < 0)
//			alt_phi2 = -pi - phi2;
//			
//	if (fabs(phi20 - alt_phi2) < fabs(phi20 - phi2))
//		phi2 = alt_phi2;
//
//	// Calculate phi1
//	double phi1_dot = 0.;
//	phi1_dot = xvvect*vect3(1,0,0);
//
//	if (cos(phi2) != 0.0)
//		phi1 = acos(vvect.x/cos(phi2)); 
//	else				
//		phi1 = acos(phi1_dot);
//						
//	if (vvect.x/cos(phi2) > 1)
//				phi1 = 0;
//			else if (vvect.x/cos(phi2) < -1)
//				if (phi10 > 0)
//					phi1 = pi;
//				else if (phi10 < 0)
//					phi1 = -pi;
//	
//	if (vvect.y < 0)
//		phi1 = -phi1;
//
//	if (vvect.norm() == 0.){
//		phi1 = phi10;
//		phi2 = phi20;}
//						
//	if (vvect.y == 0.)
//		phi1 = phi10;
//						
//	if (vvect.z == 0.)
//		phi2 = 0.;
//				
//    return;
//}
