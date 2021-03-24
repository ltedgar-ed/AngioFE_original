///////////////////////////////////////////////////////////////////////
// Filein.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"

#include <iostream>
#include <fstream>

#include "FEAngInput.h"
#include "Filein.h"
#include "angio3d.h"
#include "Elem.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Filein::Filein()
{
    
}

Filein::~Filein()
{

}

bool Filein::Input(const char* szfile, FEAngInput &input, FEMesh &mesh)
{
	FILE* fp = fopen(szfile, "rt");
	if (fp == 0) return false;

    const int buff_max = 100;										// Maximum size of the buffer
    char buffer[buff_max];											// Create the buffer
    char first_char;												// First character of the current line in the buffer
    char str_type = 'n';											// Set default string type to 'no line'
    
    while (!feof(fp))											// Read until you reach the end of the file...
    {
        str_type = 'n';													// Reset the default string type to 'no line'
		fgets(buffer, 255, fp);
                
        first_char = buffer[0];											// Determine the first character of the line in the buffer
       
        if (first_char == 37)											// If first character is a '%' symbol...
             str_type = 'c';												// ... then set string type to 'Comment'
             else if (first_char == 62)										// If first character is a '>' symbol...
             str_type = 'p';												// ... then set string type to 'Param'
                      
        if (str_type == 'p')
            read_param(input, buffer);									// If the line describes a parameter, then read in that parameter
    }

	fclose(fp);

	read_FEmesh(input, mesh);										// Use the mesh from FEBio to create the grid for angio3d

	return true;
}

//////////////////////////////////////////////////////////////////////
// Member Functions
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// read_param
///////////////////////////////////////////////////////////////////////

void Filein::read_param(FEAngInput &input, char (&buffer)[100])
{
    char pname[20] = {0};											// Create string for the parameter name
    
    sscanf(buffer,"%*s %s %*f",&pname);								// Scan the buffer for the name of the parameter
    
    set_param(input,buffer,pname);									// Set the parameter based on the name
    
    return;
}



///////////////////////////////////////////////////////////////////////
// set_param
///////////////////////////////////////////////////////////////////////
    
void Filein::set_param(FEAngInput &input, char (&buffer)[100], char (&pname)[20])
{
    //// Parameters for angio3d (In paratheneses, string identifying the parameter and an example value):
	// Branching Probability (brnch_ch 0.1)
	if (!strcmp(pname,"brnch_ch")){
        sscanf(buffer,"%*s %*s %f",&input.ibrnch_ch);
        return;}

    //  Matrix conditions (matx_cnd 0) random
	if (!strcmp(pname,"matx_cnd")){
        sscanf(buffer,"%*s %*s %i",&input.imatx_cnd);
        return;}
    
	// Initial matrix density (matx_den 3.0) mg/mL
    if (!strcmp(pname,"matx_den")){
        sscanf(buffer,"%*s %*s %lf",&input.imatx_den);
        return;}  
        
	// Number of initial fragments (nfrag 70, based on 30K frags/mL)
    if (!strcmp(pname,"nfrag")){
        sscanf(buffer,"%*s %*s %i",&input.infrag);
        return;}        
    
	// End of culture period (max_time 6.0) days
    if (!strcmp(pname,"max_time")){
        sscanf(buffer,"%*s %*s %lf",&input.imax_time);
        return;}  
    
	// Initial time step (dt 0.25) days
    if (!strcmp(pname,"dt")){
        sscanf(buffer,"%*s %*s %lf",&input.idt);
        return;}  
    
	// Anastomosis distance (anst_dst 25.0) um 
    if (!strcmp(pname,"anst_dst")){
        sscanf(buffer,"%*s %*s %lf",&input.ianst_dst);
        return;}  
        
    // Segment length adjustment scale (lngth_adj 1.0) 
	if (!strcmp(pname,"lngth_adj")){
        sscanf(buffer,"%*s %*s %lf",&input.ilngth_adj);
        return;}
    
	// Number of nodes in x-direction for autogrid (xnodes 7)
    if (!strcmp(pname,"xnodes")){
        sscanf(buffer,"%*s %*s %i",&input.ixnodes);
        return;}  
    
	// Number of nodes in y-direction for autogrid (ynodes 7)
    if (!strcmp(pname,"ynodes")){
        sscanf(buffer,"%*s %*s %i",&input.iynodes);
        return;} 
    
	// Number of nodes in z-direction for autogrid (znodes 3)
    if (!strcmp(pname,"znodes")){
        sscanf(buffer,"%*s %*s %i",&input.iznodes);
        return;} 
    
	// Alternative way of specifying the number of nodes in each direction.  Reads three int in succession indicating the number of nodes in the x-, y-, and z-
	// direction, respectively (num_nodes 7 7 3)
    if (!strcmp(pname,"num_nodes")){
        sscanf(buffer,"%*s %*s %i %i %i",&input.ixnodes,&input.iynodes,&input.iznodes);
        return;}
    
	// Read in the minimum and maximum dimension of the domain in the x-direction (xrange 0 300) um
    if (!strcmp(pname,"xrange")){
        sscanf(buffer,"%*s %*s %lf %lf",&input.ixmin,&input.ixmax);
        return;}
    
	// Read in the minimum and maximum dimension of the domain in the y-direction (yrange 0 300) um
    if (!strcmp(pname,"yrange")){
        sscanf(buffer,"%*s %*s %lf %lf",&input.iymin,&input.iymax);
        return;}
    
	// Read in the minimum and maximum dimension of the domain in the z-direction (zrange 0 300) um
    if (!strcmp(pname,"zrange")){
        sscanf(buffer,"%*s %*s %lf %lf",&input.izmin,&input.izmax);
        return;}
    
	// Specify which type of boundary condition to enforces at the faces of the domain normal to the x-axis (x_bc w) flat wall BC
    if (!strcmp(pname,"x_bc")){
        sscanf(buffer,"%*s %*s %c",&input.ixbctype);
        return;} 
    
    // Specify which type of boundary condition to enforces at the faces of the domain normal to the y-axis (y_bc w) flat wall BC
	if (!strcmp(pname,"y_bc")){
        sscanf(buffer,"%*s %*s %c",&input.iybctype);
        return;}
        
    // Specify which type of boundary condition to enforces at the faces of the domain normal to the z-axis (z_bc w) flat wall BC
	if (!strcmp(pname,"z_bc")){
        sscanf(buffer,"%*s %*s %c",&input.izbctype);
        return;}   

	// Specify whether or not to create the grid using the autogrid constructor or read in the grid from input files (igrid n) autogrid
	if (!strcmp(pname,"igrid")){
        sscanf(buffer,"%*s %*s %c",&input.igrid_in);
        return;} 

	// Specify if a composite material model is being using 
	if (!strcmp(pname,"comp_mat")){
        sscanf(buffer,"%*s %*s %i",&input.icomp_mat);
        return;}
    
	// Specify if a composite material model is being using 
	if (!strcmp(pname,"sproutf")){
        sscanf(buffer,"%*s %*s %lf",&input.isproutf);
        return;}

	// Specify if a composite material model is being using 
	if (!strcmp(pname,"tip_range")){
        sscanf(buffer,"%*s %*s %lf",&input.itip_range);
        return;}

	// Specify if a composite material model is being using 
	if (!strcmp(pname,"sub_cyc")){
        sscanf(buffer,"%*s %*s %i",&input.isub_cyc);
        return;}

	// Specify which type of boundary condition to enforces at the faces of the domain normal to the y-axis (y_bc w) flat wall BC
	if (!strcmp(pname,"gelbc")){
        sscanf(buffer,"%*s %*s %c",&input.igelbc);
        return;}

	// Specify which type of boundary condition to enforces at the faces of the domain normal to the y-axis (y_bc w) flat wall BC
	if (!strcmp(pname,"stiff_fact")){
        sscanf(buffer,"%*s %*s %lf",&input.istifffact);
        return;}

	// Specify which type of boundary condition to enforces at the faces of the domain normal to the y-axis (y_bc w) flat wall BC
	if (!strcmp(pname,"sprout_verify")){
        sscanf(buffer,"%*s %*s %i",&input.isprout_verify);
        return;}

	// Specify which type of boundary condition to enforces at the faces of the domain normal to the y-axis (y_bc w) flat wall BC
	if (!strcmp(pname,"sprout_sphere")){
        sscanf(buffer,"%*s %*s %i",&input.isp_sphere);
        return;}

		// Specify which type of boundary condition at the front edge of the gel
	if (!strcmp(pname,"front_bc")){
        sscanf(buffer,"%*s %*s %c",&input.ifrontbc);
        return;}
	
	// Specify which type of boundary condition at the right edge of the gel
	if (!strcmp(pname,"right_bc")){
        sscanf(buffer,"%*s %*s %c",&input.irightbc);
        return;}
	
	// Specify which type of boundary condition at the back edge of the gel
	if (!strcmp(pname,"back_bc")){
        sscanf(buffer,"%*s %*s %c",&input.ibackbc);
        return;}
	
	// Specify which type of boundary condition at the left edge of the gel
	if (!strcmp(pname,"left_bc")){
        sscanf(buffer,"%*s %*s %c",&input.ileftbc);
        return;}
	
	// Specify which type of boundary condition at the bottom edge of the gel
	if (!strcmp(pname,"bottom_bc")){
        sscanf(buffer,"%*s %*s %c",&input.ibottombc);
        return;}
	
	// Specify which type of boundary condition at the top edge of the gel
	if (!strcmp(pname,"top_bc")){
        sscanf(buffer,"%*s %*s %c",&input.itopbc);
        return;}

	// Specify the location of the x symmetry panel
	if (!strcmp(pname,"Sx")){
        sscanf(buffer,"%*s %*s %lf",&input.iSx);
        return;}

	// Specify the location of the y symmetry panel
	if (!strcmp(pname,"Sy")){
        sscanf(buffer,"%*s %*s %lf",&input.iSy);
        return;}

	// Specify the location of the z symmetry panel
	if (!strcmp(pname,"Sz")){
        sscanf(buffer,"%*s %*s %lf",&input.iSz);
        return;}

	// Specify the seed number for the random generator
	if (!strcmp(pname,"rseed")){
        sscanf(buffer,"%*s %*s %i",&input.irseed);
        return;}

	// Specify if branching is on or off
	if (!strcmp(pname,"branch")){
        sscanf(buffer,"%*s %*s %i",&input.ibranch);
        return;}
	
	// Specify if anastomosis is on or off
	if (!strcmp(pname,"anast")){
        sscanf(buffer,"%*s %*s %i",&input.ianast);
        return;}

	// Specify the factor for the directional sprout force
	if (!strcmp(pname,"spfact")){
        sscanf(buffer,"%*s %*s %lf",&input.ispfactor);
        return;}

	// Specify to 'flatten' fibers in z
	if (!strcmp(pname,"zfibflat")){
        sscanf(buffer,"%*s %*s %lf",&input.izfibflat);
        return;}

	// Read in weights for determing the direction of growth
    	if (!strcmp(pname,"gweights")){
        sscanf(buffer,"%*s %*s %lf %lf",&input.iweight1,&input.iweight4);
        return;}

	// Specify whether or not to create the grid using the autogrid constructor or read in the grid from input files (igrid n) autogrid
	if (!strcmp(pname,"circ")){
        sscanf(buffer,"%*s %*s %i",&input.icirc);
        return;} 

    return;
}



///////////////////////////////////////////////////////////////////////
// read_FEmesh
///////////////////////////////////////////////////////////////////////

void Filein::read_FEmesh(FEAngInput &input, FEMesh &mesh)
{
	bool new_prob = true;
		
	//// Read in the nodes from the FEBio mesh
	input.iNn = mesh.Nodes();										// Read in the total number of nodes from the FEBio mesh
	vect3 fiber;
	double radius = 0.;

	if (input.icirc == 1)
		radius = input.ixmax;

	for (int i = 0; i < input.iNn; ++i)								// Iterate through all nodes...
	{
		Node node;														// Create a new node
		node.id = i;													// Give the node it's ID 
		node.x = mesh.Node(i).m_r0.x;									// Set the current x position as FEBio initial x position
		node.y = mesh.Node(i).m_r0.y;									// Set the current y position as FEBio initial y position
		node.z = mesh.Node(i).m_r0.z;									// Set the current z position as FEBio initial z position 
		
		node.x0 = node.x;												// Set initial x position to current x position
		node.y0 = node.y;												// Set initial y position to current y position 						
		node.z0 = node.z;												// Set initial z position to current z position

		node.ecm_den = input.imatx_den;									// Set the density of the ECM at the node 
		node.ecm_den0 = node.ecm_den;									// Set the initial ECM density to the current ECM density

		fiber.x = 2*(float(rand())/RAND_MAX - 0.5); 
		fiber.y = 2*(float(rand())/RAND_MAX - 0.5); 
		fiber.z = 2*(float(rand())/RAND_MAX - 0.5); 
		
		if (input.izfibflat == 1)
			fiber.z = 0.25*fiber.z;

		if (fiber.norm() != 0.0)
			fiber = fiber/fiber.norm();

		node.collfib.x = fiber.x;
		node.collfib.y = fiber.y;
		node.collfib.z = fiber.z;
				
		if (new_prob == false){
			if ((node.x <= input.ixmin) || (node.x >= input.ixmax))			// If the node lies on a boundary normal to the x-axis...
				enforce_fiber_BCS(node, input, false);									// Enforce the fiber boudary conditions
		
			if ((node.y <= input.iymin) || (node.y >= input.iymax))			// If the node lies on a boundary normal to the y-axis...
				enforce_fiber_BCS(node, input, false);									// Enforce the fiber boundary conditions
		
			if ((node.z <= input.izmin) || (node.z >= input.izmax))			// If the node lies on a boundary normal to the z-axis...
				enforce_fiber_BCS(node, input, false);									// Enforce the fiber boundary conditions
		}

		if (input.icirc == 1){
			vect3 r;
			r.x = node.x - input.ixmax;
			r.y = node.y - input.iymax;
			
			if (fabs(r.norm() - radius) <= 1.0)
				enforce_fiber_BCS(node, input, true);}
				
		node.collfib0.reassign(node.collfib.x, node.collfib.y, node.collfib.z);	

		//if (input.imatx_cnd == 1){
		//	node.theta = 0.;
		//	node.eta = 0.;}
		
		//node.theta0 = node.theta;										// Set initial fiber orientation to current fiber orientation (theta component)
		//node.eta0 = node.eta;											// Set initial fiber orientation to current fiber orientation (eta component)

		input.inodes.push_back(node);									// Add node to the node container within the INPUT class which will be used to construct the grid
	}
	
	if (new_prob)
	{
		//// Read in element connectivity from the FEBio mesh	
		for (int d = 0; d < mesh.Domains(); d++)
		{
			FEDomain& domain = mesh.Domain(d);								// Obtain the domain from FEBio (only one domain)
			int num_elem = domain.Elements();									// Read in the total number of elements from the FEBio domain

			for (int i = 0; i < num_elem; ++i)									// Iterate through all elements...
			{
				input.iNe++;

				FEElement& FEelem = domain.ElementRef(i);						// Obtain element i from the FEBio domain

				vector<int> elem_conn;											// Create the element connectivity vector
				elem_conn.resize(9);										
				elem_conn[0] = i;										// Set the element ID number
				elem_conn[1] = FEelem.m_node[0];								// Identify Node 1 in angio3d (Node 0 in FEBio)
				elem_conn[2] = FEelem.m_node[1];								// Identify Node 2 in angio3d (Node 1 in FEBio)
				elem_conn[3] = FEelem.m_node[3];								// Identify Node 3 in angio3d (Node 3 in FEBio)
				elem_conn[4] = FEelem.m_node[2];								// Identify Node 4 in angio3d (Node 2 in FEBio)
				elem_conn[5] = FEelem.m_node[4];								// Identify Node 5 in angio3d (Node 4 in FEBio)
				elem_conn[6] = FEelem.m_node[5];								// Identify Node 6 in angio3d (Node 5 in FEBio)
				elem_conn[7] = FEelem.m_node[7];								// Identify Node 7 in angio3d (Node 7 in FEBio)								
				elem_conn[8] = FEelem.m_node[6];								// Identify Node 8 in angio3d (Node 6 in FEBio)
            
				input.ieconn.push_back(elem_conn);								// Add the element connectivity vector to the container within the INPUT classe which will be used to construct the angio3d grid 

				//// Gradient
				//double den = 0.;
				//double xpt = 0.;
				//for (int k = 1; k < 9; k++){
				//	xpt = input.inodes[elem_conn[k]].x;
				//	den = 8.0 - (7.0/7500.0)*xpt;
				//	input.inodes[elem_conn[k]].ecm_den = den;
				//	input.inodes[elem_conn[k]].ecm_den0 = den;}
				
				
				//// Microchannels
				//if (d != 0){
				//	for (int k = 1; k < 9; k++){
				//		input.inodes[elem_conn[k]].ecm_den = 3.0;
				//		input.inodes[elem_conn[k]].ecm_den0 = 3.0;
				//		
				//		vec3d cnew;
				//		if (input.inodes[elem_conn[k]].collfib.x >= 0.)
				//			cnew.x = 1.25;
				//		else
				//			cnew.x = -1.25;

				//		cnew.y = input.inodes[elem_conn[k]].collfib.y;
				//		cnew.z = input.inodes[elem_conn[k]].collfib.z;
				//		cnew.unit();
				//		input.inodes[elem_conn[k]].collfib.x = cnew.x;
				//		input.inodes[elem_conn[k]].collfib.y = cnew.y;
				//		input.inodes[elem_conn[k]].collfib.z = cnew.z;
				//	}}

				//// Interface
				//for (int k = 1; k < 9; k++){
				//	double size = 100;

				//	if ((input.inodes[elem_conn[k]].x >= 3750 - size) && (input.inodes[elem_conn[k]].x <= 3750 + size)){	
				//	input.inodes[elem_conn[k]].ecm_den = 6.0;
				//	input.inodes[elem_conn[k]].ecm_den0 = 6.0;}}

				//// Two plugs
				//for (int k = 1; k < 9; k++){
				//	vec3d r;
				//	r.x = 7500 - input.inodes[elem_conn[k]].x;
				//	r.y = 5000 - input.inodes[elem_conn[k]].y;
				//	r.z = 0.;	
				//	double d = r.unit();

				//	if (d < 4000){
				//		input.inodes[elem_conn[k]].ecm_den = 6.0;
				//		input.inodes[elem_conn[k]].ecm_den0 = 6.0;}

				//	r.x = input.inodes[elem_conn[k]].x;
				//	r.y = input.inodes[elem_conn[k]].y;
				//	r.z = 0.;	
				//	d = r.unit();

				//	if (d < 4000){
				//		input.inodes[elem_conn[k]].ecm_den = 6.0;
				//		input.inodes[elem_conn[k]].ecm_den0 = 6.0;}
				//}

				//// One plug
				//if (d != 0){
				//	for (int k = 1; k < 9; k++){
				//		vec3d r;
				//		r.x = 7500 - input.inodes[elem_conn[k]].x;
				//		r.y = 5000 - input.inodes[elem_conn[k]].y;
				//		r.z = 0.;	
				//		r.unit();

				//		vec3d cold;
				//		cold.x = input.inodes[elem_conn[k]].collfib.x;
				//		cold.y = input.inodes[elem_conn[k]].collfib.y;
				//		cold.z = 0.;
				//		
				//		vec3d cnew;

				//		cnew = cold - r*(r*cold);
				//		cnew.z = 0.;
				//		cnew.unit();

				//		input.inodes[elem_conn[k]].ecm_den = 8.0;
				//		input.inodes[elem_conn[k]].ecm_den0 = 8.0;
				//		input.inodes[elem_conn[k]].collfib.x = cnew.x;
				//		input.inodes[elem_conn[k]].collfib.y = cnew.y;
				//		input.inodes[elem_conn[k]].collfib.z = cnew.z;
				//	}}

			}
		}
	}
	else
	{
		//// read in element connectivity from the febio mesh	
		FEDomain& Domain = mesh.Domain(0);								// obtain the domain from febio (only one domain)
		input.iNe = Domain.Elements();									// read in the total number of elements from the febio domain

		for (int i = 0; i < input.iNe; ++i)									// iterate through all elements...
		{
			FEElement& FEElem = Domain.ElementRef(i);						// obtain element i from the febio domain

			vector<int> elem_conn;											// create the element connectivity vector
			elem_conn.resize(9);										
			elem_conn[0] = i;												// set the nodes id number
			elem_conn[1] = FEElem.m_node[0];								// identify node 1 in angio3d (node 0 in febio)
			elem_conn[2] = FEElem.m_node[1];								// identify node 2 in angio3d (node 1 in febio)
			elem_conn[3] = FEElem.m_node[3];								// identify node 3 in angio3d (node 3 in febio)
			elem_conn[4] = FEElem.m_node[2];								// identify node 4 in angio3d (node 2 in febio)
			elem_conn[5] = FEElem.m_node[4];								// identify node 5 in angio3d (node 4 in febio)
			elem_conn[6] = FEElem.m_node[5];								// identify node 6 in angio3d (node 5 in febio)
			elem_conn[7] = FEElem.m_node[7];								// identify node 7 in angio3d (node 7 in febio)								
			elem_conn[8] = FEElem.m_node[6];								// identify node 8 in angio3d (node 6 in febio)
            
			input.ieconn.push_back(elem_conn);								// add the element connectivity vector to the container within the input classe which will be used to construct the angio3d grid 
		}
	}
	
	for (int i = 0; i < input.iNn; i++){
		input.inodes[i].ecm_den_store.push_back(input.inodes[i].ecm_den0);
		input.inodes[i].ecm_fibril_store.push_back(input.inodes[i].collfib0);}


	//// Determine boundary faces for each element
	int BC_violate = 0;												// Create and initialize the boundary inidicator
	int BC_count = 0;												// Create and initialize the boundary element counter			

	vector<int> elem_conn;											// Create the temporary element connectivity vector

	for (int i = 0; i < input.iNe; ++i)								// Iterate through all elements...
	{
		elem_conn = input.ieconn[i];									// Obtain element connectivity vector for element i

		vector<int> eBC_violate;										// Create the boundary face vector
        eBC_violate.resize(7);
        
		eBC_violate[0] = elem_conn[0];									// Set the element ID within the boundary face vector

		// 
		// If face 1 lies on the boundary of the domain...
		if ((input.inodes[elem_conn[1]].y <= input.iymin) && (input.inodes[elem_conn[2]].y <= input.iymin) && (input.inodes[elem_conn[5]].y <= input.iymin) && (input.inodes[elem_conn[6]].y <= input.iymin))
		{
			BC_violate = 1;												// Turn on the boundary indicator
		    eBC_violate[1] = 1;											// Indicate that face 1 of this element lies at the minimum in the y-direction
		}
		
		// If face 2 lies on the boundary of the domain...
		if ((input.inodes[elem_conn[2]].x >= input.ixmax) && (input.inodes[elem_conn[4]].x >= input.ixmax) && (input.inodes[elem_conn[6]].x >= input.ixmax) && (input.inodes[elem_conn[8]].x >= input.ixmax))
		{
			BC_violate = 1;												// Turn on the boundary indicator
		    eBC_violate[2] = 1;											// Indicate that face 2 of this element lies at the maximum in the x-direction 
		}
		
		// If face 3 lies on the boundary of the domain...
		if ((input.inodes[elem_conn[3]].y >= input.iymax) && (input.inodes[elem_conn[4]].y >= input.iymax) && (input.inodes[elem_conn[7]].y >= input.iymax) && (input.inodes[elem_conn[8]].y >= input.iymax))
		{
			BC_violate = 1;												// Turn on the boundary indicator
		    eBC_violate[3] = 1;											// Indicate that face 3 of this element lies at the maximum in the y-direction
		}
		
		// If face 4 lies on the boundary of the domain...
		if ((input.inodes[elem_conn[1]].x <= input.ixmin) && (input.inodes[elem_conn[3]].x <= input.ixmin) && (input.inodes[elem_conn[5]].x <= input.ixmin) && (input.inodes[elem_conn[7]].x <= input.ixmin))
		{
			BC_violate = 1;												// Turn on the boundary indicator
		    eBC_violate[4] = 1;											// Indicate that face 4 of this element lies at the minimum in the x-direction
		}

		// If face 5 lies on the boundary of the domain...
		if ((input.inodes[elem_conn[5]].z >= input.izmax) && (input.inodes[elem_conn[6]].z >= input.izmax) && (input.inodes[elem_conn[7]].z >= input.izmax) && (input.inodes[elem_conn[8]].z >= input.izmax))
		{
			BC_violate = 1;												// Turn on the boundary indicator
		    eBC_violate[5] = 1;											// Indicate that face 5 of this element lies at the maximum in the z-direction
		}
		
		// If face 6 lies on the boundary of the domain...
		if ((input.inodes[elem_conn[1]].z <= input.izmin) && (input.inodes[elem_conn[2]].z <= input.izmin) && (input.inodes[elem_conn[3]].z <= input.izmin) && (input.inodes[elem_conn[4]].z <= input.izmin))
		{
			BC_violate = 1;												// Turn on the boundar indicator
		    eBC_violate[6] = 1;											// Indicate that face 6 of this element lies at the minimum in the z-direction
		}
		
		if (BC_violate == 1)											// If the boundary indicator is on...
		{
			input.ieBC.push_back(eBC_violate);								// Add the element's boundary face vector to the contianer in the INPUT class 				
			BC_count++;														// Increment the boundary element counter
		}

		BC_violate = 0;													// Turn off the boundary indicator
	}
	
	input.iNBC = BC_count;											// Set the total number of elements at the boundaries
		
	input.igrid_in = 'y';											// Indicate that the grid has been read in from the FEBio mesh

	return;
}



///////////////////////////////////////////////////////////////////////
// enfore_fiber_BCS
///////////////////////////////////////////////////////////////////////

void Filein::enforce_fiber_BCS(Node &node, FEAngInput &input, bool circ)
{
	if (circ == true){
		vect3 r;
		r.x = node.x - input.ixmax;
		r.y = node.y - input.iymax;

		r = r/r.norm();

		mat3 I; I.M11 = 1.; I.M22 = 1.; I.M33 = 1.;
		mat3 M; M.M11 = r.x*r.x; M.M12 = r.x*r.y; M.M13 = r.x*r.z; M.M21 = r.y*r.x; M.M22 = r.y*r.y; M.M23 = r.y*r.z; M.M31 = r.z*r.x; M.M32 = r.z*r.y; M.M33 = r.z*r.z;

		vect3 r_new; r_new = (I - M)*node.collfib; r_new = r_new/r_new.norm();

		node.collfib.x = r_new.x; node.collfib.y = r_new.y; node.collfib.z = r_new.z;

		return;
	}
		
	if (input.igelbc == 'l'){
		//Node on the front face
		if (node.y == input.iymin)
			node.collfib.y = 0;

		// Node on the back face
		//if (node.y == input.iymax)
			//node.collfib.y = 0;

		// Node on the bottom face
		if (node.z == input.izmin)
			node.collfib.z = 0;

		// Node on the top face
		//if (node.z == input.izmax)
			//node.collfib.z = 0;
	}

	if (input.igelbc == 's'){
		// Node on the right face
		//if (node.x == input.ixmax)
			//node.collfib.x = 0;

		// Node on the left face
		if (node.x == input.ixmin)
			node.collfib.x = 0;

		// Node on the bottom face
		if (node.z == input.izmin)
			node.collfib.z = 0;

		// Node on the top face
		//if (node.z == input.izmax)
			//node.collfib.z = 0;
	}

	if (input.igelbc == 'u'){
		//Node on the front face
		if (node.y == input.iymin)
			node.collfib.y = 0;

		// Node on the right face
		//if (node.x == input.ixmax)
			//node.collfib.x = 0;

		// Node on the back face
		//if (node.y == input.iymax)
			//node.collfib.y = 0;

		// Node on the left face
		if (node.x == input.ixmin)
			node.collfib.x = 0;

		// Node on the bottom face
		if (node.z == input.izmin)
			node.collfib.z = 0;

		// Node on the top face
		//if (node.z == input.izmax)
			//node.collfib.z = 0;
	}

	if (input.igelbc == 'n'){
	}

	return;
}







