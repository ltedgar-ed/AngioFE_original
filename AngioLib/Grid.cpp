///////////////////////////////////////////////////////////////////////
// Grid.cpp
///////////////////////////////////////////////////////////////////////



// Include:
#include "stdafx.h"

#include "Grid.h"
#include "Data.h"
#include "angio3d.h"
#include "Input.h"
#include "vect3.h"
#include "Elem.h"
#include "math.h"
///////////////////////////////////////////////////////////////////////
// Local Functions
///////////////////////////////////////////////////////////////////////

// GRID.factorial - Takes an integer and returns the factorial value

int factorial (int num)
{
	int result=1;
	for (int i=1; i<=num; ++i)
		result*=i;
	return result;
}



///////////////////////////////////////////////////////////////////////
// Construction/Destruction
///////////////////////////////////////////////////////////////////////
Grid::Grid(Input &input)                      // Constructor for GRID object
{                                                               
	create_grid(input);        
}


Grid::~Grid()                                                   // Destructor for GRID object
{
	//delete [] density;
}



///////////////////////////////////////////////////////////////////////
// Member Functions
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// create_grid
///////////////////////////////////////////////////////////////////////

void Grid::create_grid(Input &input)
{
	int i = 0, j = 0, k = 0;                                          // Declare book keeping indices i, j, and k
	
	coll_den = static_cast<double>(input.imatx_den);            // Determine collagen density from the input file
	load_cond = input.imatx_cnd;                                // Determine matrix condition from the input file
	
	xnodes = input.ixnodes;                                     // Number of nodes in the x direction
	ynodes = input.iynodes;                                     // Number of nodes in the y direction
	znodes = input.iznodes;                                     // Number of nodes in the z direction
    Nn = xnodes*ynodes*znodes;
    Ne = (xnodes - 1)*(ynodes - 1)*(znodes - 1);            
    
	///// Define grid dimensions
	xrange[0] = static_cast<double>(input.ixmin);               // Minimum distance in x direction 
	xrange[1] = static_cast<double>(input.ixmax);               // Maximum distance in x direction
	yrange[0] = static_cast<double>(input.iymin);               // Minimum distance in y direction
	yrange[1] = static_cast<double>(input.iymax);               // Maximum distance in y direction
	zrange[0] = static_cast<double>(input.izmin);               // Minimum distance in z direction
	zrange[1] = static_cast<double>(input.izmax);               // Maximum distance in z direction

	///// Define number of nodes in each direction
	num_nodes[0] = xnodes;                                      // Number of nodes in x direction                                      
	num_nodes[1] = ynodes;                                      // Number of nodes in y direction
	num_nodes[2] = znodes;                                      // Number of nodes in z direction
    
    
    ///// Resize arrays containing the nodal positions
	x.resize(xnodes);                                           // Nodal positions in the x direction
	y.resize(ynodes);                                           // Nodal positions in the y direction
	z.resize(znodes);                                           // Nodal positions in the z direction

	//x_bctype = input.ixbctype;
	//y_bctype = input.iybctype;
	//z_bctype = input.izbctype;
	
	frontbc = input.ifrontbc;
	rightbc = input.irightbc;
	backbc = input.ibackbc;
	leftbc = input.ileftbc;
	bottombc = input.ibottombc;
	topbc = input.itopbc;
	
	/////// Find collagen density scaling factor
	//a = -2.0;
	//b = 5.184;
	//c = 0.1823;
	//den_scale = a+b*pow(E,-c*coll_den );                        // den_scale - Determine the density scaling factor using the function defined by a, b, c	 
	
	
	/////  If the grid is input by the user...
	if (input.igrid_in == 'y'){
	    	    
	    Ne = input.iNe;
	    Nn = input.iNn;
	    nodes = input.inodes;
	    
	  //  for (i = 0; i < Nn; ++i){
	  //      nodes[i].ecm_den = coll_den;
			//nodes[i].ecm_den0 = coll_den;
   //         /*nodes[i].theta = float(rand())/RAND_MAX*2*pi-pi;
   //         nodes[i].eta = float(rand())/RAND_MAX*2*pi-pi;*/}	    
	  //  	    
	    int n1, n2, n3, n4, n5, n6, n7, n8;
	    
	    for (i = 0; i < Ne; ++i){
	        Elem elem;
	        elem.elem_num = input.ieconn[i][0];
	        n1 = input.ieconn[i][1];
	        n2 = input.ieconn[i][2];
	        n3 = input.ieconn[i][3];
	        n4 = input.ieconn[i][4];
	        n5 = input.ieconn[i][5];
	        n6 = input.ieconn[i][6];
	        n7 = input.ieconn[i][7];
	        n8 = input.ieconn[i][8];
	            
	        elem.n1 = &nodes[n1];
	        elem.n2 = &nodes[n2];
	        elem.n3 = &nodes[n3];
	        elem.n4 = &nodes[n4];
	        elem.n5 = &nodes[n5];
	        elem.n6 = &nodes[n6];
	        elem.n7 = &nodes[n7];
	        elem.n8 = &nodes[n8];
	        
			ebin.push_back(elem);}   
	    
	    
	    int elem_num = 0;
	    
	    for (i = 0; i < input.iNBC; ++i){
	        elem_num = input.ieBC[i][0];
	        
	        if (input.ieBC[i][1] == 1){
	            ebin[elem_num].f1.BC = true;
	            //ebin[elem_num].f1.bc_type = y_bctype;}
				ebin[elem_num].f1.bc_type = frontbc;}
	            
	        if (input.ieBC[i][2] == 1){
	            ebin[elem_num].f2.BC = true;
	            //ebin[elem_num].f2.bc_type = x_bctype;}
	            ebin[elem_num].f2.bc_type = rightbc;}

	        if (input.ieBC[i][3] == 1){
	            ebin[elem_num].f3.BC = true;
	            //ebin[elem_num].f3.bc_type = y_bctype;}
	            ebin[elem_num].f3.bc_type = backbc;}

	        if (input.ieBC[i][4] == 1){
	            ebin[elem_num].f4.BC = true;
	            //ebin[elem_num].f4.bc_type = x_bctype;}
				ebin[elem_num].f4.bc_type = leftbc;}
	            
	        if (input.ieBC[i][5] == 1){
	            ebin[elem_num].f5.BC = true;
	            //ebin[elem_num].f5.bc_type = z_bctype;}
	            ebin[elem_num].f5.bc_type = topbc;}

	        if (input.ieBC[i][6] == 1){
	            ebin[elem_num].f6.BC = true;    
	            //ebin[elem_num].f6.bc_type = z_bctype;}}
				ebin[elem_num].f6.bc_type = bottombc;}}
   } 
	            
	            
	
	///// ... else construct the default regular hex grid.
	else{
	    
	    
	    // Create the nodal positions in the x, y, z direction
	    for (i = 0; i < num_nodes[0]; ++i)
		    x[i] = xrange[0] + float(i)/(num_nodes[0]-1)*(xrange[1]-xrange[0]);
    		
	    for (j = 0; j < num_nodes[1]; ++j)
		    y[j] = yrange[0] + float(j)/(num_nodes[1]-1)*(yrange[1]-yrange[0]);
    	
	    for (k = 0; k < num_nodes[2]; ++k)
		    z[k] = zrange[0] + float(k)/(num_nodes[2]-1)*(zrange[1]-zrange[0]);
    	  
	    
	        
	    // Define arrays theta and eta, which contain information about the orientation of collagen fibers
	    for (i = 0; i < xnodes; i++){
	        theta.resize(xnodes);
	        eta.resize(xnodes);
	        for (j = 0; j < ynodes; j++){
	            theta[i].resize(ynodes);
	            eta[i].resize(ynodes);
	            for (k = 0; k < znodes; k++){
	            theta[i][j].assign(znodes,0.);
	            eta[i][j].assign(znodes,0.);}}}
	    
		for (i = 0; i < xnodes; i++){
	        coll_fib_x.resize(xnodes);
	        coll_fib_y.resize(xnodes);
			coll_fib_z.resize(xnodes);
	        for (j = 0; j < ynodes; j++){
	            coll_fib_x[i].resize(ynodes);
	            coll_fib_y[i].resize(ynodes);
				coll_fib_z[i].resize(ynodes);
	            for (k = 0; k < znodes; k++){
					coll_fib_x[i][j].assign(znodes,0.);
					coll_fib_y[i][j].assign(znodes,0.);
					coll_fib_z[i][j].assign(znodes,0.);}}}

		// Define collagen fiber orientation
	    switch (load_cond)                                          // Switch between the two loading opitions:
	    {
	        case 0:                                                 //  Case 1: Randomly generate fiber orientation between -pi and pi 
		        for (i = 0; i < num_nodes[0]; ++i)
		        {
			        for (j = 0; j < num_nodes[1]; ++j)
			        {
				        for (k=0; k < num_nodes[2]; ++k)
				        {
					        if (i < num_nodes[0]-1 && j < num_nodes[1]-1 && k < num_nodes[2]-1)
					        {
						        theta[i][j][k] = float(rand())/RAND_MAX*2*pi-pi;                        
						        eta[i][j][k] = float(rand())/RAND_MAX*2*pi-pi;
								coll_fib_x[i][j][k] = 2*(float(rand())/RAND_MAX - 0.5); 
								coll_fib_y[i][j][k] = 2*(float(rand())/RAND_MAX - 0.5); 
								coll_fib_z[i][j][k] = 2*(float(rand())/RAND_MAX - 0.5); 
								//coll_fib_z[i][j][k] = 0;
							}
					        else if (i == num_nodes[0]-1 && j < num_nodes[1]-1 && k < num_nodes[2]-1)
					        {
						        theta[i][j][k] = theta[0][j][k];
						        eta[i][j][k] = eta[0][j][k];						
							}
					        else if (j == num_nodes[1]-1 && i < num_nodes[0]-1 && k < num_nodes[2]-1)
					        {
						        theta[i][j][k] = theta[i][0][k];
						        eta[i][j][k] = eta[i][0][k];
					        }
					        else if (k == num_nodes[2]-1 && i < num_nodes[0]-1 && j < num_nodes[1]-1)
					        {
						        theta[i][j][k] = theta[i][j][0];
						        eta[i][j][k] = eta[i][j][0];
					        }
					        else
					        {
						        theta[i][j][k] = theta[0][0][0];
						        eta[i][j][k] = eta[0][0][0];
					        }
				        }
			        }
		        }
		        break;
    	
	        case 1:                                                 // Case 2: Orient collagen fibers to align along the x-axis
		        for (i = 0; i < num_nodes[0]; ++i)
		        {
			        for (j = 0; j < num_nodes[1]; ++j)
			        {
				        for (k=0; k< num_nodes[2]; ++k)
				        {
					        theta[i][j][k] = 0;
					        eta[i][j][k] = 0;
				        }
			        }
		        }
		        break;
    		    
	        case 2:                                                 // Case 2: Orient collagen fibers to align along the x-axis
		        for (i = 0; i < num_nodes[0]; ++i)
		        {
			        for (j = 0; j < num_nodes[1]; ++j)
			        {
				        for (k=0; k< num_nodes[2]; ++k)
				        {
					        theta[i][j][k] = pi/2;
					        eta[i][j][k] = 0;
				        }
			        }
		        }
		        break;
	    }
    	       
	    
		// Create the list of nodal coordinates
	    vect3 fiber;

		int nnumber = 0;
	    nodes.resize(Nn);
    	
	    for (i = 0; i < num_nodes[0]; ++i){
	        for (j = 0; j < num_nodes[1]; ++j){
	            for (k = 0; k < num_nodes[2]; ++k){ 
	            nodes[nnumber].id = nnumber;
	            nodes[nnumber].x = x[i];
	            nodes[nnumber].y = y[j];
	            nodes[nnumber].z = z[k];
				
				nodes[nnumber].x0 = nodes[nnumber].x;
				nodes[nnumber].y0 = nodes[nnumber].y;
				nodes[nnumber].z0 = nodes[nnumber].z;

				nodes[nnumber].theta = theta[i][j][k];
				nodes[nnumber].eta = eta[i][j][k];	

				fiber.x = 2*(float(rand())/RAND_MAX - 0.5); 
				fiber.y = 2*(float(rand())/RAND_MAX - 0.5); 
				fiber.z = 2*(float(rand())/RAND_MAX - 0.5); 
				
				if (fiber.norm() != 0.0)
					fiber = fiber/fiber.norm();

				nodes[nnumber].collfib.x = fiber.x;
				nodes[nnumber].collfib.y = fiber.y;
				nodes[nnumber].collfib.z = fiber.z;
				//nodes[nnumber].collfib.z = 0.5*fiber.z;
				nodes[nnumber].collfib0.reassign(fiber.x, fiber.y, fiber.z);

				coll_fib_x[i][j][k] = fiber.x;
				coll_fib_y[i][j][k] = fiber.y;
				coll_fib_z[i][j][k] = fiber.z;

				nodes[nnumber].ecm_den = coll_den;
				nodes[nnumber].ecm_den0 = coll_den;
				
				//nodes[nnumber].ecm_den = input.icoll_den[nnumber];
				//nodes[nnumber].ecm_den0 = input.icoll_den[nnumber];

				nnumber++;}}} 	   
		
		// Create the bin that contains all the elements within the grid
	    int ecount = 0;
	        	
	    for (i = 0; i < num_nodes[0] - 1; ++i){
	        for (j = 0; j < num_nodes[1] - 1; ++j){
	            for (k = 0; k < num_nodes[2] - 1; ++k){         
	                		
					Elem elem;
	                Node testnode;
					
					testnode.x = x[i];
	                testnode.y = y[j];
	                testnode.z = z[k];
	                elem.n1 = &nodes[findnnumber(testnode)];
	                
					testnode.x = x[i+1];
	                testnode.y = y[j];
	                testnode.z = z[k];
	                elem.n2 = &nodes[findnnumber(testnode)];

					testnode.x = x[i];
	                testnode.y = y[j+1];
	                testnode.z = z[k];
	                elem.n3 = &nodes[findnnumber(testnode)];

					testnode.x = x[i+1];
	                testnode.y = y[j+1];
	                testnode.z = z[k];
	                elem.n4 = &nodes[findnnumber(testnode)];
               
					testnode.x = x[i];
	                testnode.y = y[j];
	                testnode.z = z[k+1];
	                elem.n5 = &nodes[findnnumber(testnode)];
	                
					testnode.x = x[i+1];
	                testnode.y = y[j];
	                testnode.z = z[k+1];
	                elem.n6 = &nodes[findnnumber(testnode)];

					testnode.x = x[i];
	                testnode.y = y[j+1];
	                testnode.z = z[k+1];
	                elem.n7 = &nodes[findnnumber(testnode)];

					testnode.x = x[i+1];
	                testnode.y = y[j+1];
	                testnode.z = z[k+1];
	                elem.n8 = &nodes[findnnumber(testnode)];
    	            
	                if (i == 0){
	                    elem.f4.BC = true;
	                    elem.f4.bc_type = leftbc;}
    	                
	                if (i == num_nodes[0] - 2){
	                    elem.f2.BC = true;
	                    elem.f2.bc_type = rightbc;}
    	                
	                if (j == 0){
	                    elem.f1.BC = true;
	                    elem.f1.bc_type = frontbc;}
    	                
	                if (j == num_nodes[1] - 2){
	                    elem.f3.BC = true;
	            	    elem.f3.bc_type = backbc;}
    	            	
	                if (k == 0){
	                    elem.f6.BC = true;
	                    elem.f6.bc_type = topbc;}
    	                
	                if (k == num_nodes[2] - 2){
	                    elem.f5.BC = true;
	                    elem.f5.bc_type = bottombc;}
    	                            
	                elem.elem_num = ecount;
	                ebin.push_back(elem);
	                ++ecount;}}}	
	}    

	return;
}

///////////////////////////////////////////////////////////////////////
// findnnumber
///////////////////////////////////////////////////////////////////////

int Grid::findnnumber(Node node)
{
    int nnumber = -1;
    Node testnode;
    
    for (int i = 0; i < Nn; ++i){
        testnode = nodes[i];
        
        if ((testnode.x == node.x) && (testnode.y == node.y) && (testnode.z == node.z)){
            nnumber = i;
            break;}
    }      
       
    return nnumber;
}



///////////////////////////////////////////////////////////////////////
// ijk
///////////////////////////////////////////////////////////////////////

// GRID.ijk - Determines the i,j,k nodes that describes the element that contains a given position in global coorindates
//      Input:  - Position in global coorindates (xpt, ypt, zpt)
//              - Reference to I, J, K integers
//              - DATA structure
//
//      Output: - None (operates on references)

void Grid::ijk(double xpt, double ypt, double zpt, int &I, int &J, int &K, Data &data)
{
    // I, J, and K describe which grid element the vessel tip is in. 
    // In particular, they refer to the lower nodes of the grid in the x, y, and z direction respectively.
	// The 8 nodes of any hex element are (I, J, K), (I+1, J, K), (I, J+1, K), (I+1, J+1, K),
    //                                    (I, J, K+1), (I+1, J, K+1), (I, J+1, K+1), and (I+1, J+1, K+1).
    
	for (int i = 1; i < num_nodes[0]; ++i)
        if (x[i] > xpt){
             I = i - 1;
             break;}
             
    for (int j = 1; j < num_nodes[1]; ++j)
        if (y[j] > ypt){
             J = j - 1;
             break;}
             
    for (int k = 1; k < num_nodes[2]; ++k)
        if (z[k] > zpt){
             K = k - 1;
             break;}     
    
    I=(int)floor((xpt-xrange[0])/data.dx);                           // I, J, and K describe which grid element the vessel tip is in
	J=(int)floor((ypt-yrange[0])/data.dy);                           // In particular, they refer to the lower nodes of the grid in the
	K=(int)floor((zpt-zrange[0])/data.dz);                           // x, y, and z direction respectively.
	                                                            // The 8 nodes of the element are (I, J, K), (I+1, J, K), (I, J+1, K), (I+1, J+1, K),
	                                                            //                              (I, J, K+1), (I+1, J, K+1), (I, J+1, K+1), and (I+1, J+1, K+1).

	if (I >= num_nodes[0]-2)
		I = num_nodes[0]-2;
    else if (I<0)
		I=0;

	if (J >= num_nodes[1]-2)
		J = num_nodes[1]-2;
	else if (J<0)
		J=0;
	
	if (K > num_nodes[2]-2)
		K = num_nodes[2]-2;
	else if (K<0)
		K=0;
                  
    return;
}



///////////////////////////////////////////////////////////////////////
// findelem
///////////////////////////////////////////////////////////////////////

int Grid::findelem(double xpt, double ypt, double zpt)
{
    int noBC = 0;
    int elem_num = -1;
    double bb_eps = 0.1;
	double eps = 0;
    
    double xmin, xmax, ymin, ymax, zmin, zmax = {0.0};
    double xix, xiy, xiz = {0.0};
    Elem elem;
    
    for (int i = 0; i < Ne; i++)
    {
        elem = ebin[i];
        
        noBC = 0;
                        
        xmin = elem.bb_xmin()*(1. - bb_eps);
        xmax = elem.bb_xmax()*(1. + bb_eps);
		if (xmin == 0.0)
			xmin = -1.0*bb_eps*xmax;
		
		ymin = elem.bb_ymin()*(1. - bb_eps);
        ymax = elem.bb_ymax()*(1. + bb_eps);
        if (ymin == 0.0)
			ymin = -1.0*bb_eps*ymax;

		zmin = elem.bb_zmin()*(1. - bb_eps);
        zmax = elem.bb_zmax()*(1. + bb_eps);
        if (zmin == 0.0)
			zmin = -1.0*bb_eps*zmax;

        
		if ((xpt >= xmin) && (xpt <= xmax)){
            if ((ypt >= ymin) && (ypt <= ymax)){
                if ((zpt >= zmin) && (zpt <= zmax)){
                    natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                    
                    if ((fabs(xix) <= (1.0 + eps)) && (fabs(xiy) <= (1.0 + eps)) && (fabs(xiz) <= (1.0 + eps))){
                        elem_num = i;
						return elem_num;}}}}
	}	  
           
    return elem_num;
}

int Grid::enhanced_findelem(double xpt, double ypt, double zpt)
{
    int noBC = 0;
    int elem_num = 0;
    double bb_eps = 0.001;
	double eps = 0.000001;
    
    double xmin, xmax, ymin, ymax, zmin, zmax = {0.0};
    double xix, xiy, xiz = {0.0};
    Elem elem;
    
    for (int i = 0; i < Ne; i++)
    {
        elem = ebin[i];
        
        noBC = 0;
        
        if ((elem.f1.BC == true) || (elem.f2.BC == true) || (elem.f3.BC == true) || (elem.f4.BC == true) || (elem.f5.BC == true) || (elem.f6.BC == true))
            noBC = 1;
                
        xmin = elem.bb_xmin()*(1 - bb_eps);
        xmax = elem.bb_xmax()*(1 + bb_eps);
		if (xmin == 0.0)
			xmin = -1.0*bb_eps*xmax;
		
		ymin = elem.bb_ymin()*(1 - bb_eps);
        ymax = elem.bb_ymax()*(1 + bb_eps);
        if (ymin == 0.0)
			ymin = -1.0*bb_eps*ymax;

		zmin = elem.bb_zmin()*(1 - bb_eps);
        zmax = elem.bb_zmax()*(1 + bb_eps);
        if (zmin == 0.0)
			zmin = -1.0*bb_eps*zmax;

        if (noBC == 0){
            
            if ((xpt >= xmin) && (xpt <= xmax)){
                if ((ypt >= ymin) && (ypt <= ymax)){
                    if ((zpt >= zmin) && (zpt <= zmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                    
                        if ((fabs(xix) <= 1.0 + eps) && (fabs(xiy) <= 1.0 + eps) && (fabs(xiz) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}}
                            
        else if (noBC == 1){
        
            if ((xpt >= xmin) && (xpt <= xmax)){
                if ((ypt >= ymin) && (ypt <= ymax)){
                    if ((zpt >= zmin) && (zpt <= zmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                    
                        if ((fabs(xix) <= 1.0 + eps) && (fabs(xiy) <= 1.0 + eps) && (fabs(xiz) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}
                            
            // face 1
            if (elem.f1.BC == true){
                if ((xpt >= xmin) && (xpt <= xmax)){
                    if (ypt <= ymin){
                        if ((zpt >= zmin) && (zpt <= zmax)){
                            natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                        
                            if ((fabs(xix) <= 1.0 + eps) && (xiy <= -1.0) && (fabs(xiz) <= 1.0 + eps)){
                                elem_num = i;
                                break;}}}}}
            
            // face 2
            if (elem.f2.BC == true){
                if (xpt >= xmax){
                    if ((ypt >= ymin) && (ypt <= ymax)){
                        if ((zpt >= zmin) && (zpt <= zmax)){
                            natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                        
                            if ((xix >= 1.0) && (fabs(xiy) <= 1.0 + eps) && (fabs(xiz) <= 1.0 + eps)){
                                elem_num = i;
                                break;}}}}}
                                   
            // face 3
            if (elem.f3.BC == true){
                if ((xpt >= xmin) && (xpt <= xmax)){
                    if (ypt >= ymax){
                        if ((zpt >= zmin) && (zpt <= zmax)){
                            natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                        
                            if ((fabs(xix) <= 1.0 + eps) && (xiy >= 1.0) && (fabs(xiz) <= 1.0 + eps)){
                                elem_num = i;
                                break;}}}}}
            
            // face 4
            if (elem.f4.BC == true){
                if (xpt <= xmin){
                    if ((ypt >= ymin) && (ypt <= ymax)){
                        if ((zpt >= zmin) && (zpt <= zmax)){
                            natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                        
                            if ((xix <= -1.0) && (fabs(xiy) <= 1.0 + eps) && (fabs(xiz) <= 1.0 + eps)){
                                elem_num = i;
                                break;}}}}}
            
            // face 5
            if (elem.f5.BC == true){
                if ((xpt >= xmin) && (xpt <= xmax)){
                    if ((ypt >= ymin) && (ypt <= ymax)){
                        if (zpt >= zmax){
                            natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                        
                            if ((fabs(xix) <= 1.0 + eps) && (fabs(xiy) <= 1.0 + eps) && (xiz >= 1.0)){
                                elem_num = i;
                                break;}}}}}
            
            // face 6
            if (elem.f6.BC == true){
                if ((xpt >= xmin) && (xpt <= xmax)){
                    if ((ypt >= ymin) && (ypt <= ymax)){
                        if (zpt <= zmin){
                            natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);
                        
                            if ((fabs(xix) <= 1.0 + eps) && (fabs(xiy) <= 1.0 + eps) && (xiz <= -1.0)){
                                elem_num = i;
                                break;}}}}}
            
            // 1 and 2
            if ((elem.f1.BC == true) && (elem.f2.BC == true)){
                if ((ypt <= ymin) && (xpt >= xmax)){
                    if ((zpt >= zmin) && (zpt <= zmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiy <= -1.0) && (xix >= 1.0)) && (fabs(xiz) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}} 
            
            // 2 and 3
            if ((elem.f2.BC == true) && (elem.f3.BC == true)){
                if ((xpt >= xmax) && (ypt >= ymax)){
                    if ((zpt >= zmin) && (zpt <= zmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xix >= 1.0) && (xiy >= 1.0)) && (fabs(xiz) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}                
            
            // 3 and 4
            if ((elem.f3.BC == true) && (elem.f4.BC == true)){
                if ((ypt >= ymax) && (xpt <= xmin)){
                    if ((zpt >= zmin) && (zpt <= zmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiy >= 1.0) && (xix <= -1.0)) && (fabs(xiz) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}
            
            // 4 and 1
            if ((elem.f4.BC == true) && (elem.f1.BC == true)){
                if ((xpt <= xmin) && (ypt <= ymin)){
                    if ((zpt >= zmin) && (zpt <= zmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xix <= -1.0) && (xiy <= -1.0)) && (fabs(xiz) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}
            
            // 5 and 1
            if ((elem.f5.BC == true) && (elem.f1.BC == true)){
                if ((zpt >= zmax) && (ypt <= ymin)){
                    if ((xpt >= xmin) && (xpt <= xmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiz >= 1.0) && (xiy <= -1.0)) && (fabs(xix) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}} 
            
            // 5 and 2
            if ((elem.f5.BC == true) && (elem.f2.BC == true)){
                if ((zpt >= zmax) && (xpt >= xmax)){
                    if ((ypt >= ymin) && (ypt <= ymax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiz >= 1.0) && (xix >= 1.0)) && (fabs(xiy) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}           
            
            // 5 and 3
            if ((elem.f5.BC == true) && (elem.f3.BC == true)){
                if ((zpt >= zmax) && (ypt >= ymax)){
                    if ((xpt >= xmin) && (xpt <= xmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiz >= 1.0) && (xiy >= 1.0)) && (fabs(xix) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}} 
                        
            // 5 and 4
            if ((elem.f5.BC == true) && (elem.f4.BC == true)){
                if ((zpt >= zmax) && (xpt <= xmin)){
                    if ((ypt >= ymin) && (ypt <= ymax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiz >= 1.0) && (xix <= -1.0)) && (fabs(xiy) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}       
            // 6 and 1
            if ((elem.f6.BC == true) && (elem.f1.BC == true)){
                if ((zpt <= zmin) && (ypt <= ymin)){
                    if ((xpt >= xmin) && (xpt <= xmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiz <= -1.0) && (xiy <= -1.0)) && (fabs(xix) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}} 
            
            // 6 and 2
            if ((elem.f6.BC == true) && (elem.f2.BC == true)){
                if ((zpt <= zmin) && (xpt >= xmax)){
                    if ((ypt >= ymin) && (ypt <= ymax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiz <= -1.0) && (xix >= 1.0)) && (fabs(xiy) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}           
            
            // 6 and 3
            if ((elem.f6.BC == true) && (elem.f3.BC == true)){
                if ((zpt <= zmin) && (ypt >= ymax)){
                    if ((xpt >= xmin) && (xpt <= xmax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiz <= -1.0) && (xiy >= 1.0)) && (fabs(xix) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}} 
                        
            // 6 and 4
            if ((elem.f6.BC == true) && (elem.f4.BC == true)){
                if ((zpt <= zmin) && (xpt <= xmin)){
                    if ((ypt >= ymin) && (ypt <= ymax)){
                        natcoord(xix, xiy, xiz, xpt, ypt, zpt, i);    
        
                        if (((xiz <= -1.0) && (xix <= -1.0)) && (fabs(xiy) <= 1.0 + eps)){
                            elem_num = i;
                            break;}}}}
                                            
        }                   
    }      
           
    return elem_num;
}










///////////////////////////////////////////////////////////////////////
// natcoordinates
///////////////////////////////////////////////////////////////////////

// GRID.natcoorindates - Converts a global position into a position in natural coorindates for a particular element
//      Input:  - Position in global coorindates (xpt, ypt, zpt)
//              - Nodes describing the current element (I, J, K)
//              - Reference to position in natural coorindates (xix, xiy, xiz)
//
//      Output: - None (operates on references)

void Grid::natcoordinates(double &xix, double &xiy, double &xiz, double xpt, double ypt, double zpt, int I, int J, int K)
{
    if (xpt > xrange[1])
       	xpt = xrange[1];
    else if (xpt < xrange[0])
    	xpt = xrange[0];
	
    if (ypt > yrange[1])
    	ypt = yrange[1];
    else if (ypt < yrange[0])
    	ypt = yrange[0];
	
    if (zpt > zrange[1])
    	zpt = zrange[1];
    else if (zpt < zrange[0])
    	zpt = zrange[0];
    
    ///// Transform global position to natural coorindates (ranging from -1 to 1):
    xix = (xpt - (x[I] + x[I+1])/2)/((x[I+1] - x[I])/2); 
    xiy = (ypt - (y[J] + y[J+1])/2)/((y[J+1] - y[J])/2);
    xiz = (zpt - (z[K] + z[K+1])/2)/((z[K+1] - z[K])/2);
        
    return;
}



///////////////////////////////////////////////////////////////////////
// natcoord
///////////////////////////////////////////////////////////////////////

void Grid::natcoord(double &xix, double &xiy, double &xiz, double xpt, double ypt, double zpt, int elem_num)
{
    vect3 F;
    mat3 Jmat;
    vect3 E;
    vect3 dE;
    vect3 newE;
    
    double err = 1;
    double tol = 1e-5;
    
    Elem elem = ebin[elem_num];
    
    double dN1[3] = {0};
    double dN2[3] = {0};
    double dN3[3] = {0};
    double dN4[3] = {0};
    double dN5[3] = {0};
    double dN6[3] = {0};
    double dN7[3] = {0};
    double dN8[3] = {0};
    double shapeF[8] = {0};
       
    int iter = 0;
	int max_iter = 6;
	
	while ((err > tol) && (iter < max_iter)){
        xix = E.x;
        xiy = E.y;
        xiz = E.z;
        
        shapefunctions(shapeF, xix, xiy, xiz);    
        shapefun_d1(dN1, xix, xiy, xiz, 1);
        shapefun_d1(dN2, xix, xiy, xiz, 2);
        shapefun_d1(dN3, xix, xiy, xiz, 3);
        shapefun_d1(dN4, xix, xiy, xiz, 4);
        shapefun_d1(dN5, xix, xiy, xiz, 5);
        shapefun_d1(dN6, xix, xiy, xiz, 6);
        shapefun_d1(dN7, xix, xiy, xiz, 7);
        shapefun_d1(dN8, xix, xiy, xiz, 8);
                
        F.x = xpt - (shapeF[0]*(*elem.n1).x + shapeF[1]*(*elem.n2).x + shapeF[2]*(*elem.n3).x + shapeF[3]*(*elem.n4).x + shapeF[4]*(*elem.n5).x + shapeF[5]*(*elem.n6).x + shapeF[6]*(*elem.n7).x + shapeF[7]*(*elem.n8).x);          
        F.y = ypt - (shapeF[0]*(*elem.n1).y + shapeF[1]*(*elem.n2).y + shapeF[2]*(*elem.n3).y + shapeF[3]*(*elem.n4).y + shapeF[4]*(*elem.n5).y + shapeF[5]*(*elem.n6).y + shapeF[6]*(*elem.n7).y + shapeF[7]*(*elem.n8).y); 
        F.z = zpt - (shapeF[0]*(*elem.n1).z + shapeF[1]*(*elem.n2).z + shapeF[2]*(*elem.n3).z + shapeF[3]*(*elem.n4).z + shapeF[4]*(*elem.n5).z + shapeF[5]*(*elem.n6).z + shapeF[6]*(*elem.n7).z + shapeF[7]*(*elem.n8).z); 
               
        Jmat.M11 = -(dN1[0]*(*elem.n1).x + dN2[0]*(*elem.n2).x + dN3[0]*(*elem.n3).x + dN4[0]*(*elem.n4).x + dN5[0]*(*elem.n5).x + dN6[0]*(*elem.n6).x + dN7[0]*(*elem.n7).x + dN8[0]*(*elem.n8).x);
        Jmat.M12 = -(dN1[1]*(*elem.n1).x + dN2[1]*(*elem.n2).x + dN3[1]*(*elem.n3).x + dN4[1]*(*elem.n4).x + dN5[1]*(*elem.n5).x + dN6[1]*(*elem.n6).x + dN7[1]*(*elem.n7).x + dN8[1]*(*elem.n8).x);
        Jmat.M13 = -(dN1[2]*(*elem.n1).x + dN2[2]*(*elem.n2).x + dN3[2]*(*elem.n3).x + dN4[2]*(*elem.n4).x + dN5[2]*(*elem.n5).x + dN6[2]*(*elem.n6).x + dN7[2]*(*elem.n7).x + dN8[2]*(*elem.n8).x);
        
        Jmat.M21 = -(dN1[0]*(*elem.n1).y + dN2[0]*(*elem.n2).y + dN3[0]*(*elem.n3).y + dN4[0]*(*elem.n4).y + dN5[0]*(*elem.n5).y + dN6[0]*(*elem.n6).y + dN7[0]*(*elem.n7).y + dN8[0]*(*elem.n8).y);
        Jmat.M22 = -(dN1[1]*(*elem.n1).y + dN2[1]*(*elem.n2).y + dN3[1]*(*elem.n3).y + dN4[1]*(*elem.n4).y + dN5[1]*(*elem.n5).y + dN6[1]*(*elem.n6).y + dN7[1]*(*elem.n7).y + dN8[1]*(*elem.n8).y);
        Jmat.M23 = -(dN1[2]*(*elem.n1).y + dN2[2]*(*elem.n2).y + dN3[2]*(*elem.n3).y + dN4[2]*(*elem.n4).y + dN5[2]*(*elem.n5).y + dN6[2]*(*elem.n6).y + dN7[2]*(*elem.n7).y + dN8[2]*(*elem.n8).y);
        
        Jmat.M31 = -(dN1[0]*(*elem.n1).z + dN2[0]*(*elem.n2).z + dN3[0]*(*elem.n3).z + dN4[0]*(*elem.n4).z + dN5[0]*(*elem.n5).z + dN6[0]*(*elem.n6).z + dN7[0]*(*elem.n7).z + dN8[0]*(*elem.n8).z);
        Jmat.M32 = -(dN1[1]*(*elem.n1).z + dN2[1]*(*elem.n2).z + dN3[1]*(*elem.n3).z + dN4[1]*(*elem.n4).z + dN5[1]*(*elem.n5).z + dN6[1]*(*elem.n6).z + dN7[1]*(*elem.n7).z + dN8[1]*(*elem.n8).z);
        Jmat.M33 = -(dN1[2]*(*elem.n1).z + dN2[2]*(*elem.n2).z + dN3[2]*(*elem.n3).z + dN4[2]*(*elem.n4).z + dN5[2]*(*elem.n5).z + dN6[2]*(*elem.n6).z + dN7[2]*(*elem.n7).z + dN8[2]*(*elem.n8).z);
        
        dE = -Jmat.invert()*F;
        newE = E + dE;
        
        err = dE.norm();
        E = newE;
		++iter;}
        
    xix = E.x;
    xiy = E.y;
    xiz = E.z;
    
    return;
}



///////////////////////////////////////////////////////////////////////
// shapefunctions
///////////////////////////////////////////////////////////////////////

// GRID.shapefunctions - Determines the shape function values for a given position in natural coordinates
//      Input:  - Position in natural coordinates (xix, xiy, xiz)
//              - Reference to array containing the 8 nodal shape function values (shapeF)
//
//      Output: - None (operates on references)

void Grid::shapefunctions(double (&shapeF)[8], double xix, double xiy, double xiz)
{
    shapeF[0] = ((1-xix)*(1-xiy)*(1-xiz))/8;                    // Shape function for node I, J, K
    shapeF[1] = ((1+xix)*(1-xiy)*(1-xiz))/8;                    // Shape function for node I+1, J, K
    shapeF[2] = ((1-xix)*(1+xiy)*(1-xiz))/8;                    // Shape function for node I, J+1, K
    shapeF[3] = ((1+xix)*(1+xiy)*(1-xiz))/8;                    // Shape function for node I+1, J+1, K
    shapeF[4] = ((1-xix)*(1-xiy)*(1+xiz))/8;                    // Shape function for node I, J, K+1
    shapeF[5] = ((1+xix)*(1-xiy)*(1+xiz))/8;                    // Shape function for node I+1, J, K+1
    shapeF[6] = ((1-xix)*(1+xiy)*(1+xiz))/8;                    // Shape function for node I, J+1, K+1
    shapeF[7] = ((1+xix)*(1+xiy)*(1+xiz))/8;                    // Shape function for node I+1, J+1, K+1
        
    return;
}



///////////////////////////////////////////////////////////////////////
// shapefun_d1
///////////////////////////////////////////////////////////////////////

// GRID.shapefun_d1 - For a given node (1 - 8), this function returns the value of the derivative of the shape function
//                    with respect to x, y, and z.
//      Input:  - Position in natural coordinates (xix, xiy, xiz)
//              - Reference to array containing the 3 first-order derivative shape function values (shapeF)
//              - Node number (1 - 8)
//
//      Output: - None (operates on references)

void Grid::shapefun_d1(double (&dshapeF)[3], const double xix, const double xiy, const double xiz, int node)
{
    if (node == 1)
    {
        dshapeF[0] = -(1 - xiy)*(1 - xiz)/8.;
        dshapeF[1] = -(1 - xix)*(1 - xiz)/8.;
        dshapeF[2] = -(1 - xix)*(1 - xiy)/8.;
    }
    
    if (node == 2)
    {
        dshapeF[0] = (1 - xiy)*(1 - xiz)/8.;
        dshapeF[1] = -(1 + xix)*(1 - xiz)/8.;
        dshapeF[2] = -(1 + xix)*(1 - xiy)/8.;
    }    
    
    if (node == 3)
    {
        dshapeF[0] = -(1 + xiy)*(1 - xiz)/8.;
        dshapeF[1] = (1 - xix)*(1 - xiz)/8.;
        dshapeF[2] = -(1 - xix)*(1 + xiy)/8.;
    }    
    
    if (node == 4)
    {
        dshapeF[0] = (1 + xiy)*(1 - xiz)/8.;
        dshapeF[1] = (1 + xix)*(1 - xiz)/8.;
        dshapeF[2] = -(1 + xix)*(1 + xiy)/8.;
    }    
    
    if (node == 5)
    {
        dshapeF[0] = -(1 - xiy)*(1 + xiz)/8.;
        dshapeF[1] = -(1 - xix)*(1 + xiz)/8.;
        dshapeF[2] = (1 - xix)*(1 - xiy)/8.;
    }    

    if (node == 6)
    {
        dshapeF[0] = (1 - xiy)*(1 + xiz)/8.;
        dshapeF[1] = -(1 + xix)*(1 + xiz)/8.;
        dshapeF[2] = (1 + xix)*(1 - xiy)/8.;
    }    
        
    if (node == 7)
    {
        dshapeF[0] = -(1 + xiy)*(1 + xiz)/8.;
        dshapeF[1] = (1 - xix)*(1 + xiz)/8.;
        dshapeF[2] = (1 - xix)*(1 + xiy)/8.;
    }   

    if (node == 8)
    {
        dshapeF[0] = (1 + xiy)*(1 + xiz)/8.;
        dshapeF[1] = (1 + xix)*(1 + xiz)/8.;
        dshapeF[2] = (1 + xix)*(1 + xiy)/8.;
    }   

    return;
}



///////////////////////////////////////////////////////////////////////
// nattoglobal
///////////////////////////////////////////////////////////////////////

void Grid::nattoglobal(double &xpt, double &ypt, double &zpt, double xix, double xiy, double xiz, int elem_num)
{
    Elem elem = ebin[elem_num];
    double shapeF[8];
    
    shapefunctions(shapeF, xix, xiy, xiz);

    xpt = shapeF[0]*(*elem.n1).x + shapeF[1]*(*elem.n2).x + shapeF[2]*(*elem.n3).x + shapeF[3]*(*elem.n4).x + shapeF[4]*(*elem.n5).x + shapeF[5]*(*elem.n6).x + shapeF[6]*(*elem.n7).x + shapeF[7]*(*elem.n8).x;
    ypt = shapeF[0]*(*elem.n1).y + shapeF[1]*(*elem.n2).y + shapeF[2]*(*elem.n3).y + shapeF[3]*(*elem.n4).y + shapeF[4]*(*elem.n5).y + shapeF[5]*(*elem.n6).y + shapeF[6]*(*elem.n7).y + shapeF[7]*(*elem.n8).y;
    zpt = shapeF[0]*(*elem.n1).z + shapeF[1]*(*elem.n2).z + shapeF[2]*(*elem.n3).z + shapeF[3]*(*elem.n4).z + shapeF[4]*(*elem.n5).z + shapeF[5]*(*elem.n6).z + shapeF[6]*(*elem.n7).z + shapeF[7]*(*elem.n8).z;
    
    return;
}



///////////////////////////////////////////////////////////////////////
// elem_find_neighbor
///////////////////////////////////////////////////////////////////////

int Grid::elem_find_neighbor(int elem_num,int neighbor_id)
{
	int elem_neighbor = -1;
	
	if (elem_num == -1)
		return elem_neighbor;

	Elem elem;
	Elem neighbor;
	
	elem = ebin[elem_num];

	// neighbor_ids:
	// 1 - Fronty (neighbor on face with normal -y)
	// 2 - Righty (neighbor on face with normal +x)
	// 3 - Backy  (neighbor on face with noraml +y)
	// 4 - Lefty  (neighbor on face with normal -x)
	// 5 - Toppy  (neighbor on face with normal +z)
	// 6 - Bottomy (neighbor on face with normal -z)
	
	for (int i = 0; i < Ne; i++)
	{
		if (i != elem_num)
		{
			neighbor = ebin[i];

			// Find Fronty
			if (neighbor_id == 1){
				if ((elem.n1 == neighbor.n3) && (elem.n2 == neighbor.n4) && (elem.n5 == neighbor.n7) && (elem.n6 == neighbor.n8)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Righty
			if (neighbor_id == 2){
				if ((elem.n2 == neighbor.n1) && (elem.n4 == neighbor.n3) && (elem.n6 == neighbor.n5) && (elem.n8 == neighbor.n7)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Backy
			if (neighbor_id == 3){
				if ((elem.n3 == neighbor.n1) && (elem.n4 == neighbor.n2) && (elem.n7 == neighbor.n5) && (elem.n8 == neighbor.n6)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Lefty
			if (neighbor_id == 4){
				if ((elem.n1 == neighbor.n2) && (elem.n3 == neighbor.n4) && (elem.n5 == neighbor.n6) && (elem.n7 == neighbor.n8)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Toppy
			if (neighbor_id == 5){
				if ((elem.n5 == neighbor.n1) && (elem.n6 == neighbor.n2) && (elem.n7 == neighbor.n3) && (elem.n8 == neighbor.n4)){
					elem_neighbor = i;
					return elem_neighbor;}}

			// Find Bottomy
			if (neighbor_id == 5){
				if ((elem.n1 == neighbor.n5) && (elem.n2 == neighbor.n6) && (elem.n3 == neighbor.n7) && (elem.n4 == neighbor.n8)){
					elem_neighbor = i;
					return elem_neighbor;}}
		}
	}

	return elem_neighbor;
}




double Grid::find_density_scale(double coll_den)
{
	double den_scale;
	double a = -0.016;
	double b = 5.1605;
	double c = 0.5112;
	den_scale = a + b*pow(E, -c*coll_den );                        // den_scale - Determine the density scaling factor using the function defined by a, b, c

	if (coll_den == 3.0)
		den_scale = 1.0;

	return den_scale;
}











///////////////////////////////////////////////////////////////////////
// find_density
///////////////////////////////////////////////////////////////////////

//void Grid::find_density(Segment &seg, Grid &grid, Data &data)
//{
//	// Current technique - march along segment distributing current segment fraction area to cell containing leading
//	// tip of fraction
//
//	// local variable declaration
//	double NN = 10; //divide segment into this many pieces
//	double V[3]; //displacement vector between segment tips
//	double position[3];
//	double natx,naty,natz;
//	double ss; //Step size
//	double unit_vec[3];
//	int i;
//	int Ilow, Jlow, Ihigh, Jhigh, Zlow, Zhigh;
//	double normV;
//	double area_vec[3], rho;
//	double S1, S2, S3, S4, S5, S6, S7, S8; //nodal shape functions
//	
//	if (seg.tip[1] == 1)
//	{
//		V[0] = seg.x[1] - seg.x[0];
//		V[1] = seg.y[1] - seg.y[0];
//		V[2] = seg.z[1] - seg.z[0];
//	}
//	else
//	{
//		V[0] = seg.x[0] - seg.x[1];
//		V[1] = seg.y[0] - seg.y[1];
//		V[2] = seg.z[0] - seg.z[1];
//	}
//
//	normV = vec_norm(V);
//
//	if (normV <= 0)
//		return;
//	
//
//	ss = normV/5;
//	unit_vec[0] = V[0]/normV;
//	unit_vec[1] = V[1]/normV;
//	unit_vec[2] = V[2]/normV;
//	if (seg.tip[1] == 1)
//	{
//		position[0] = seg.x[0];
//		position[1] = seg.y[0];
//		position[2] = seg.z[0];
//	}
//	else
//	{
//		position[0] =seg.x[1];
//		position[1] = seg.y[1];
//		position[2] = seg.z[1];
//	}
//
//	if (position[0] < grid.xrange[0])
//		position[0] = grid.xrange[0];
//	else if (position[0] > grid.xrange[1])
//		position[0] = grid.xrange[1];
//	if (position[1] < grid.yrange[0])
//		position[1] = grid.yrange[0];
//	else if (position[1] > grid.yrange[1])
//		position[1] = grid.yrange[1];
//	if (position[2] < grid.zrange[0])
//		position[2] = grid.zrange[0];
//	else if (position[2] > grid.zrange[1])
//		position[2] = grid.zrange[1];
//
//	for (i=0; i<NN; ++i)
//	{
//		//Increment current tip position
//		position[0] = position[0] + ss*unit_vec[0];
//		position[1] = position[1] + ss*unit_vec[1];
//		position[2] = position[2] + ss*unit_vec[2];
//		
//		//Determine what cell fraction tip is in
//		Ilow = floor((position[0]-grid.xrange[0])/data.dx);
//		Jlow = floor((position[1]-grid.yrange[0])/data.dy);
//		Zlow = floor((position[2]-grid.zrange[0])/data.dz);
//		Ihigh = Ilow+1;
//		Jhigh = Jlow+1;
//		Zhigh = Zlow+1;
//
//		
//		// Border Control
//
//		if (Ilow > grid.num_nodes[0]-2)
//		{
//			Ilow = grid.num_nodes[0]-2;
//			Ihigh = Ilow+1;
//			return;
//		}
//		else if (Ilow<0)
//		{
//			Ilow = 0;
//			Ihigh = 1;
//			return;
//		}
//
//		if (Jlow > grid.num_nodes[1]-2)
//		{
//			Jlow = grid.num_nodes[1]-2;
//			Jhigh = Jlow+1;
//			return;
//		}
//		else if (Jlow<0)
//		{
//			Jlow = 0;
//			Jhigh = 1;
//			return;
//		}
//		
//		if (Zlow > grid.num_nodes[2]-2)
//		{
//			Zlow = grid.num_nodes[2]-2;
//			Zhigh = Zlow+1;
//			return;
//		}
//		else if (Zlow<0)
//		{
//			Zlow = 0;
//			Zhigh = 1;
//			return;
//		} 
//		//Get displacement vector for current segment fraction
//		area_vec[0] = position[0]-(position[0]-ss*unit_vec[0]);
//		area_vec[1] = position[1]-(position[1]-ss*unit_vec[1]);
//		area_vec[2] = position[2]-(position[2]-ss*unit_vec[2]);
//		
//		//Determine density contribution of fraction to cell
//		rho = pi*vec_norm(area_vec)*data.vessel_width*data.vessel_width/4/(data.dx*data.dy*data.dz);
//
//		double xix, xiy, xiz;
//        double shapeF[8];
//		
//		// Convert to natural coordinates -1 <= Xi <= +1
//        grid.natcoordinates(xix, xiy, xiz, position[0], position[1], position[2], Ilow, Jlow, Zlow);        
//        
//        // Obtain shape function weights
//        grid.shapefunctions(shapeF, xix, xiy, xiz);
//		
//		
//		grid.density[Ilow][Jlow][Zlow] = grid.density[Ilow][Jlow][Zlow] + rho*shapeF[0];
//		grid.density[Ilow][Jhigh][Zlow] = grid.density[Ilow][Jhigh][Zlow] + rho*shapeF[1];
//		grid.density[Ihigh][Jlow][Zlow] = grid.density[Ihigh][Jlow][Zlow] + rho*shapeF[2];
//		grid.density[Ihigh][Jhigh][Zlow] = grid.density[Ihigh][Jhigh][Zlow] + rho*shapeF[3];
//		grid.density[Ilow][Jlow][Zhigh] = grid.density[Ilow][Jlow][Zhigh] + rho*shapeF[4];
//		grid.density[Ilow][Jhigh][Zhigh] = grid.density[Ilow][Jhigh][Zhigh] + rho*shapeF[5];
//		grid.density[Ihigh][Jlow][Zhigh] = grid.density[Ihigh][Jlow][Zhigh] + rho*shapeF[6];
//		grid.density[Ihigh][Jhigh][Zhigh] = grid.density[Ihigh][Jhigh][Zhigh] + rho*shapeF[7];
//	}
//	return;
//}

///////////////////////////////////////////////////////////////////////
// smooth_density
///////////////////////////////////////////////////////////////////////

//void Grid::smooth_density(double alpha[3][3][3], double density[xnodes][ynodes][znodes], int num_nodes[3])
//{
	 //density smoothing
		//int m,n,p;
		//int i,j,k;
		//double sum;
		//int indexi, indexj, indexk;
		//for (m=0; m< num_nodes[0]; ++m)
		//{
		//	for (n=0; n< num_nodes[1]; ++n)
		//	{
		//		for (p=0; p< num_nodes[2]; ++p)
		//		{
		//		sum = 0;
		//		for (i=0; i<3; ++i)
		//		{
		//			for (j=0; j<3; ++j)
		//			{
		//				for (k=0; k<3; ++k)
		//				{
		//					if (m< num_nodes[0]-1 && m>0 && n< num_nodes[1]-1 && n>0 && p< num_nodes[2]-1 && p>0)
		//						sum +=  density[m-1+i][n-1+j][p-1+k]*alpha[i][j][k];
		//					else if (m==0 && n< num_nodes[1] && n>0 && p< num_nodes[2] && p>0)
		//					{
		//						if (m-1+i == -1)
		//							indexi =  num_nodes[0]-1;
		//						else
		//							indexi = m-1+i;
		//						sum +=  density[indexi][n-1+j][p-1+k]*alpha[i][j][k];
		//					}
		//					else if (m==  num_nodes[0]-1 && n< num_nodes[1]-1 && n>0 && p< num_nodes[2]-1 && p>0)
		//					{
		//						if (m-1+i ==  num_nodes[0])
		//							indexi = 0;
		//						else
		//							indexi = m-1+i;
		//						sum +=  density[indexi][n-1+j][p-1+k]*alpha[i][j][k];
		//					}
		//					else if (n==0 && m< num_nodes[0]-1 && m>0 && p< num_nodes[2]-1 && p>0)
		//					{
		//						if (n-1+j == -1)
		//							indexj =  num_nodes[1]-1;
		//						else
		//							indexj = n-1+j;
		//						sum +=  density[m-1+i][indexj][p-1+k]*alpha[i][j][k];
		//					}
		//					else if (n== num_nodes[1]-1 && m< num_nodes[0]-1 && m>0 && p< num_nodes[2]-1 && p>0)
		//					{
		//						if (n-1+j ==  num_nodes[1])
		//							indexj = 0;
		//						else
		//							indexj = n-1+j;
		//						sum +=  density[m-1+i][indexj][p-1+k]*alpha[i][j][k];
		//					}
		//					else if (p==0 && m< num_nodes[0]-1 && m>0 && n< num_nodes[1]-1 && n>0)
		//					{
		//						if (p-1+k == -1)
		//							indexk =  num_nodes[2]-1;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[m-1+1][n-1+j][indexk]*alpha[i][j][k];
		//					}
		//					else if (p== num_nodes[2]-1 && m< num_nodes[0]-1 && m>0 && n< num_nodes[1]-1 && n>0)
		//					{
		//						if (p-1+k ==  num_nodes[2])
		//							indexk = 0;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[m-1+1][n-1+j][indexk]*alpha[i][j][k];
		//					}
		//					else if (m==0 && n==0 && p< num_nodes[2]-1 && p>0)
		//					{
		//						if (m-1+i == -1)
		//							indexi =  num_nodes[0]-1;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j == -1)
		//							indexj =  num_nodes[1]-1;
		//						else
		//							indexj = n-1+j;
		//						sum +=  density[indexi][indexj][p-1+k]*alpha[i][j][k];
		//					}
		//					else if (m== num_nodes[0]-1 && n== num_nodes[1]-1 && p< num_nodes[2]-1 && p>0)
		//					{
		//						if (m-1+i ==  num_nodes[0])
		//							indexi = 0;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j ==  num_nodes[1])
		//							indexj = 0;
		//						else
		//							indexj = n-1+j;
		//						sum +=  density[indexi][indexj][p-1+k]*alpha[i][j][k];
		//					}
		//					else if (m==0 && n< num_nodes[1]-1 && n>0 && p==0)
		//					{
		//						if (m-1+i == -1)
		//							indexi =  num_nodes[0]-1;
		//						else
		//							indexi = m-1+i;
		//						if (p-1+k == -1)
		//							indexk =  num_nodes[2]-1;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][n-1+j][indexk]*alpha[i][j][k];
		//					}
		//					else if (m== num_nodes[0]-1 && n< num_nodes[1]-1 && n>0 && p== num_nodes[2]-1)
		//					{
		//						if (m-1+i ==  num_nodes[0])
		//							indexi = 0;
		//						else
		//							indexi = m-1+i;
		//						if (p-1+k ==  num_nodes[2])
		//							indexk = 0;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][n-1+j][indexk]*alpha[i][j][k];
		//					}

		//					else if (m< num_nodes[0]-1 && m>0 && n==0 && p==0)
		//					{
		//						if (n-1+j == -1)
		//							indexj =  num_nodes[1]-1;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k == -1)
		//							indexk =  num_nodes[2]-1;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[m-1+i][indexj][indexk]*alpha[i][j][k];
		//					}
		//					else if (m< num_nodes[0]-1 && m>0 && n== num_nodes[1]-1 && p== num_nodes[2]-1)
		//					{
		//						if (n-1+j ==  num_nodes[1])
		//							indexj = 0;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k ==  num_nodes[2])
		//							indexk = 0;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[m-1+i][indexj][indexk]*alpha[i][j][k];
		//					}
		//					else if (m==0 && n==0 && p==0)
		//					{
		//						if (m-1+i == -1)
		//							indexi =  num_nodes[0]-1;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j == -1)
		//							indexj =  num_nodes[1]-1;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k == -1)
		//							indexk =  num_nodes[2]-1;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][indexj][indexk]*alpha[i][j][k];
		//					}
		//					else if (m== num_nodes[0]-1 && n== num_nodes[1]-1 && p== num_nodes[2]-1)
		//					{
		//						if (m-1+i ==  num_nodes[0])
		//							indexi = 0;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j ==  num_nodes[1])
		//							indexj = 0;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k ==  num_nodes[2])
		//							indexk = 0;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][indexj][indexk]*alpha[i][j][k];
		//					}
		//					else if (m==0 && n== num_nodes[1]-1 && p== num_nodes[2]-1)
		//					{
		//						if (m-1+i == -1)
		//							indexi =  num_nodes[0]-1;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j ==  num_nodes[1])
		//							indexj = 0;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k ==  num_nodes[2])
		//							indexk = 0;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][indexj][indexk]*alpha[i][j][k];
		//					}
		//					
		//					else if (m==0 && n==0 && p== num_nodes[2]-1)
		//					{
		//						if (m-1+i == -1)
		//							indexi =  num_nodes[0]-1;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j == -1)
		//							indexj =  num_nodes[1]-1;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k ==  num_nodes[2])
		//							indexk = 0;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][indexj][indexk]*alpha[i][j][k];
		//					}
		//					else if (m==0 && n== num_nodes[1]-1 && p==0)
		//					{
		//						if (m-1+i == -1)
		//							indexi =  num_nodes[0]-1;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j ==  num_nodes[1])
		//							indexj = 0;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k == -1)
		//							indexk =  num_nodes[2]-1;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][indexj][indexk]*alpha[i][j][k];
		//					}
		//					else if (m== num_nodes[0]-1 && n==0 && p== num_nodes[2]-1)
		//					{
		//						if (m-1+i ==  num_nodes[0])
		//							indexi = 0;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j == -1)
		//							indexj =  num_nodes[1]-1;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k ==  num_nodes[2])
		//							indexk = 0;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][indexj][indexk]*alpha[i][j][k];
		//					}
		//					else if (m== num_nodes[0]-1 && n== num_nodes[1]-1 && p==0)
		//					{
		//						if (m-1+i ==  num_nodes[0])
		//							indexi = 0;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j ==  num_nodes[1])
		//							indexj = 0;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k == -1)
		//							indexk =  num_nodes[2]-1;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][indexj][indexk]*alpha[i][j][k];
		//					}
		//					else if (m== num_nodes[0]-1 && n==0 && p==0)
		//					{
		//						if (m-1+i ==  num_nodes[0])
		//							indexi = 0;
		//						else
		//							indexi = m-1+i;
		//						if (n-1+j == -1)
		//							indexj =  num_nodes[1]-1;
		//						else
		//							indexj = n-1+j;
		//						if (p-1+k == -1)
		//							indexk =  num_nodes[2]-1;
		//						else
		//							indexk = p-1+k;
		//						sum +=  density[indexi][indexj][indexk]*alpha[i][j][k];
		//					}

		//				}
		//			}
		//		}
		//		}
		//		 smoothed_dens[m][n][p] = sum;
		//	}
		//}
		//for (m=0; m< num_nodes[0]; ++m)
		//{
		//	for (n=0; n< num_nodes[1]; ++n)
		//	{
		//		for (p=0; p< num_nodes[2]; ++p)
		//		{
		//			 density[m][n][p] =  smoothed_dens[m][n][p];
		//		}
		//
		//	}
		//}
//}//end density smoothing


///////////////////////////////////////////////////////////////////////
// Gradient
///////////////////////////////////////////////////////////////////////
//void Grid::Gradient(Grid &grid, Data &data)
//{
	//int i,j,k;
	//for (i=0;i<grid.num_nodes[0];++i)
	//{
	//	for (j=0;j<grid.num_nodes[1];++j)
	//	{
	//		for (k=0;k<grid.num_nodes[2];++k)
	//		{
	//			if (i == 0) //use forward difference
	//				grid.drho_dx[i][j][k] = (grid.density[i+1][j][k]-grid.density[i][j][k])/data.dx;
	//			else if (i == grid.num_nodes[0] - 1) //use backward difference
	//				grid.drho_dx[i][j][k] = (grid.density[i][j][k]-grid.density[i-1][j][k])/data.dx;
	//			else //use central difference
	//				grid.drho_dx[i][j][k] = (grid.density[i+1][j][k]-grid.density[i-1][j][k])/(2.*data.dx);

	//			if (j == 0) //use forward difference
	//				grid.drho_dy[i][j][k] = (grid.density[i][j+1][k]-grid.density[i][j][k])/data.dy;
	//			else if (j == grid.num_nodes[1] - 1) //use backward difference
	//				grid.drho_dy[i][j][k] = (grid.density[i][j][k]-grid.density[i][j-1][k])/data.dy;
	//			else //use central difference
	//				grid.drho_dy[i][j][k] = (grid.density[i][j+1][k]-grid.density[i][j-1][k])/(2.*data.dy);

	//			if (k == 0) //use forward difference
	//				grid.drho_dz[i][j][k] = (grid.density[i][j][k+1]-grid.density[i][j][k])/data.dz;
	//			else if (k == grid.num_nodes[2] - 1) //use backward difference
	//				grid.drho_dz[i][j][k] = (grid.density[i][j][k]-grid.density[i][j][k-1])/data.dz;
	//			else //use central difference
	//				grid.drho_dz[i][j][k] = (grid.density[i][j][k+1]-grid.density[i][j][k-1])/(2.*data.dz);
	//		}
	//	}
	//}
//	return;
//}

///////////////////////////////////////////////////////////////////////
// Clear density
///////////////////////////////////////////////////////////////////////
//void Grid::clear_density(Grid &grid)
//{
	/*int i,j,k;
	for (i=0; i<xnodes; ++i)
	{
		for (j=0; j<ynodes; ++j)
		{
			for (k=0; k<znodes; ++k)
			{
				grid.density[i][j][k] = 0;
			}
		}
	}*/
//	return;
//}