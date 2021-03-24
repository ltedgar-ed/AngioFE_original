///////////////////////////////////////////////////////////////////////
// Filein.cpp
///////////////////////////////////////////////////////////////////////



#include "StdAfx.h"
#include "Input.h"
#include "Filein.h"
#include "angio3d.h"
#include "Elem.h"
#include <iostream>
#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Filein::Filein(Input &input)
{
    in_file.open(filename);

    const int buff_max = 100;                   // Maximum size of the buffer
    char buffer[buff_max];                      // Create the buffer
    char first_char;
    char str_type = 'n';                        // Set default string type to 'no line'
    
    while (! in_file.eof() )
    {
        str_type = 'n';                         // Reset the default string type to 'no line'
        in_file.getline (buffer,255,'\n');      // Read in the next line of the input file into the buffer  
                
        first_char = buffer[0];                 // Determine the first character of the line in the buffer
       
        if (first_char == 37)                   // If first character is a '%' symbol...
             str_type = 'c';                        // ... then set string type to 'Comment'
             else if (first_char == 62)         // If first character is a '>' symbol...
             str_type = 'p';                        // ... then set string type to 'Param'
                      
        if (str_type == 'p')
            read_param(input, buffer);                 // If the line describes a parameter, then read in that parameter
    }

    if (input.igrid_in == 'y'){
        read_nodes(input);
        read_econn(input);
        read_eBC(input);
        //read_collfib(input);
        //read_disp(input);
        }    

	//read_coll_den(input);
}



Filein::~Filein()
{

}



//////////////////////////////////////////////////////////////////////
// Member Functions
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// read_param
///////////////////////////////////////////////////////////////////////

void Filein::read_param(Input &input, char (&buffer)[100])
{
    char pname[20] = {0};                       // Create string for the parameter name
    
    sscanf(buffer,"%*s %s %*f",&pname);         // Scan the buffer for the name of the parameter
    
    set_param(input,buffer,pname);                    // Set the parameter based on the name
    
    return;
}



///////////////////////////////////////////////////////////////////////
// set_param
///////////////////////////////////////////////////////////////////////
    
void Filein::set_param(Input &input, char (&buffer)[100], char (&pname)[20])
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
    
	// Read in weights for determing the direction of growth
    if (!strcmp(pname,"gweights")){
        sscanf(buffer,"%*s %*s %lf %lf",&input.iweight1,&input.iweight4);
        return;}

    return;
}



///////////////////////////////////////////////////////////////////////
// read_nodes
///////////////////////////////////////////////////////////////////////

void Filein::read_nodes(Input &input)
{
    const int buff_max = 100;                   // Maximum size of the buffer
    char buffer[buff_max];                      // Create the buffer
            
    gr_file.open(gridname);
      
    int count = 0;
    int node_id = 0; float xpt = 0.; float ypt = 0.; float zpt = 0.;
    int readin = 0;
    
    while (! gr_file.eof() )
    {
        gr_file.getline (buffer,255,'\n');      
        readin = sscanf(buffer,"%i %f %f %f", &node_id, &xpt, &ypt, &zpt); 
        
        if (readin == 4){
            Node node;
            node.id = node_id;
            node.x = (double)xpt;
            node.y = (double)ypt;
            node.z = (double)zpt;
            input.inodes.push_back(node);
            count++;}
    }
    
    input.iNn = count;
        
    return;
}



///////////////////////////////////////////////////////////////////////
// read_coll_den
///////////////////////////////////////////////////////////////////////

void Filein::read_coll_den(Input &input)
{
    const int buff_max = 100;                   // Maximum size of the buffer
    char buffer[buff_max];                      // Create the buffer
            
    gr_file.open(coll_den_name);
      
    int count = 0;
    float coll_den;
    int readin = 0;
    
    while (! gr_file.eof() )
    {
        gr_file.getline (buffer,255,'\n');      
        readin = sscanf(buffer,"%f", &coll_den); 
        
        if (readin == 1){
            input.icoll_den.push_back(coll_den);
            count++;}
    }
    
    return;
}


///////////////////////////////////////////////////////////////////////
// read_econn
///////////////////////////////////////////////////////////////////////

void Filein::read_econn(Input &input)
{
    const int buff_max = 100;                   // Maximum size of the buffer
    char buffer[buff_max];                      // Create the buffer
            
    ec_file.open(econname);
      
    int count = 0;
    int elem_num = 0; int n1 = 0; int n2 = 0; int n3 = 0; int n4 = 0; int n5 = 0; int n6 = 0; int n7 = 0; int n8 = 0;
    int readin = 0;
    
    while (! ec_file.eof() )
    {
        ec_file.getline (buffer,255,'\n');      
        readin = sscanf(buffer,"%i %i %i %i %i %i %i %i %i", &elem_num, &n1, &n2, &n3, &n4, &n5, &n6, &n7, &n8); 
        
        if (readin == 9){
            vector<int> elem_conn;
            elem_conn.resize(9);
            elem_conn[0] = elem_num;
            elem_conn[1] = n1;
            elem_conn[2] = n2;
            elem_conn[3] = n3;
            elem_conn[4] = n4;
            elem_conn[5] = n5;
            elem_conn[6] = n6;
            elem_conn[7] = n7;
            elem_conn[8] = n8;
            
            input.ieconn.push_back(elem_conn);
            count++;}
    }
    
    input.iNe = count;
        
    return;
}



///////////////////////////////////////////////////////////////////////
// read_eBC
///////////////////////////////////////////////////////////////////////

void Filein::read_eBC(Input &input)
{
    const int buff_max = 100;                   // Maximum size of the buffer
    char buffer[buff_max];                      // Create the buffer
            
    bc_file.open(eBCname);
    
    int count = 0;
    int elem_num = 0; int BC1 = 0; int BC2 = 0; int BC3 = 0; int BC4 = 0; int BC5 = 0; int BC6 = 0;
    
     
    while (! bc_file.eof() )
    {
        bc_file.getline (buffer,255,'\n');      
        sscanf(buffer,"%i %i %i %i %i %i %i", &elem_num, &BC1, &BC2, &BC3, &BC4, &BC5, &BC6); 
        
        vector<int> eBC_violate;
        eBC_violate.resize(7);
        eBC_violate[0] = elem_num;
        eBC_violate[1] = BC1;
        eBC_violate[2] = BC2;
        eBC_violate[3] = BC3;
        eBC_violate[4] = BC4;
        eBC_violate[5] = BC5;
        eBC_violate[6] = BC6;
                
        input.ieBC.push_back(eBC_violate);
        count++;
    }
    
    input.iNBC = count - 1;
        
    return;
}



///////////////////////////////////////////////////////////////////////
// read_collfib
///////////////////////////////////////////////////////////////////////

void Filein::read_collfib(Input &input)
{
    const int buff_max = 100;                   // Maximum size of the buffer
    char buffer[buff_max];                      // Create the buffer
            
    fib_file.open(fibname);
    
    int count = 0;
    int nnum = 0; double itheta = 0.; double ieta = 0.;
    
         
    while (! fib_file.eof() )
    {
        fib_file.getline (buffer,255,'\n');      
        sscanf(buffer,"%i %lf %lf", &nnum, &itheta, &ieta); 
        input.inodes[nnum].theta = itheta;
        input.inodes[nnum].eta = ieta;

        count++;
    }
            
    return;
}


///////////////////////////////////////////////////////////////////////
// read_nodes
///////////////////////////////////////////////////////////////////////

void Filein::read_disp(Input &input)
{
    const int buff_max = 100;                   // Maximum size of the buffer
    char buffer[buff_max];                      // Create the buffer
            
    u_file.open(uname);
      
    int count = 0;
    int nnum = 0; float ux = 0.; float uy = 0.; float uz = 0.;
    int readin = 0;
    
    while (! u_file.eof() )
    {
        u_file.getline (buffer,255,'\n');      
        readin = sscanf(buffer,"%i %f %f %f", &nnum, &ux, &uy, &uz); 
        
        input.inodes[nnum].u.x = ux;
        input.inodes[nnum].u.y = uy;
        input.inodes[nnum].u.z = uz;

    }
    
            
    return;
}