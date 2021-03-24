///////////////////////////////////////////////////////////////////////
// Fileout.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Fileout.h"
#include "Segment.h"
#include "Culture.h"
#include "Data.h"
#include "angio3d.h"
#include "Grid.h"
#include "Elem.h"
#include "Angio.h"
#include <iostream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Fileout::Fileout()
{
    logstream.open("out_log.ang");
    stream3 = fopen("tracking.ang","wt");   // tracking.ang: time step, model time, total length in culture, number of branches in culture
}

Fileout::~Fileout()
{
    logstream.close();
}

///////////////////////////////////////////////////////////////////////
// Member Functions
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// timestart
///////////////////////////////////////////////////////////////////////

void Fileout::timestart()
{
    time(&start);                                               // Start the timer
    return;
}



///////////////////////////////////////////////////////////////////////
// writeTracking
///////////////////////////////////////////////////////////////////////

void Fileout::writeTracking(Data &data)
{
    fprintf(stream3,"%-12.7f %-12.7f %-12.7f %-5i\n",data.dt,data.t,data.total_length,data.num_branches);   // Write to tracking.ang
    
    return;
}



///////////////////////////////////////////////////////////////////////
// closeTracking
///////////////////////////////////////////////////////////////////////

void Fileout::closeTracking()
{
    fclose(stream3);                                                        // Close stream to 'tracking.ang' (stream3) 
    
    return;
}



///////////////////////////////////////////////////////////////////////
// printStatus
///////////////////////////////////////////////////////////////////////

void Fileout::printStatus(Data &data)
{
    cout << endl << "Time: " << data.t << endl;                             // Print out current time to user
	//cout << "dt: " << data.dt << endl;
    cout << "Segments: " << data.nsegs << endl;                             // Print out current number of segments to user
	cout << "Total Length: " << data.total_length << endl;                  // Print out the current total length to user
	cout << "Branch Points: " << data.num_branches << endl;                 // Print out the current number of branches to user
	cout << "Anastomoses: " << data.num_anastom << endl << endl;            // Print out the current number of anastomoses to user
    
    logstream << endl << "Time: " << data.t << endl;                        // Print out current time to log file
	//logstream << "dt: " << data.dt << endl;
    logstream << "Segments: " << data.nsegs << endl;                        // Print out current number of segments to log file
	logstream << "Total Length: " << data.total_length << endl;             // Print out the current total length to log file
	logstream << "Branch Points: " << data.num_branches << endl;            // Print out the current number of branches to log file
	logstream << "Anastomoses: " << data.num_anastom << endl << endl;       // Print out the current number of anastomoses to log file
        
    return;
}



///////////////////////////////////////////////////////////////////////
// dataout
///////////////////////////////////////////////////////////////////////

void Fileout::dataout(Angio &angio, Data &data, Grid &grid)
{
    writeData(angio.frag);                                    // Create and write to 'data.ang'
	
    writeGrid(data, grid);                              // Create and write to 'grid.ang'
   
    writeNodes(data, grid);
    
    writeEconn(data, grid);
    
    writeBC(grid);
    
    //writeGrad(data, grid);                              // Create and write to 'grad.ang'
    
    //writeAngle(frag);                                   // Create and write to 'angle.ang'

    printtime();                                        // Display the run-time to the user

	writeSegConn(angio.frag);

    return;
}



///////////////////////////////////////////////////////////////////////
// writeData
///////////////////////////////////////////////////////////////////////

void Fileout::writeData(list<Segment> &frag)
{
	list<Segment>::iterator it;
	
	///// File output: 'data.ang' /////
	
	FILE *stream;                                                           // Open stream to 'data.ang' (stream)
	stream = fopen("out_data.ang","wt");                                        // data.ang: Store 3D coordinates of begining and end of each vessel segment
	                                                                        // as well as total length of the segment
	
	fprintf(stream,"%-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-5s %-5s\n","Time","X1","Y1","Z1","X2","Y2","Z2","Length","Vess","Label");  // Write column labels to data.ang
	
	for (it = frag.begin(); it != frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	{
		fprintf(stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i %-5.2i\n",it->TofBirth,it->x[0],it->y[0],it->z[0],it->x[1],it->y[1],it->z[1],it->length,it->seg_num,it->label);  // Write to data.ang
	    //fprintf(stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i %-5.2i\n",it->TofBirth,it->x[0],it->y[0],it->z[0],it->x[1],it->y[1],it->z[1],it->length,it->sprout,it->label);  // Write to data.ang
	}
	fclose(stream);                                                         // Close stream to 'data.ang' (stream)
	
	return;
}



///////////////////////////////////////////////////////////////////////
// writeGrid
///////////////////////////////////////////////////////////////////////

void Fileout::writeGrid(Data &data, Grid &grid)
{
    /// File output: 'grid.ang' /////
     
	FILE *stream2;                                                          // Open stream to 'grid.ang' (stream2)                                                                   
	stream2 = fopen("out_grid.ang","wt");                                       // grid.ang: Store 3D coordinates of grid nodes, collagen fiber orientation at the nodes,
	                                                                        // and vessel density at the nodes.
	int i; int j; int k;                                                                  
	
	for (i=0; i < grid.num_nodes[0]; ++i)                                 // Iterate through nodes in the x-direction (i)
		for (j=0; j < grid.num_nodes[1]; ++j)                             // Iterate through nodes in the y-direction (j)
			for (k=0; k < grid.num_nodes[2]; ++k)                         // Iterate through nodes in the z-direction (k)
				//fprintf(stream2,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f \n", grid.x[i], grid.y[j], grid.z[k], grid.theta[i][j][k], grid.eta[i][j][k]);   // Write to grid.ang
				fprintf(stream2,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n", grid.x[i], grid.y[j], grid.z[k], grid.coll_fib_x[i][j][k], grid.coll_fib_y[i][j][k], grid.coll_fib_z[i][j][k]);   // Write to grid.ang
	fclose(stream2);                                                        // Close stream to 'grid.ang' (stream2)
	
	return;
}



///////////////////////////////////////////////////////////////////////
// writeNodes
///////////////////////////////////////////////////////////////////////

void Fileout::writeNodes(Data &data, Grid &grid)
{
    /// File output: 'nodes.ang' /////
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_nodes.ang","wt");                                       
	Node node;
	
	for (int i = 0; i < grid.Nn; ++i){
	    node = grid.nodes[i];
	    fprintf(stream2, "%-5.2i %-12.7f %-12.7f %-12.7f\n", node.id, node.x, node.y, node.z);
	}  
	                                                                      	
	fclose(stream2);                                                        
	
	return;
}



///////////////////////////////////////////////////////////////////////
// writeEconn
///////////////////////////////////////////////////////////////////////

void Fileout::writeEconn(Data &data, Grid &grid)
{
    /// File output: 'econn.ang' /////
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_econn.ang","wt");                                       
	Elem elem;
	
	for (int i = 0; i < grid.Ne; ++i){
	    elem = grid.ebin[i];
	    fprintf(stream2, "%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i \n", elem.elem_num, (*elem.n1).id, (*elem.n2).id, (*elem.n3).id, (*elem.n4).id, (*elem.n5).id, (*elem.n6).id,  (*elem.n7).id,  (*elem.n8).id);
	}  
	                                                                      	
	fclose(stream2);                                                        
	
	return;
}


///////////////////////////////////////////////////////////////////////
// writeGrad
///////////////////////////////////////////////////////////////////////

void Fileout::writeGrad(Data &data, Grid &grid)
{
    /// File output: 'grad.ang' /////

    FILE *stream4;                                                          // Open stream to 'grad.ang' (stream4)
    stream4 = fopen("out_grad.ang","wt");                                       // grad.ang: Store 3D coordinates of grid nodes, collagen fiber orientation at the nodes,
	                                                                        // vessel density at the nodes, and value of the density gradient at the nodes.
    int i; int j; int k;
    
    for (i=0; i < grid.num_nodes[0]; ++i)                             
      for (j=0; j < grid.num_nodes[1]; ++j) 
        for (k=0; k < grid.num_nodes[2]; ++k)
            //fprintf(stream4,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n",i*data.dx,j*data.dy,k*data.dz,grid.density[i][j][k],grid.drho_dx[i][j][k],grid.drho_dy[i][j][k],grid.drho_dz[i][j][k]); // Write to grad.ang    
    
    fclose(stream4); 
    
    return;
}



///////////////////////////////////////////////////////////////////////
// writeAngle
///////////////////////////////////////////////////////////////////////

void Fileout::writeAngle(list<Segment> &frag)
{
 //   list<Segment>::iterator it;
 //   
 //   /// File output: 'angle.ang' /////
	//	
	//FILE *stream5;                                                          // Open stream to 'angle.ang' (stream5)
	//stream5 = fopen("out_angle.ang","wt");                                      // angle.ang: Store the vessel label, time of birth, sprout marker,
	//                                                                        // angle in the x-y plane, and angle off the x-y plane for each segment
	//                                                                        
	//for (it = frag.begin(); it != frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	//{
	//	fprintf(stream5,"%-5.2i %15.7f %-5.1i %-12.7f %-12.7f \n", it->label, it->TofBirth, it->sprout, it->phi1, it->phi2);    // Write to angle.ang
	//}
	//
	//fclose(stream5);                                                        // Close stream to 'angle.ang' (stream5)
	
	return;
}



///////////////////////////////////////////////////////////////////////
// writeBC
///////////////////////////////////////////////////////////////////////

void Fileout::writeBC(Grid &grid)
{
    /// File output: 'eBC.ang' /////
     
	FILE *stream2;                                                                                                                           
	stream2 = fopen("out_eBC.ang","wt");                                       
	Elem elem;
	
	int BC_violate[6] = {0};
	
	for (int i = 0; i < grid.Ne; ++i){
	    elem = grid.ebin[i];
	    
	    for (int j = 0; j < 6; ++j)
	        BC_violate[j] = 0;
	    
	    if ((elem.f1.BC == true) || (elem.f2.BC == true) || (elem.f3.BC == true) || (elem.f4.BC == true) || (elem.f5.BC == true) || (elem.f6.BC == true)){
	        if (elem.f1.BC == true)
	            BC_violate[0] = 1;
	        if (elem.f2.BC == true)
	            BC_violate[1] = 1;
	        if (elem.f3.BC == true)
	            BC_violate[2] = 1;
	        if (elem.f4.BC == true)
	            BC_violate[3] = 1;
	        if (elem.f5.BC == true)
	            BC_violate[4] = 1;
	        if (elem.f6.BC == true)
	            BC_violate[5] = 1;   
	        
	        fprintf(stream2, "%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i %-5.2i \n", elem.elem_num, BC_violate[0], BC_violate[1], BC_violate[2], BC_violate[3], BC_violate[4], BC_violate[5]);}
	}  
	                                                                      	
	fclose(stream2);                                                        
	
	return;

}


///////////////////////////////////////////////////////////////////////
// printtime
///////////////////////////////////////////////////////////////////////

void Fileout::printtime()
{
    time(&stop);                                                // Stop the timer
	
	t_seconds = (double) difftime(stop, start);                 // Calculate the simulation time in seconds
	
	cout << endl << "Simulation time: " << t_seconds << " seconds (" 
	    << floor(t_seconds/60) << " minutes)." << endl << endl;                // Show the user how long the simulation took (in seconds)
    
    logstream << endl << "Simulation time: " << t_seconds << " seconds (" << floor(t_seconds/60) << " minutes)." << endl << endl;  
	    
    return;
}

void Fileout::printrandseed(int randseed)
{
	logstream << endl << "Rand seed:" << randseed << endl << endl;
}


void Fileout::writeSegConn(list<Segment> &frag)
{
	list<Segment>::iterator it;
	
	///// File output: 'out_seg_conn.ang' /////
	
	FILE *stream;                                                           // Open stream to 'out_seg_conn.ang' (stream)
	stream = fopen("out_seg_conn.ang","wt");                                        
	
	fprintf(stream,"%-5s %-5s %-5s %-5s %-5s\n","SNum","Node1","     ","Node2","     ");  // Write column labels to data.ang

	for (it = frag.begin(); it != frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	{
		fprintf(stream,"%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i\n",it->seg_num,it->seg_conn[0][0],it->seg_conn[0][1],it->seg_conn[1][0],it->seg_conn[1][1]);  // Write to seg_conn.ang
	
	}
	
	fclose(stream);                                                         // Close stream to 'data.ang' (stream)
	
	return;
}