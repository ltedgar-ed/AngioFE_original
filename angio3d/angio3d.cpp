///////////////////////// Angiogenesis Growth Model /////////////////////////
////////////////////// Written by Lowell Taylor Edgar ///////////////////////
///////////////////////////// University of Utah ////////////////////////////
////////////////////////////// Copyright 2010 ///////////////////////////////
////////////// Original model by Scott Sibole and Jim Guilkey ///////////////

/////////////////////////////////////////////////////////////////////////////
// angio3d.cpp
/////////////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include <iostream>
#include "assert.h"

#include "angio3d.h"
#include "Angio.h"
#include "Filein.h"
#include "Fileout.h"
#include "Input.h"
#include "Grid.h"
#include "Data.h"

#include "mat3.h"
#include "vect3.h"

using namespace std;



/////////////////////////////// BEGIN MAIN /////////////////////////////////

int main(int nargs, char* argv[])                               // Main executable function
{
    time_t ranseed = time(NULL);
	srand((unsigned int)ranseed);
    //srand(1234567890);

////////////////////////////// DECLARE OBJECTS ///////////////////////////////
    
    Fileout fileout;                                            // Declare and initialize FILEOUT object 'fileout' 
		
	fileout.timestart();                                        // Start the timer	
		
	Input input;                                                // Declare and initialize INPUT object 'input'
	
	Filein filein(input);
    
    Angio angio(input);                                          // Declare and initialize ANGIO object 'angio'	
                                              
	
////////////////////// SEED INITIAL VESSEL FRAGMENTS ///////////////////////
	
	angio.seedFrag();											// Seed the initial fragments within the simulation domain
	
	//angio.find_active_tips();
	
	angio.initBranch();											// Handle branching within inital fragments
    
    angio.updateTotalLength();									// Update the total vascular length within the simulation   
    
    //fileout.writeTracking(data);                              // Write to 'tracking.ang'

/////////////////////////////// ANGIOGENESIS ///////////////////////////////

    while (angio.data.t < angio.data.maxt)                                  // While culture time is less than the maximum time
	{
		///// Update time /////
    
        angio.updateTime();                                 // Determine the current time step and update time
        angio.updateLength();        		

		///// Vessel Growth /////
	    
    	angio.Growth();                               // Elongate active segments found within a_frag
    	
		angio.kill_dead_segs();
		
		angio.updateTotalLength();                          // Update the total vascular length within the simulation
	    
		angio.find_active_tips();

		//angio.find_newborns();

		//fileout.writeTracking(data);                            // Write to 'tracking.ang'
	    
	    fileout.printStatus(angio.data);                              // Print the status of the simulation to the user
	    
	    
	    ///// End Time Step /////
	    
	}                                                           // End 'while' loop once max time is reached
	
	
//////////////////////////// WRITE OUTPUT FILES ////////////////////////////

	//angio.displace_vessels(grid);
	
	//angio.removeErrors();
	
	//frag = angio.returnFragList();                              // Create the complete list of Segments 'frag'
	
	//fileout.closeTracking();                                      // Close data stream to 'tracking.ang'
                    
    fileout.dataout(angio, angio.data, angio.grid);

    //cout << endl << ranseed << endl;
	
 //////////////////////////////// END MAIN ////////////////////////////////
	
	return 0;
}