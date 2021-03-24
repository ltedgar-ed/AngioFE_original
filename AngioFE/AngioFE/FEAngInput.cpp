///////////////////////////////////////////////////////////////////////
// FEAngInput.cpp
///////////////////////////////////////////////////////////////////////



#include "StdAfx.h"
#include "FEAngInput.h"
#include "Angio.h"
#include "angio3d.h"
#include <iostream>
#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEAngInput::FEAngInput() : Input()
{   
	// Set parameters to default values
   	isproutf = 1.0;					// Sprout force magnitude	
	itip_range = 250.;				// Sprout force range
	ispfactor = 0.;					// Sprout force directional factor
	
	isub_cyc = 2;					// Number of subgrowth cycles
	igelbc = 'u';					// Gel boundary conditions ('u' unconstrained)
	istifffact = 1.0;				// Displacement stiffness factor
	isprout_verify = 0;				// Sprout verification problem flag
	isp_sphere = 0;					// Switch between local directional (0), local isotropic (1), and global isotropic (2) sprout froce representations

	iSx = 0.;						// Location of the x-symmetry plane
	iSy = 0.;						// Location of the y-symmetry plane
	iSz = 0.;						// Location of the z-symmetry plane

	irseed = 0;						// Input random seed number
			
	ibranch = 1;					// Branching flag
	ianast = 1;						// Anastomosis flag

	izfibflat = 0;

	icirc = 0;

	iNe = 0;
}

FEAngInput::~FEAngInput()
{

}


