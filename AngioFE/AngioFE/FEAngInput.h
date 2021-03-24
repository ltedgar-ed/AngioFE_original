///////////////////////////////////////////////////////////////////////
// FEAngInput.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The INPUT class reads in the input file and stores the prescribed parameter values used.
///////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <vector>

#include "Input.h"

class Node;

using namespace std;

class FEAngInput : public Input
{
public:
	FEAngInput();
	virtual ~FEAngInput();
			  
	// AngioFE params
	double isproutf;					// Sprout force magnitude
	double itip_range;					// Sprout force range
	double ispfactor;					// Sprout force directional factor
	
	int isub_cyc;						// Number of subgrowth cycles
	char igelbc;						// Boundary conditions for the gel (LAC, SAC, UNC)
	double istifffact;

	int isprout_verify;					// Flag for the sprout verification problem
    int isp_sphere;						// Flag for sprout force representations

	double iSx;							// Location of the symmetry plane along the x-direction
	double iSy;							// Location of the symmetry plane along the y-direction
	double iSz;							// Location of the symmetry plane along the z-direction

	int irseed;							// Seed number for the random number generator

	int ibranch;						// Flag to turn branching on/off
	int ianast;							// Flag to turn anastomosis on/off

	int izfibflat;

	int icirc;
};
