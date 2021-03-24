///////////////////////////////////////////////////////////////////////
// Filein.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The FILEIN class reads in the input file that contain parameters
// for the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include <fstream>

#include "FECore/FEModel.h"
#include "FEAngInput.h"
#include "mat3.h"

using namespace std;

class Filein
{
// Public functions:
public:
	// Constructor
	Filein();

	// Destructor
	virtual ~Filein();    

	// input
	bool Input(const char* szfile, FEAngInput &input, FEMesh &mesh);
    
	// read_param - Read a parameter from the current line within the buffer
    void read_param(FEAngInput &input, char (&buffer)[100]);

	// set_param - Set a certain parameter within the INPUT class based on what's read from the input file
    void set_param(FEAngInput &input, char (&buffer)[100], char (&pname)[20]);		
	
	// read_FEmesh - Read in the nodes and element connectivity from the FEBio mesh in order to construct the angio3d grid
	void read_FEmesh(FEAngInput &input, FEMesh &mesh);

	// enforce_fiber_BCS - If a node lies on a boundary face, adjust the collagen fiber orientation at that node accordingly
	void enforce_fiber_BCS(Node &node, FEAngInput &input, bool circ);
};
