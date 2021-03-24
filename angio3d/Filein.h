///////////////////////////////////////////////////////////////////////
// Filein.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// The FILEIN class reads in the input files that contain parameters
// for the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include <fstream>
#include "Input.h"

using namespace std;

//const char filename[13] = "in_input.ang";
const char filename[15] = "in_inputBC.ang";
//const char filename[16] = "in_inputLAC.ang";
//
const char gridname[14] = "in_nodes.ang";
const char econname[13] = "in_econn.ang";
const char eBCname[11] = "in_eBC.ang";
const char fibname[15] = "in_collfib.ang";
const char uname[12] = "in_disp.ang";
const char coll_den_name[16] = "in_coll_den.ang";

//const char gridname[15] = "sac_nodesi.ang";
//const char econname[14] = "sac_econn.ang";
//const char eBCname[12] = "sac_eBC.ang";
//const char fibname[16] = "sac_collfib.ang";
//const char uname[13] = "sac_disp.ang";


class Filein
{
public:
	Filein(Input &input);
	virtual ~Filein();    
    
    void read_param(Input &input, char (&buffer)[100]);
    void set_param(Input &input, char (&buffer)[100], char (&pname)[20]);		
	void read_nodes(Input &input);
	void read_coll_den(Input &input);
	void read_econn(Input &input);
	void read_eBC(Input &input);
	void read_collfib(Input &input);
	void read_disp(Input &input);

public:
    ifstream in_file;
    ifstream gr_file;
    ifstream ec_file;
    ifstream bc_file;
    ifstream fib_file;
    ifstream u_file;
};