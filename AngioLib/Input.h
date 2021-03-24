///////////////////////////////////////////////////////////////////////
// Input.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// The INPUT class reads in the input file and stores the parameter
// values used in the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include <iostream>
#include <vector>

class Node;

using namespace std;

class Input
{
public:
	Input();
	virtual ~Input();
			  
    double ibrnch_ch;
    int imatx_cnd;
    double imatx_den;
    int icomp_mat;

    int infrag;
    double imax_time;
    double idt;
    double ianst_dst;  
    double ilngth_adj;
    
    int ixnodes;
    int iynodes;
    int iznodes;  
    
    double ixmin;
    double ixmax;
    double iymin;
    double iymax;
    double izmin;
    double izmax;
    
    char ixbctype;                              // Boundary condition types:
    char iybctype;                              //  'w' - flat wall
    char izbctype;                              //  'p' - periodic 
                                                //  'b' - bouncy wall
    char ifrontbc;
	char irightbc;
	char ibackbc;
	char ileftbc;
	char ibottombc;
	char itopbc;
	
	char igrid_in;
    
    int iNn;                                            
    vector<Node> inodes;
    
    int iNe;
    vector<vector<int> > ieconn;
    
    int iNBC;
    vector<vector<int> > ieBC;
    
	vector<double> icoll_den;

	double iweight1;
	double iweight4;

};
