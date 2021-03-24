///////////////////////////////////////////////////////////////////////
// Input.cpp
///////////////////////////////////////////////////////////////////////



#include "stdafx.h"
#include "Input.h"
#include "Angio.h"
#include "angio3d.h"
#include <iostream>
#include <fstream>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Input::Input()
{   
    ibrnch_ch = 0.1;
    imatx_cnd = 0;
    imatx_den = 3.0;
    icomp_mat = 0;

    infrag = 3;
    imax_time = 0.0;
    idt = 0.25;
    ianst_dst = 75.0;  
    ilngth_adj = 1.0;
    
	iNn = 1;
	iNe = 1;
	
    ixnodes = 0;
    iynodes = 0;
    iznodes = 0;
    
    ixbctype = 'w';
    iybctype = 'w';
    izbctype = 'w'; 
    
	ifrontbc = 'w';
	irightbc = 'w';
	ibackbc = 'w';
	ileftbc = 'w';
	ibottombc = 'w';
	itopbc = 'w';
	
    igrid_in = 'n';

    iweight1 = 0.;
    iweight4 = 0.;	
}

Input::~Input()
{

}


