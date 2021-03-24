///////////////////////////////////////////////////////////////////////
// Elem.cpp
///////////////////////////////////////////////////////////////////////



// Include:
#include "stdafx.h"
#include "Elem.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double Elem::bb_xmin()
{
    double exmin = bb_xmax();
    
    if ((*n1).x < exmin)
        exmin = (*n1).x;
        
    /*if ((*n2).x < exmin)
        exmin = (*n2).x;*/
        
    if ((*n3).x < exmin)
        exmin = (*n3).x;
        
    /*if ((*n4).x < exmin)
        exmin = (*n4).x;*/
        
    if ((*n5).x < exmin)
        exmin = (*n5).x;   
    
    /*if ((*n6).x < exmin)
        exmin = (*n6).x;*/
        
    if ((*n7).x < exmin)
        exmin = (*n7).x;
        
    /*if ((*n8).x < exmin)
        exmin = (*n8).x;*/

    return exmin;
}

double Elem::bb_xmax()
{
    double exmax = 0.0;
    
   /* if ((*n1).x > exmax)
        exmax = (*n1).x;*/
        
    if ((*n2).x > exmax)
        exmax = (*n2).x;
        
    /*if ((*n3).x > exmax)
        exmax = (*n3).x;*/
        
    if ((*n4).x > exmax)
        exmax = (*n4).x;
        
    /*if ((*n5).x > exmax)
        exmax = (*n5).x;  */ 
    
    if ((*n6).x > exmax)
        exmax = (*n6).x;
        
    /*if ((*n7).x > exmax)
        exmax = (*n7).x;*/
        
    if ((*n8).x > exmax)
        exmax = (*n8).x;

    return exmax;
}

double Elem::bb_ymin()
{
    double eymin = bb_ymax();
    
    if ((*n1).y < eymin)
        eymin = (*n1).y;
        
    if ((*n2).y < eymin)
        eymin = (*n2).y;
        
    /*if ((*n3).y < eymin)
        eymin = (*n3).y;
        
    if ((*n4).y < eymin)
        eymin = (*n4).y;*/
        
    if ((*n5).y < eymin)
        eymin = (*n5).y;   
    
    if ((*n6).y < eymin)
        eymin = (*n6).y;
        
    /*if ((*n7).y < eymin)
        eymin = (*n7).y;
        
    if ((*n8).y < eymin)
        eymin = (*n8).y;*/

    return eymin;
}

double Elem::bb_ymax()
{
    double eymax = 0.0;
    
    /*if ((*n1).y > eymax)
        eymax = (*n1).y;
        
    if ((*n2).y > eymax)
        eymax = (*n2).y;*/
        
    if ((*n3).y > eymax)
        eymax = (*n3).y;
        
    if ((*n4).y > eymax)
        eymax = (*n4).y;
        
    /*if ((*n5).y > eymax)
        eymax = (*n5).y;   
    
    if ((*n6).y > eymax)
        eymax = (*n6).y;*/
        
    if ((*n7).y > eymax)
        eymax = (*n7).y;
        
    if ((*n8).y > eymax)
        eymax = (*n8).y;

    return eymax;
}
 
double Elem::bb_zmin()
{
    double ezmin = bb_zmax();
    
    if ((*n1).z < ezmin)
        ezmin = (*n1).z;
        
    if ((*n2).z < ezmin)
        ezmin = (*n2).z;
        
    if ((*n3).z < ezmin)
        ezmin = (*n3).z;
        
    if ((*n4).z < ezmin)
        ezmin = (*n4).z;
        
    /*if ((*n5).z < ezmin)
        ezmin = (*n5).z;   
    
    if ((*n6).z < ezmin)
        ezmin = (*n6).z;
        
    if ((*n7).z < ezmin)
        ezmin = (*n7).z;
        
    if ((*n8).z < ezmin)
        ezmin = (*n8).z;*/

    return ezmin;
}

double Elem::bb_zmax()
{
    double ezmax = 0.0;
    
    /*if ((*n1).z > ezmax)
        ezmax = (*n1).z;
        
    if ((*n2).z > ezmax)
        ezmax = (*n2).z;
        
    if ((*n3).z > ezmax)
        ezmax = (*n3).z;
        
    if ((*n4).z > ezmax)
        ezmax = (*n4).z;*/
        
    if ((*n5).z > ezmax)
        ezmax = (*n5).z;   
    
    if ((*n6).z > ezmax)
        ezmax = (*n6).z;
        
    if ((*n7).z > ezmax)
        ezmax = (*n7).z;
        
    if ((*n8).z > ezmax)
        ezmax = (*n8).z;

    return ezmax;
}       

double Elem::ib_xmin()
{
    double exmax = (*n1).x;
    
    if ((*n3).x > exmax)
        exmax = (*n3).x;
        
    if ((*n5).x > exmax)
        exmax = (*n5).x;
        
    if ((*n7).x > exmax)
        exmax = (*n7).x;
    
    return exmax;
}


double Elem::ib_xmax()
{
    double exmin = (*n2).x;
    
    if ((*n4).x < exmin)
        exmin = (*n4).x;
        
    if ((*n6).x < exmin)
        exmin = (*n6).x;
        
    if ((*n8).x < exmin)
        exmin = (*n8).x;
    
    return exmin;
}

double Elem::ib_ymin()
{
    double eymax = (*n1).y;
    
    if ((*n2).y > eymax)
        eymax = (*n2).y;
        
    if ((*n5).y > eymax)
        eymax = (*n5).y;
        
    if ((*n6).y > eymax)
        eymax = (*n6).y;
    
    return eymax;
}

double Elem::ib_ymax()
{
    double eymin = (*n3).y;
    
    if ((*n4).y < eymin)
        eymin = (*n4).y;
        
    if ((*n7).y < eymin)
        eymin = (*n7).y;
        
    if ((*n8).y < eymin)
        eymin = (*n8).y;
    
    return eymin;
}

double Elem::ib_zmin()
{
    double ezmax = (*n1).z;
    
    if ((*n2).z > ezmax)
        ezmax = (*n2).z;
        
    if ((*n3).z > ezmax)
        ezmax = (*n3).z;
        
    if ((*n4).z > ezmax)
        ezmax = (*n4).z;
    
    return ezmax;
}

double Elem::ib_zmax()
{
    double ezmin = (*n5).z;
    
    if ((*n6).z < ezmin)
        ezmin = (*n6).z;
        
    if ((*n7).z < ezmin)
        ezmin = (*n7).z;
        
    if ((*n8).z < ezmin)
        ezmin = (*n8).z;
    
    return ezmin;
}