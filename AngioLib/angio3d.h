///////////////////////////////////////////////////////////////////////
// angio3d.h
///////////////////////////////////////////////////////////////////////


#pragma once

#include <math.h>

//using namespace std;

const double pi = 3.14159265;                                   // Global constant pi
const double E = 2.718281828;                                   // Global constant e (Euler's number)


///////////////////////// LINEAR ALGEBRA FUNCTIONS /////////////////////////

inline double vec_dot(double a[3], double b[3])
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

inline double vec_norm(double a[3])
{
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

inline int sign(float n) { return (n>=0.f?1:-1); }



///////////////////////////// OTHER FUNCTIONS /////////////////////////////

inline double distance(double x0, double xf, double y0, double yf, double z0, double zf)
{
    return sqrt(pow(x0-xf,2.0)+pow(y0-yf,2.0)+pow(z0-zf,2.0));
}

// Insert function that constricts the angle into the domain +/- pi
inline double ang_dom(double angle)
{
    if (angle != 0.0)
    {
        while(1)
        {
            if (angle == pi)
                return angle;
            
            if (angle == -pi)
                return angle;
            
            if (angle > pi)
                angle = angle - 2*pi;
    
            if (angle < -pi)
                angle = angle + 2*pi;
        
            if ((angle <= pi) && (angle >= -pi))
                return angle;
        }
    }
    
    return angle;
}