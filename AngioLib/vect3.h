///////////////////////////////////////////////////////////////////////
// vect3.h
///////////////////////////////////////////////////////////////////////

//// 3x1 vector class

#pragma once

#include <math.h>

class vect3
{
  public:
	//// vect3 - Contents of the vector
	double x, y, z;
        
  public:
	//// vect3 - Constructors
	vect3() : x(0.), y(0.), z(0.) {}
	vect3(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
	vect3(double (&iarray)[3]) : x(iarray[0]), y(iarray[1]), z(iarray[2]) {}
	
	
	//// vect3 - Operators
	vect3 operator + (const vect3& r) const { return vect3(x + r.x, y + r.y, z + r.z); }
	vect3 operator - (const vect3& r) const { return vect3(x - r.x, y - r.y, z - r.z); }

	vect3 operator * (double a) { return vect3(x*a, y*a, z*a); }
	vect3 operator / (double a) { return vect3(x/a, y/a, z/a); }

	vect3& operator += (const vect3& r) { x += r.x; y += r.y; z += r.z; return (*this); }
	vect3& operator -= (const vect3& r) { x -= r.x; y -= r.y; z -= r.z; return (*this); }

	vect3& operator *= (double a) { x*=a; y*=a; z*=a; return (*this); }
	vect3& operator /= (double a) { x/=a; y/=a; z/=a; return (*this); }

	vect3 operator - () { return vect3(-x, -y, -z); }

	bool operator == (const vect3& r)
	{
		if ((x == r.x) && (y == r.y) && (z == r.z))
			return true;
		else
			return false;
	}
		
	//// vect3 - Dot product
	double operator * (const vect3& r) { return (x*r.x + y*r.y + z*r.z); }

	
	//// vect3 - Cross product
	vect3 operator ^ (const vect3& r) { return vect3(y*r.z - z*r.y, z*r.x - x*r.z, x*r.y - y*r.x); }

	
	//// vect3 - Normalize the vector
	double unit()
	{
		double d = sqrt(x*x + y*y + z*z);
		if (d != 0) { x/=d; y/=d; z/=d; }
		return d;
	}
    
    
    //// vect3 - Norm
	double norm() { return sqrt(x*x + y*y + z*z); }
	
	
	//// vect3 - Reassign
	void reassign(double ix, double iy, double iz)
	{
	    x = ix;
	    y = iy;
	    z = iz;
	    
	    return;
	}
	
};
