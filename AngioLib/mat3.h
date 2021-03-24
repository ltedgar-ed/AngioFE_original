///////////////////////////////////////////////////////////////////////
// mat3.h
///////////////////////////////////////////////////////////////////////

//// 3x3 matrix class

#pragma once

#include "vect3.h"

class mat3
{
  public:
    //// mat3 - Contents of the matrix
    double M11, M12, M13, M21, M22, M23, M31, M32, M33;
    
  public:
    //// mat3 - Constructors
    mat3() : M11(0.), M12(0.), M13(0.), M21(0.), M22(0.), M23(0.), M31(0.), M32(0.), M33(0.)  {}
    //mat3(vector< vector<double> >& idata) : M11(idata[0][0]), M12(idata[0][1]), M13(idata[0][2]), M21(idata[1][0]), M22(idata[1][1]), M23(idata[1][2]), M31(idata[2][0]), M32(idata[2][1]), M33(idata[2][2]) {}
    mat3(double (&iarray)[9]) : M11(iarray[0]), M12(iarray[1]), M13(iarray[2]), M21(iarray[3]), M22(iarray[4]), M23(iarray[5]), M31(iarray[6]), M32(iarray[7]), M33(iarray[8]) {}
    mat3(double iM11, double iM12, double iM13, double iM21, double iM22, double iM23, double iM31, double iM32, double iM33) : M11(iM11), M12(iM12), M13(iM13), M21(iM21), M22(iM22), M23(iM23), M31(iM31), M32(iM32), M33(iM33)  {} 
    
    
    //// mat3 - Operators
    mat3 operator + (const mat3& R) const { return mat3(M11 + R.M11, M12 + R.M12, M13 + R.M13, M21 + R.M21, M22 + R.M22, M23 + R.M23, M31 + R.M31, M32 + R.M32, M33 + R.M33); }
    mat3 operator - (const mat3& R) const { return mat3(M11 - R.M11, M12 - R.M12, M13 - R.M13, M21 - R.M21, M22 - R.M22, M23 - R.M23, M31 - R.M31, M32 - R.M32, M33 + R.M33); }
    
    mat3 operator * (double a) { return mat3(M11*a, M12*a, M13*a, M21*a, M22*a, M23*a, M31*a, M32*a, M33*a); }
    mat3 operator / (double a) { return mat3(M11/a, M12/a, M13/a, M21/a, M22/a, M23/a, M31/a, M32/a, M33/a); }  
        
    vect3 operator * (const vect3& v) const { return vect3(M11*v.x + M12*v.y + M13*v.z, M21*v.x + M22*v.y + M23*v.z, M31*v.x + M32*v.y + M33*v.z); }  
    mat3 operator * (const mat3& R) const { return mat3(M11*R.M11 + M12*R.M21 + M13*R.M31, M11*R.M12 + M12*R.M22 + M13*R.M32, M11*R.M13 + M12*R.M23 + M13*R.M33, M21*R.M11 + M22*R.M21 + M23*R.M31, M21*R.M12 + M22*R.M22 + M23*R.M32, M21*R.M13 + M22*R.M23 + M23*R.M33, M31*R.M11 + M32*R.M21 + M33*R.M31, M31*R.M12 + M32*R.M22 + M33*R.M32, M31*R.M13 + M32*R.M23 + M33*R.M33);}
    
    mat3 operator - () { return mat3(-M11, -M12, -M13, -M21, -M22, -M23, -M31, -M32, -M33); }
    	
    
    //// mat3 - Determinant
    double det()
    {
        double Mdet = M11*(M22*M33 - M23*M32) - M12*(M21*M33 - M23*M31) + M13*(M21*M32 - M22*M31);
        return Mdet;
    }
    
    
    //// mat3 - Invert
    mat3 invert()
    {
        mat3 inverse;      
        double Mdet = det();
        
        if (Mdet != 0)
            inverse.reassign((1/Mdet)*(M22*M33 - M23*M32), (1/Mdet)*(M13*M32 - M12*M33), (1/Mdet)*(M12*M23 - M13*M22), (1/Mdet)*(M23*M31 - M21*M33), (1/Mdet)*(M11*M33 - M13*M31),  (1/Mdet)*(M13*M21 - M11*M23),  (1/Mdet)*(M21*M32 - M22*M31),  (1/Mdet)*(M12*M31 - M11*M32),  (1/Mdet)*(M11*M22 - M12*M21));
            
        return inverse;
    }  
    
    //// mat3 - Transpose
	mat3 transpose()
	{
		mat3 transpose;

		transpose.M11 = M11;
		transpose.M22 = M22;
		transpose.M33 = M33;
		transpose.M12 = M21;
		transpose.M13 = M31;
		transpose.M21 = M12;
		transpose.M23 = M32;
		transpose.M31 = M13;
		transpose.M32 = M23;

		return transpose;
	}




    //// mat3 - Reassign
    void reassign(double iM11, double iM12, double iM13, double iM21, double iM22, double iM23, double iM31, double iM32, double iM33)
    {
        M11 = iM11;
        M12 = iM12;
        M13 = iM13;
        M21 = iM21;
        M22 = iM22;
        M23 = iM23;
        M31 = iM31;
        M32 = iM32;
        M33 = iM33;
    
        return;
    }
        
        
    //// mat3 - solve
    vect3 solve(vect3 b)
    {
        vect3 x;
        double Mdet = det();
        
        if (Mdet != 0)
        {
            mat3 Ainv = invert();      
            x = Ainv*b;
        }
        
        return x;
    }       
        
};      