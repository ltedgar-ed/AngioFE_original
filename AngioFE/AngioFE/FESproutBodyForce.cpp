///////////////////////////////////////////////////////////////////////
// FESproutBodyForce.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "FECore/FECoreKernel.h"
#include "FESproutBodyForce.h"
#include "FECore/FEModel.h"
#include "FEBioMech/FEElasticMaterial.h"
#include <cmath>

extern FECoreKernel* pFEBio;


//////////////////////////////////////////////////////////////////////
// Declare parameters
//////////////////////////////////////////////////////////////////////

BEGIN_PARAMETER_LIST(FESproutBodyForce, FEBodyForce);
	ADD_PARAMETER(m_a, FE_PARAM_DOUBLE, "a");
	ADD_PARAMETER(m_b, FE_PARAM_DOUBLE, "b");
//	ADD_PARAMETER(m_rc, FE_PARAM_VEC3D, "rc");
//	ADD_PARAMETER(m_sprout, FE_PARAM_VEC3D, "sprout");
	ADD_PARAMETER(m_inode, FE_PARAM_INT, "node");
	ADD_PARAMETER(m_brigid, FE_PARAM_BOOL, "rigid");
END_PARAMETER_LIST();


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FESproutBodyForce::FESproutBodyForce(FEModel* pfem) : FEBodyForce(pfem)
{
	m_factor = 1.0;
	m_pel = 0; 
	m_brigid = true; 
	m_inode = -1; 

	// Create symmetry vectors
	sym_planes[0] = 0; sym_planes[1] = 0; sym_planes[2] = 0; sym_planes[3] = 0; sym_planes[4] = 0; sym_planes[5] = 0; sym_planes[6] = 0;
	Sx = 0.; Sy = 0.; Sz = 0.;
	
	sym_vects[0][0] = 1.; sym_vects[0][1] = 0.; sym_vects[0][2] = 0.;
	sym_vects[1][0] = 0.; sym_vects[1][1] = 1.; sym_vects[1][2] = 0.;
	sym_vects[2][0] = 0.; sym_vects[2][1] = 0.; sym_vects[2][2] = 1.;
	sym_vects[3][0] = 1.; sym_vects[3][1] = 1.; sym_vects[3][2] = 0.;
	sym_vects[4][0] = 1.; sym_vects[4][1] = 0.; sym_vects[4][2] = 1.;
	sym_vects[5][0] = 0.; sym_vects[5][1] = 1.; sym_vects[5][2] = 1.;
	sym_vects[6][0] = 1.; sym_vects[6][1] = 1.; sym_vects[6][2] = 1.;

	sym_on = false;
}

///////////////////////////////////////////////////////////////////////
// Member Functions:
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
// FESproutBodyForce - force
//      Calculate the force at a partciular material point for each vessel sprout
///////////////////////////////////////////////////////////////////////

vec3d FESproutBodyForce::force(FEMaterialPoint& mp)
{
	//profiler.time_in(2);
	
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();		// Get the material point data
	vec3d x = pt.m_rt;															// Get the position of the material point 

	vec3d f(0,0,0);																// Initialize the force vector 
	
	int NS = (int) m_sp.size();													// Get the number of sprouts
	
#pragma omp parallel for shared(f)
	for (int i=0; i<NS; ++i)													// For each sprout...
	{
		SPROUT& sp = m_sp[i];														// Get the sprout
		
		if (sp.active)																// If the sprout is active...
		{
			vec3d r = sp.rc - x;														// Calculate the vector r between x and the sprout position				
			double l = r.unit();														// Get the length of r
			sp.sprout.unit();															// Normalize the sprout direction vector
			double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

			double g = m_a*(pow(cos(theta/2),m_factor))*exp(-m_b*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation
			//double g = m_a*exp(-m_b*l);

			vec3d fi = -r*g;
			if (sym_on == true)															// If symmetry is turned on, apply symmetry
				MirrorSym(x, fi, sp);

#pragma omp critical
 			f += fi;																	// Add results to the force vector
		}
	}
	
	//profiler.time_out();

	return f;
}


///////////////////////////////////////////////////////////////////////
// FESproutBodyForce - ApplySym
//      Determine if symmetry is turned on, if so create the symmetry vectors
///////////////////////////////////////////////////////////////////////

void FESproutBodyForce::ApplySym()
{
	if (Sx != 0)															// Turn on x symmetry
		sym_planes[0] = 1;

	if (Sy != 0)															// Turn on y symmetry
		sym_planes[1] = 1;

	if (Sz != 0)															// Turn on z symmetry
		sym_planes[2] = 1;

	if ((Sx != 0) && (Sy != 0))												// Turn on x and y symmetry
		sym_planes[3] = 1;

	if ((Sx != 0) && (Sz != 0))												// Turn on x and z symmetry
		sym_planes[4] = 1;

	if ((Sy != 0) && (Sz != 0))												// Turn on y and z symmetry
		sym_planes[5] = 1;

	if ((Sx != 0) && (Sy != 0) && (Sz != 0))								// Turn on x y and z symmetry
		sym_planes[6] = 1;
	
	if (sym_planes[0] + sym_planes[1] + sym_planes[2] + sym_planes[3] + sym_planes[4] + sym_planes[5] + sym_planes[6] != 0)
		sym_on = true;														

	return;
}


///////////////////////////////////////////////////////////////////////
// FESproutBodyForce - MirrorSym
//      Calculate force due to mirrored vessels at a particular material point at position x
///////////////////////////////////////////////////////////////////////

void FESproutBodyForce::MirrorSym(vec3d x, vec3d &f, SPROUT sp)
{
	sym.x = Sx; sym.y = Sy; sym.z = Sz;										// Set the position of the symmetry planes
	vec3d fsym;																// Mirrored sprout force vector
	vec3d sprout_vect;														// Sprout direction vector
	vec3d sym_v;															// Symmetry vector

	for (int i = 0; i < 7; i++){											// For each of the possible symmetry planes
		if (sym_planes[i] == 1){												// If that symmetry plane is turned on...
			sym_v.x = sym_vects[i][0]; sym_v.y = sym_vects[i][1]; sym_v.z = sym_vects[i][2];	// Obtain the symmetry vector

			vec3d r = sp.rc - x;													// Draw the vector r from the sprout location to the position of the material point
			
			r.x += 2*sym_v.x*(sym.x - sp.rc.x); r.y += 2*sym_v.y*(sym.y - sp.rc.y); r.z += 2*sym_v.z*(sym.z - sp.rc.z);		// Find r for the mirrored vessel sprout
			double l = r.unit();													// Find the length of r
			
			sprout_vect.x = sp.sprout.x; sprout_vect.y = sp.sprout.y; sprout_vect.z = sp.sprout.z;	// Set the sprout direction vector 
			sprout_vect.unit();														// Normalize the sprout direction vector
			//sp.sprout.unit();
			
			if (m_factor != 0){														// If a directional sprout force is being used...
				switch (i) {
				case 0:
					sprout_vect.x = -sprout_vect.x; break;								// Mirror across x
				case 1:
					sprout_vect.y = -sprout_vect.y; break;								// Mirror across y
				case 2:
					sprout_vect.z = -sprout_vect.z; break;								// Mirror across z
				case 3:
					sprout_vect.x = -sprout_vect.x; sprout_vect.y = -sprout_vect.y; break;	// Mirror across x and y
				case 4:
					sprout_vect.x = -sprout_vect.x; sprout_vect.z = -sprout_vect.z; break;	// Mirror across x and z
				case 5:
					sprout_vect.y = -sprout_vect.y; sprout_vect.z = -sprout_vect.z; break;	// Mirror across y and z
				case 6:
					sprout_vect.x = -sprout_vect.x; sprout_vect.y = -sprout_vect.y; sprout_vect.z = -sprout_vect.z; break;}	// Mirror across x y and z
			}
	
			double theta = acos(sprout_vect*r);										// Calculate theta, the angle between r and the sprout direction vector

			double g = m_a*(pow(cos(theta/2),m_factor))*exp(-m_b*l);				// Calculate the magnitude of the sprout force using the localized directional sprout force equation			
			//double g = m_a*exp(-m_b*l);
			
			if ((g != g) || (r.x != r.x) || (r.y != r.y) || (r.z != r.z)){			// If the mirrored force vector isn't realy...
				g = 0.; r.x = 0.; r.y = 0.; r.z = 0.;}									// Set it to zero

			fsym += -r*g;															// Add the results to the force symmetry vector
		}
	}

	f += fsym;															// Add the symmetry results to the force vector

	return;
}


///////////////////////////////////////////////////////////////////////
// FESproutBodyForce - stiffness
//      
///////////////////////////////////////////////////////////////////////

mat3ds FESproutBodyForce::stiffness(FEMaterialPoint &mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	vec3d x = pt.m_rt;

	mat3ds k; k.zero();
	
	//int NS = (int) m_sp.size();
	
	//for (int i=0; i<NS; ++i)
	//{
	//	SPROUT& sp = m_sp[i];
	//	if (sp.active)
	//	{
	//		vec3d r = x - sp.rc;
	//		double l = r.unit();
	//		sp.sprout.unit();
	//		double theta = acos(sp.sprout*r);	
	//
	//		if (l != 0.0)
	//		{
	//			double g = m_a*(pow(cos(theta/2),m_factor))*exp(-m_b*l);
	//			//double g = m_a*exp(-m_b*l);

	//			mat3ds rxr = dyad(r);
	//			mat3ds I = mat3dd(1.0);

	//			//k += (rxr*m_b - (I - rxr)/l)*g;
	//		}
	//	}
	//}

	return k;
}


///////////////////////////////////////////////////////////////////////
// FESproutBodyForce - Serialize
//      
///////////////////////////////////////////////////////////////////////

void FESproutBodyForce::Serialize(DumpFile &ar)
{
	if (ar.IsSaving())
	{
		ar << m_a << m_b;
		ar << m_inode << m_brigid;
	}
	else
	{
		ar >> m_a >> m_b;
		ar >> m_inode >> m_brigid;
	}
}


///////////////////////////////////////////////////////////////////////
// FESproutBodyForce - Init
//      
///////////////////////////////////////////////////////////////////////

bool FESproutBodyForce::Init()
{
/*
	if (m_inode == -1)
	{
		if (!m_brigid)
		{
			// find the element in which point r0 lies
			FEMesh& m = m_pfem->m_mesh;
			m_pel = m.FindSolidElement(m_rc, m_rs);
		}
		else m_pel = 0;
	}
	else 
	{
		FEMesh& m = m_pfem->m_mesh;
		m_rc = m.Node(m_inode).m_r0;
	}
	*/
	return true;
}



///////////////////////////////////////////////////////////////////////
// FESproutBodyForce - Update
//      
///////////////////////////////////////////////////////////////////////

void FESproutBodyForce::Update()
{
/*	if (m_inode == -1)
	{
		if (m_pel)
		{
			FEMesh& m = m_pfem->m_mesh;
			vec3d x[8];
			for (int i=0; i<8; ++i) x[i] = m.Node(m_pel->m_node[i]).m_rt;

			double* r = m_rs;
			double H[8];
			H[0] = 0.125*(1 - r[0])*(1 - r[1])*(1 - r[2]);
			H[1] = 0.125*(1 + r[0])*(1 - r[1])*(1 - r[2]);
			H[2] = 0.125*(1 + r[0])*(1 + r[1])*(1 - r[2]);
			H[3] = 0.125*(1 - r[0])*(1 + r[1])*(1 - r[2]);
			H[4] = 0.125*(1 - r[0])*(1 - r[1])*(1 + r[2]);
			H[5] = 0.125*(1 + r[0])*(1 - r[1])*(1 + r[2]);
			H[6] = 0.125*(1 + r[0])*(1 + r[1])*(1 + r[2]);
			H[7] = 0.125*(1 - r[0])*(1 + r[1])*(1 + r[2]);

			m_rc = vec3d(0,0,0);
			for (int i=0; i<8; ++i) m_rc += x[i]*H[i];
		}
	}
	else
	{
		FEMesh& m = m_pfem->m_mesh;
		m_rc = m.Node(m_inode).m_rt;
	}
*/
}


