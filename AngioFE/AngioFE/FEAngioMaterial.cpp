#include "StdAfx.h"
#include "FEAngioMaterial.h"
#include <FECore/FEModel.h>
#include <FECore/FESolidDomain.h>
#include <cmath>
#include <Elem.h>

//-----------------------------------------------------------------------------
BEGIN_PARAMETER_LIST(FEAngioMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_a, FE_PARAM_DOUBLE, "a");
	ADD_PARAMETER(m_b, FE_PARAM_DOUBLE, "b");
	ADD_PARAMETER(m_N, FE_PARAM_DOUBLE, "N");
	ADD_PARAMETER(m_s, FE_PARAM_VEC3D , "sprout");
END_PARAMETER_LIST();

//-----------------------------------------------------------------------------
FEAngioMaterial::FEAngioMaterial(FEModel* pfem) : FEElasticMaterial(pfem)
{
	scale = 1.0;

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

//-----------------------------------------------------------------------------
void FEAngioMaterial::Init()
{
	// add the user sprouts
	for (int i=0; i<m_suser.size(); ++i)
	{
		AddSprout(m_suser[i], vec3d(0,0,0));
	}
	m_suser.clear();
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::SetParameter(FEParam& p)
{
	if (strcmp(p.m_szname, "sprout") == 0)
	{
		m_suser.push_back(m_s);
	}
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::AddSprout(const vec3d& r, const vec3d& t)
{
	SPROUT s;
	s.bactive = true;
	s.sprout = t;
	UpdateSprout(s, r);
	m_spr.push_back(s);
}

//-----------------------------------------------------------------------------
void FEAngioMaterial::UpdateSprout(SPROUT& s, const vec3d& r)
{
	FEMesh& mesh = GetFEModel()->GetMesh();
	
	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));
	s.pel = dom.FindElement(r, s.r);
}

//-----------------------------------------------------------------------------
vec3d FEAngioMaterial::CurrentSproutPosition(FEAngioMaterial::SPROUT& p)
{
	double r = p.r[0], s = p.r[1], t = p.r[2];
	double H[8];
	H[0] = 0.125*(1.0 - r)*(1.0 - s)*(1.0 - t);
	H[1] = 0.125*(1.0 + r)*(1.0 - s)*(1.0 - t);
	H[2] = 0.125*(1.0 + r)*(1.0 + s)*(1.0 - t);
	H[3] = 0.125*(1.0 - r)*(1.0 + s)*(1.0 - t);
	H[4] = 0.125*(1.0 - r)*(1.0 - s)*(1.0 + t);
	H[5] = 0.125*(1.0 + r)*(1.0 - s)*(1.0 + t);
	H[6] = 0.125*(1.0 + r)*(1.0 + s)*(1.0 + t);
	H[7] = 0.125*(1.0 - r)*(1.0 + s)*(1.0 + t);

	FEMesh& mesh = GetFEModel()->GetMesh();
	FEElement& el = *p.pel;
	vec3d rc(0,0,0);
	for (int i=0; i<8; ++i) rc += mesh.Node(el.m_node[i]).m_rt*H[i];

	return rc;
}

//-----------------------------------------------------------------------------
mat3ds FEAngioMaterial::SproutStress(const vec3d& y)
{
	double den_scale = 1.0;
	den_scale = findDenScale(y.x, y.y, y.z);

	// loop over all sprout tips
	mat3ds s(0.0);
	int NS = Sprouts();

//#pragma omp parallel for shared(s)
	for (int i=0; i<NS; ++i)
	{
		SPROUT& sp = m_spr[i];
		
		if (sp.pel && sp.bactive)
		{
			vec3d x = CurrentSproutPosition(sp);

			vec3d r = y - x;
			double l = r.unit();

			sp.sprout.unit();															// Normalize the sprout direction vector
			
			double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

			double p = den_scale*scale*m_a*(pow(cos(theta/2),m_N))*exp(-m_b*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation

			//double p = m_a*exp(-m_b*l);

			mat3ds si = dyad(r)*p;

			if (sym_on == true)															// If symmetry is turned on, apply symmetry
				MirrorSym(y, si, sp, den_scale);

//#pragma omp critical
			s += si;
		}
	}
	return s;
}

//-----------------------------------------------------------------------------
mat3ds FEAngioMaterial::Stress(FEMaterialPoint& mp)
{
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
	return SproutStress(pt.m_rt);
}

//-----------------------------------------------------------------------------
tens4ds FEAngioMaterial::Tangent(FEMaterialPoint& mp)
{
	tens4ds C(0.0);
	return C;
}

//=============================================================================
BEGIN_PARAMETER_LIST(FEPressureMaterial, FEElasticMaterial)
	ADD_PARAMETER(m_p, FE_PARAM_DOUBLE, "p");
END_PARAMETER_LIST();

mat3ds FEPressureMaterial::Stress(FEMaterialPoint& pt)
{
	mat3dd I(1.0);
	return I*m_p;
}

tens4ds FEPressureMaterial::Tangent(FEMaterialPoint& pt)
{
	return tens4ds(0.0);
}

//============================================================================

///////////////////////////////////////////////////////////////////////
// FEAngioMaterial - ApplySym
//      Determine if symmetry is turned on, if so create the symmetry vectors
///////////////////////////////////////////////////////////////////////

void FEAngioMaterial::ApplySym()
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
// FEAngioMaterial - MirrorSym
//      Calculate force due to mirrored vessels at a particular material point at position x
///////////////////////////////////////////////////////////////////////

void FEAngioMaterial::MirrorSym(vec3d y, mat3ds &s, SPROUT sp, double den_scale)
{
	sym.x = Sx; sym.y = Sy; sym.z = Sz;										// Set the position of the symmetry planes
	mat3ds ssym; ssym.zero();
	vec3d sprout_vect;														// Sprout direction vector
	vec3d sym_v;															// Symmetry vector

	for (int i = 0; i < 7; i++){											// For each of the possible symmetry planes
		if (sym_planes[i] == 1){												// If that symmetry plane is turned on...
			sym_v.x = sym_vects[i][0]; sym_v.y = sym_vects[i][1]; sym_v.z = sym_vects[i][2];	// Obtain the symmetry vector

			vec3d x = CurrentSproutPosition(sp);
			vec3d r = y - x;																	// Draw the vector r from the sprout location to the position of the material point
			
			//r.x += 2*sym_v.x*(sym.x - x.x); r.y += 2*sym_v.y*(sym.y - x.y); r.z += 2*sym_v.z*(sym.z - x.z);		// Find r for the mirrored vessel sprout
			r.x = r.x + sym_v.x*sym.x; r.y = r.y + sym_v.y*sym.y; r.z = r.z + sym_v.z*sym.z;
			double l = r.unit();													// Find the length of r
			
			sprout_vect.x = sp.sprout.x; sprout_vect.y = sp.sprout.y; sprout_vect.z = sp.sprout.z;	// Set the sprout direction vector 
			sprout_vect.unit();														// Normalize the sprout direction vector
			
			if (m_N != 0){														// If a directional sprout force is being used...
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
	
			double theta = acos(sp.sprout*r);											// Calculate theta, the angle between r and the sprout vector

			double p = den_scale*scale*m_a*(pow(cos(theta/2),m_N))*exp(-m_b*l);					// Calculate the magnitude of the sprout force using the localized directional sprout force equation
			
			if ((p != p) || (r.x != r.x) || (r.y != r.y) || (r.z != r.z)){			// If the mirrored force vector isn't real...
				p = 0.; r.x = 0.; r.y = 0.; r.z = 0.;}									// Set it to zero

			ssym += dyad(r)*p;
		}
	}


	s += ssym;															// Add the symmetry results to the force vector

	return;
}

double FEAngioMaterial::findDenScale(double xpt, double ypt, double zpt)
{
	double coll_den = 0.0;
    double den_scale = 1.0;

    double xix, xiy, xiz;
    double shapeF[8];

	int elem_num = pgrid->findelem(xpt, ypt, zpt);
    
    if (elem_num < 0)
		return den_scale;
		
	Elem elem;
    elem = pgrid->ebin[elem_num];
        
    // Convert to natural coordinates -1 <= Xi <= +1
    pgrid->natcoord(xix, xiy, xiz, xpt, ypt, zpt, elem_num);        
    
    // Obtain shape function weights
    pgrid->shapefunctions(shapeF, xix, xiy, xiz);
    
	coll_den = shapeF[0]*(*elem.n1).ecm_den + shapeF[1]*(*elem.n2).ecm_den + shapeF[2]*(*elem.n3).ecm_den + shapeF[3]*(*elem.n4).ecm_den + shapeF[4]*(*elem.n5).ecm_den + shapeF[5]*(*elem.n6).ecm_den + shapeF[6]*(*elem.n7).ecm_den + shapeF[7]*(*elem.n8).ecm_den;
    
	den_scale = pgrid->find_density_scale(coll_den);

    if (den_scale < 0.)
		den_scale = 0.;
	
	return den_scale;
}