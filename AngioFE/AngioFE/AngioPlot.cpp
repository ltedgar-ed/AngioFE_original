#include "StdAfx.h"
#include "AngioPlot.h"
#include "FEAngio.h"
#include "FECore/FESolidDomain.h"

extern FEAngio* pfeangio;

//-----------------------------------------------------------------------------
bool FEPlotAngioStress::Save(FEDomain& d, vector<float>& a)
{
	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s(0.0);
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint& pt = *(el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>());
			mat3ds sj = (pfeangio ? pfeangio->m_pmat->SproutStress(pt.m_rt) : mat3ds(0.));
			s += sj;
		}
		s /= (double) nint;

		a.push_back((float) s.xx());
		a.push_back((float) s.yy());
		a.push_back((float) s.zz());
		a.push_back((float) s.xy());
		a.push_back((float) s.yz());
		a.push_back((float) s.xz());
	}
	return true;
};

//-----------------------------------------------------------------------------
bool FEPlotAngioEffectiveStress::Save(FEDomain& d, vector<float>& a)
{
	FESolidDomain& dom = dynamic_cast<FESolidDomain&>(d);
	int NE = dom.Elements();
	for (int i=0; i<NE; ++i)
	{
		FESolidElement& el = dom.Element(i);
		int nint = el.GaussPoints();
		mat3ds s(0.0);
		for (int j=0; j<nint; ++j)
		{
			FEElasticMaterialPoint& pt = *(el.GetMaterialPoint(j)->ExtractData<FEElasticMaterialPoint>());
			mat3ds sj = (pfeangio ? pfeangio->m_pmat->SproutStress(pt.m_rt) : mat3ds(0.));
			s += pt.m_s - sj;
		}
		s /= (double) nint;

		a.push_back((float) s.xx());
		a.push_back((float) s.yy());
		a.push_back((float) s.zz());
		a.push_back((float) s.xy());
		a.push_back((float) s.yz());
		a.push_back((float) s.xz());
	}
	return true;
};
