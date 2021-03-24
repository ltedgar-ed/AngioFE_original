// AngioFE.cpp : Defines the exported functions for the DLL application.
//

#include "StdAfx.h"

#include "FECore/FECoreKernel.h"
#include "AngioFETask.h"
#include "FECore/FECoreFactory.h"
#include "FESproutBodyForce.h"
#include "FEAngioMaterial.h"
#include "AngioPlot.h"

FECoreKernel* pFEBio;

//-----------------------------------------------------------------------------
FEPluginFactory_T<AngioFETask               , FETASK_ID    > angiofe_task_factory ("angio" );
FEPluginFactory_T<FESproutBodyForce         , FEBODYLOAD_ID> angio_sprout_factory ("sprout");
FEPluginFactory_T<FEAngioMaterial           , FEMATERIAL_ID> angio_mat_factory    ("angio" );
FEPluginFactory_T<FEPressureMaterial        , FEMATERIAL_ID> pressure_mat_factory ("pressure" );
FEPluginFactory_T<FEPlotAngioStress         , FEPLOTDATA_ID> plot_angio_stress    ("angio stress");
FEPluginFactory_T<FEPlotAngioEffectiveStress, FEPLOTDATA_ID> plot_angio_eff_stress("angio effective stress" );

//-----------------------------------------------------------------------------
FECORE_EXPORT unsigned int GetSDKVersion()
{
	return FE_SDK_VERSION;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT  void PluginInitialize(FECoreKernel& febio)
{
	pFEBio = &febio;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT int PluginNumClasses()
{
	return 6;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT FECoreFactory * PluginGetFactory(int i)
{
	if      (i==0) return &angiofe_task_factory;
	else if (i==1) return &angio_sprout_factory;
	else if (i==2) return &angio_mat_factory;
	else if (i==3) return &pressure_mat_factory;
	else if (i==4) return &plot_angio_stress;
	else if (i==5) return &plot_angio_eff_stress;
	else return 0;
}

//-----------------------------------------------------------------------------
FECORE_EXPORT void PluginCleanup()
{

}
