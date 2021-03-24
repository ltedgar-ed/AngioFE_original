// AngioFERun.cpp

////////////////////////////////////////////////////////////////////////////////////////////////////
////																							////            
////			 ANGIOFE - A plugin for FEBio that simulates the deformation of the             ////		 
////                    matrix produced by microvessels during angiogenesis.                    ////
////																							////
////							      By Lowell Taylor Edgar                                    ////  
////    						     University of Utah, 2011        							////	
////																							////			
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "time.h"
#include <list>

#include "AngioFETask.h"
#include "Filein.h"
#include "Fileout.h"

#include "Grid.h"
#include "FEAngInput.h"
#include "Data.h"
#include "Segment.h"
#include "Angio.h"
#include "FECore/FEModel.h"
#include "FEAngio.h"
#include "FESproutBodyForce.h"
#include "FECore/FECoreKernel.h"
#include "FEBioMech/FEElasticMixture.h"
#include "FEBioMech/FEViscoElasticMaterial.h"
#include "FEAngioMaterial.h"
#include <omp.h>

extern FECoreKernel* pFEBio;

FEAngio* pfeangio = 0;

using namespace std;

//-----------------------------------------------------------------------------
// find the angio material component
FEAngioMaterial* FindAngioMaterial(FEMaterial* pm)
{
	// first, extract the elastic component
	FEElasticMaterial* pme = dynamic_cast<FEElasticMaterial*>(pm->GetElasticMaterial());
	if (pme == 0) return 0;

	//// for a visco-material we have to go one level deeper
	//FEViscoElasticMaterial* pmv = dynamic_cast<FEViscoElasticMaterial*>(pme);
	//if (pmv) pme = pmv->GetBaseMaterial();

	// an angio material must be used inside a mixture
	FEElasticMixture* pmx = dynamic_cast<FEElasticMixture*>(pme);
	if (pmx == 0) return 0;

	// loop over all the mixture materials until we find the angio material
	for (int i=0; i<pmx->Materials(); ++i)
	{
		FEAngioMaterial* pma = dynamic_cast<FEAngioMaterial*>(pmx->GetMaterial(i));
		if (pma) return pma;
	}

	// no luck, so report the sad news
	return 0;
}

//-----------------------------------------------------------------------------
bool AngioFETask::Init(const char* inpfile)
{
	//// FEBIO - Retrieve the FE model
	FEModel& fem = *GetFEModel();	// Get the femodel
	FEMesh& mesh = fem.GetMesh();	// Get the mesh						

	//// ANGIO3D - Initialize growth model	
	FEAngInput input;													
	Filein filein;
	if (filein.Input(inpfile, input, mesh) == false) return false;	// Read in the angio3d input file
	
	time_t ranseed; ranseed = time(NULL);							// Create the seed for the random number generator
	if (input.irseed != 0){ ranseed = input.irseed;}				// If specified in the input file, read in the seed number
	
	//ranseed = 1234567890;											// Verification - Use predetermined seed number
	srand((unsigned int)ranseed);									// Seed the random number generator											
								
	pfeangio = new FEAngio(input);									// Create the FEAngio class
	FEAngio& feangio = *pfeangio;											
	feangio.fileout.timestart();
	
	//// FEBIO - Set parameters for the finite element model
	FEAnalysis* pstep = fem.GetCurrentStep();
	pstep->m_dt = pstep->m_dt0;							// Set the current time step to the initial time step
	pstep->m_ntimesteps = 0;									// Set the number of time steps to 0
	pstep->m_iteopt = 100;									// Set the optimal iterations to 100
	pstep->m_maxretries = 10;									// Set the maximum retries for the autotimestepper to 10
	pstep->m_dt = feangio.FE_time_step;

//--> SAM
	bool bmat = true; // use material approach or body-force approach
	FEElasticMixture* pmix = dynamic_cast<FEElasticMixture*>(m_pfem->GetMaterial(0)->GetElasticMaterial());
	if (pmix == 0) bmat = false;
	else
	{
		FEAngioMaterial* pm = dynamic_cast<FEAngioMaterial*>(pmix->GetMaterial(1));
		if (pm == 0) bmat = false;
		else feangio.m_pmat = pm;
	}

	FEAngioMaterial* pma = FindAngioMaterial(m_pfem->GetMaterial(0));
	if (pma == 0) bmat = false;
	else 
	{
		//pma->scale = (1.0/4.0)*0.001;

		if (input.iSx != 0.){ pma->Sx = input.iSx;}						// If there's an x-symmetry plane, set the location of the x-symmetry plane 
		if (input.iSy != 0.){ pma->Sy = input.iSy;}						// If there's an y-symmetry plane, set the location of the y-symmetry plane
		if (input.iSz != 0.){ pma->Sz = input.iSz;}						// If there's an z-symmetry plane, set the location of the z-symmetry plane

		pma->ApplySym();

		feangio.m_pmat = pma;
		feangio.update_sprout_stress_scaling();
	}
	
	if (bmat == false)
	{
//<-- SAM
		FESproutBodyForce* pbf = new FESproutBodyForce(m_pfem);			// Define the one-and-only bodyforce	
		feangio.m_pbf = pbf; // --- SAM ---
	
		m_pfem->AddBodyLoad(pbf);										// Add the body force to the FEmodel
		FEParameterList& pl = pbf->GetParameterList();					// Get the body force's parameter list
		FEParam* pa = pl.Find("a"); assert(pa);							// Get the magnitude parameter
		FEParam* pb = pl.Find("b"); assert(pb);							// Get the range parameter
	
		pa->value<double>() = (1.0/4.0)*0.001*feangio.sproutf;			// Set the mangnitude parameter using the input file
		pb->value<double>() = feangio.b;								// Set the range parameter using the input file
		pbf->m_factor = input.ispfactor;								// Set the directional factor parameter using the input file
	
		if (input.isp_sphere == 1) pbf->m_factor = 0.;					// If using local isotropic sprout force, set directional factor to 0

		if (input.isp_sphere == 2){										// If using global isotropic sprout froce, set directional factor and range to 0 
			pbf->m_factor = 0.; 
			pb->value<double>() = 0;}					

		//cout << pbf->m_factor << endl;									// Print out the directional factor parameter

		if (input.iSx != 0.){ pbf->Sx = input.iSx;}						// If there's an x-symmetry plane, set the location of the x-symmetry plane 
		if (input.iSy != 0.){ pbf->Sy = input.iSy;}						// If there's an y-symmetry plane, set the location of the y-symmetry plane
		if (input.iSz != 0.){ pbf->Sz = input.iSz;}						// If there's an z-symmetry plane, set the location of the z-symmetry plane

		pbf->ApplySym();												// Apply any symmetry to the bodyforce class
	}

	if (input.icirc)
		feangio.circ_gel();

	//// ANGIO3D - Seed initial vessel frags
	feangio.fileout.printrandseed((int)ranseed);					// Print out the seed number for the random generator
	feangio.initialize_grid_volume();								// Initialize the volume of each element in the grid
	feangio.seedFrag();												// Create initial microvessel fragments                               
	feangio.kill_dead_segs();										// Remove buggy segments
	feangio.find_active_tips();										// Update the active tip container

	if (input.isprout_verify == 1)									// If specified, run the sprout verification problem
		feangio.setup_sprout_verification();
	
	
	//// FEBIO - Initialize the FE model
	feangio.m_pmat->pgrid = &feangio.grid;
	feangio.FE_state++;												// Update the FE state

	// Do the model initialization
	if (fem.Init() == false) return false;							// Initialize the FE model

	return true;
}

//-----------------------------------------------------------------------------
bool AngioFETask::Run()
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();
	FEAngio& feangio = *pfeangio;
	FEAnalysis* pstep = fem.GetCurrentStep();

	feangio.save_vessel_state();									// Output microvessel state
	feangio.save_time();											// Output time information
	feangio.fileout.writeCollFib(feangio.grid, true);				// Output initial collagen fiber orientation

	//// FEBIO - Solve for initial step
	feangio.FE_state++;												// Update the FE state
	fem.Solve();													// Solve the FE model
	feangio.update_grid(mesh);										// Update the position of grid nodes within angio3d using the solution from FEBio
	feangio.update_ECM();											// Update the ECM based on the solution from FEBio (collagen fiber orientation and matrix density)
	feangio.update_grid_volume();									// Update the volume for each element in the grid
	feangio.displace_vessels();										// Displace the microvessels using the solution from FEBio
	pstep->m_dt = pstep->m_dt0;
	feangio.FE_time_step = pstep->m_dt;

	//// ANGIO3D - Apply sprout forces at active growth tips
	feangio.apply_sprout_forces(fem, 1, 0.5);						// Apply sprout forces to the mesh
	feangio.adjust_mesh_stiffness(fem);								// Adjust the stiffness of the mesh based on microvessel volume
	feangio.save_vessel_state();									// Output microvessel state
	//feangio.save_bdy_forces(fem);									// Output body force state
	feangio.save_time();											// Output time

	//// ANGIO3D - Simulate angiogenesis
	feangio.initBranch();											// Handle branching within inital fragments
    feangio.updateTotalLength();									// Update the total vascular length within the simulation   
	
	while (feangio.data.t < feangio.data.maxt)                      // While culture time is less than the maximum time
	{	
		feangio.updateTime();											// Determine the current time step and update time

		feangio.updateLength();											// Determine length of new segments for this growth step        		

		feangio.m_pmat->pgrid = &feangio.grid;
		
		feangio.Growth(fem);											// Growth (Elongation, Branching, and Anastomosis)
    	
		feangio.kill_dead_segs();										// Remove buggy segments
		
		feangio.updateTotalLength();									// Update total vascular length
	    
		feangio.find_active_tips();										// Locate all active growth tips
		
		feangio.adjust_mesh_stiffness(fem);								// Uncomment this if not using a composite material model.  This will update the stiffness of any element that contains microvessel segements 
		
		feangio.Subgrowth(feangio.sub_cycles, fem);						// Discretize the growth step over N steps in quasi-time to produce smoother results (default N = 2)
		
		if (feangio.feconv == false)									// If the FE model did not converge, end the simulation
			break;
		
		feangio.update_ECM();											// Update the ECM based on the solution from FEBio (collagen fiber orientation and matrix density)

		feangio.update_grid_volume();									// Re-calculate the volume of each elment after the deformation using the Jacobian								

		feangio.save_time();											// Output time information	
		
		feangio.fileout.printStatus(feangio.data);						// Print the status of angio3d to the user    
	}
	
	//// ANGIO3D - Generate output files
	feangio.output_params();										// Output parameters for simulation (sproutf, tip_range, phi_stiff_factor)
	
	//feangio.assemble_sprout_nodes();								// Assemble the sprout locations into a structure				
	
	//feangio.fileout.printsproutnodes(feangio.sprout_nodes);		// Output the sprout locations
	
	feangio.fileout.dataout(feangio);								// Output data file
		
	feangio.fileout.writeCollFib(feangio.grid, false);				// Output final collagen fiber orientation

	feangio.fileout.writeECMDen(feangio.grid);						// Output final matrix density

	feangio.fileout.writeSegConn(feangio.frag);						// Output the segment connectivity data
	
	feangio.fileout.writeECMDenStore(feangio.grid);

	feangio.fileout.writeECMFibrilStore(feangio.grid);

	//feangio.fileout.writeECMDenGrad(feangio.grid);				// Output the ECM density gradient, I don't think this works right...

	feangio.save_cropped_vessels();									// Output the vessels within a selected region of interest

	return true;					
}
