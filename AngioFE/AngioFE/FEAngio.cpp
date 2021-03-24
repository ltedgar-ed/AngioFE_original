///////////////////////////////////////////////////////////////////////
// FEAngio.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "FEAngio.h"
#include "Segment.h"
#include "FESproutBodyForce.h"
#include "FECore/FECoreKernel.h"
#include "FEBioMech/FEElasticMaterial.h"
#include "FEBioMech/FEElasticMixture.h"

extern FECoreKernel* pFEBio;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FEAngio::FEAngio(FEAngInput &input) : Angio(input)
{																				
	sproutf = input.isproutf;													// Read in the magnitude of the body force produced by the sprout tips from the input class
	tip_range = input.itip_range;												// Read in the range in exponetial function that produces body force produced by sprout tips (um) from the input class
	b = 1./tip_range;
	sub_cycles = input.isub_cyc;												// Read in the number of sub-growth cycles from the input class											
	
	phi_stiff_factor = input.istifffact;										// Read in the stiffness factor, which scales the amount of mesh displacement that's sent to the microvessels, from the input class
	disp_vess = true;

	yes_branching = true;														// Flag to turn branching on/off
	yes_anast = true;															// Flag to turn anastomosis on/off
	
	if (input.ibranch == 0){ yes_branching = false; data.branch_chance = 0.0;}	// Read in the branching flag from the input class
	if (input.ianast == 0) yes_anast = false;									// Read in the anastomosis flag from the input class

	// If running the sprout verification problem...
	if (input.isprout_verify == 1){
		sprout_verify = true;
		yes_branching = false;										
		yes_anast = false;
		kill_off = true;
		disp_vess = true;
		data.length_adjust = 1e-10;
		updateLength();
	}

	total_bdyf = 0;													// Body force counter
	
	FE_state = 0;													// State counter to count the number of solved FE states
	FE_time_step = 0.5;												// Time step within FEBio

	stream = fopen("out_vess_state.ang","wt");						// Open the stream for the vessel state data file		
	bf_stream = fopen("out_bf_state.ang","wt");						// Open the stream for the body force state data file
	
	time_stream = fopen("out_time.ang","wt");						// Open the stream for the time and state data file
	time_write_headers = true;										// Set the time and state data file to write the headers on its first output

	comp_mat = input.icomp_mat;										// Determine if the composite constitutive model is being used

	feconv = true;													// 

	m_pmat = 0;
}


FEAngio::~FEAngio()
{
	fclose(stream);													// Close the stream for the vessel state data file
	fclose(bf_stream);												// Close the stream for the body force state data file
	fclose(time_stream);											// Close the stream for the time and state data file
}

///////////////////////////////////////////////////////////////////////
// Member Functions:
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
// FEAngio - Growth
//      Vessel elongation is represented by the addition of new line segments at the locations of the active sprouts.
///////////////////////////////////////////////////////////////////////

void FEAngio::Growth(FEModel& fem)
{
    int k;															// Iterator for the segment tips
    list<Segment>::iterator it;										// Declare iterator through list container 'it'
    	
	//// Elongate the active vessel tips
    //#pragma omp parallel for
	for (it = frag.begin(); it != frag.end(); ++it)					// Iterate through each vessel segment...
	{
		for (k=0; k<=1; ++k)											// Iterate through each tip of the vessel segment...
		{
	        if (it->tip[k] != 0)											// If tip is active (not 0)...
			{
		        Segment seg;													// Declare SEGMENT object 'seg'
				seg = cult.createNewSeg(it,grid,data,k,frag);					// Create new vessel segment at the current tip existing segment 

				frag.push_front (seg);											// Append new segment at the top of the list 'frag'
				
				++data.nsegs;													// Iterate the total segment counter
			
			}
            
	    }
		    
		if (yes_branching == true)										// If branching is turned on...
			Branch(it, fem);												// Determine if the segment forms a new branch
		    
    }

	if (yes_anast == true)											// If anatomosis is turned on...
		Fuse();															// Determine which segments form anatomoses

	   
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - Branch
//		Branching is modeled as a random process, determine is the segment passed to this function forms a branch
///////////////////////////////////////////////////////////////////////

void FEAngio::Branch(list<Segment>::iterator it, FEModel& fem)
{
    int k;
    //list<Segment>::iterator it;                                         // Declare iterator through list container 'it'
    
    // Generate a random number between 0 and 1. If that number is less than the branching probability at
    // the current time, or initial branching is activated (data.ini_branch = true), then branching occurs
	
    double den_scale = 1.0;												// Declare the density scaling factor

	double xpt = (it->x[1] + it->x[0])/2;								// Set the branch point to be half-way along the length of the segment
	double ypt = (it->y[1] + it->y[0])/2;
	double zpt = (it->z[1] + it->z[0])/2;
	
	den_scale = cult.findDenScale(xpt, ypt, zpt, grid);					// Set the density scaling factor
	

	if ( float(rand())/RAND_MAX < den_scale*data.dt*data.branch_chance/data.t || (it->init_branch == true) )	// The segment generates a random number, if this number is less than the branching probability, or if the segment has an initial branch, the form a branch
    {
	    if ((it->BCdead == 0) && (it->anast == 0))                  // Segments that have encountered a boundary condition or formed an anastomoses may not form a branch
		{                                                           
	        
			Segment seg;                                            // Declare SEGMENT object 'seg'
			data.num_branches = data.num_branches + 1;              // Iterate the total number of branches +1
			data.branch = true;                                     // Branching flag set to 'true.' This tells the program that the
			                                                        // new vessel segment being created is arising from a branch      
    		it->init_branch = false;								// Turn off the initial branch flag
    				
			if (float(rand())/RAND_MAX < 0.5)                       // Randomly determine which node of the parent segment forms the branch
			    k = 0;
		    else
			    k = 1;
    							
			it->tip[k] = sign(0.5f - float(rand())/RAND_MAX);        // Randomly assign the branch to grow as +1 or -1
    				
			seg = cult.createNewSeg(it,grid,data,k,frag);           // Create the new vessel segment

			//if (sprout_verify == true){
			//	seg.uvect.x = 0;
			//	seg.uvect.y = 1;
			//	seg.uvect.z = 0;}

			//if (seg.tip[k] == 1)
			//		seg_nodes.push_back(vect3(seg.x[1],seg.y[1],seg.z[1]));
			//	else if (seg.tip[k] == -1)
			//		seg_nodes.push_back(vect3(seg.x[0],seg.y[0],seg.z[0]));
			
            data.num_vessel = data.num_vessel + 1;					// Iterate the vessel counter
		    seg.vessel = data.num_vessel;							// Set the segments ID number
    				                    
			it->Recent_branch = 1;
			seg.Recent_branch = 1;                                  // Indicate that this segment just branched, prevents future segments from branching again too soon
    		
			create_branching_force(seg, fem);						// Create a new sprout force for the branch

			frag.push_front (seg);                                  // Append new segment at the top of the list 'frag'
			data.branch = false;                                    // Turn off branching flag once branching algorithm is complete
        }
	}                                                              
	

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - Subgrowth
//		Step through growth through a series of sub-steps to produce smoother model results
///////////////////////////////////////////////////////////////////////

void FEAngio::Subgrowth(int sub_steps, FEModel& fem)
{
	double subgrowth_scale = 0.;									// The subgrowth scaling factor					
	//double mag;
	list<list<Segment>::iterator >::iterator tip_it;				// Iterator through the list of active segment tips
	list<Segment>::iterator frag_it;								// Iterator through the segment container

	FEMesh& mesh = fem.GetMesh();										// Obtain the FE mesh

	update_sprout_stress_scaling();

	for (int k = 1; k <= sub_steps; k++)							// Iterate through the number of substeps...
	{
		subgrowth_scale = ((double)k/(double)sub_steps);				// Update the value of the subgrowth scaling factor
		
		//if (sprout_verify == true)
		//	subgrowth_scale = 0;

		for (tip_it = active_tips.begin(); tip_it != active_tips.end(); ++tip_it)	// For each of the active growth tips...
		{
			Segment& seg = (*(*tip_it));											// Dereference the tip iterator to obtain the active segment

			// Step growth for the active segment
			if (seg.tip[0] == -1){
				seg.x[0] = seg.x[1] + subgrowth_scale*seg.length*seg.uvect.x;     
				seg.y[0] = seg.y[1] + subgrowth_scale*seg.length*seg.uvect.y;     
				seg.z[0] = seg.z[1] + subgrowth_scale*seg.length*seg.uvect.z;}

			if (seg.tip[1] == 1){
				seg.x[1] = seg.x[0] + subgrowth_scale*seg.length*seg.uvect.x;     
				seg.y[1] = seg.y[0] + subgrowth_scale*seg.length*seg.uvect.y;     
				seg.z[1] = seg.z[0] + subgrowth_scale*seg.length*seg.uvect.z;}
		}

		update_body_forces(fem, 1.0);									// Update the positions of the body forces

		FE_state++;														// Increase the FE state

		// Reset some parameters for FEBio
		FEAnalysis* pstep = fem.GetCurrentStep();
		pstep->m_dt = pstep->m_dt0;
		pstep->m_ntimesteps = 0;
		pstep->m_iteopt = 100;
		pstep->m_maxretries = 10;

		//profiler.time_in(1);
		feconv = fem.Solve();											// Solve the FE problem
		//profiler.time_out();

		if (feconv == false)											// If the solver did not converge...
			return;

		//mag = sprout_mag();
		//FELoadCurve* plc = fem.GetLoadCurve(1);							// Obtain the load curve
		//plc->SetPoint(0,((double)FE_state - 1.0)*FE_time_step,0);		// Update the load curve
		//plc->SetPoint(1,((double)FE_state)*FE_time_step,mag);

		update_grid(mesh);												// Update the growth model mesh using the FE solution
		displace_vessels();												// Update microvessel position and orientation using the displacement field from the FE solution
		
		save_vessel_state();											// Save the current vessel state
		//save_bdy_forces(fem);											// Save the current positions of the body forces (doesn't work with new method of storing sprout forces)
	}

	
	for (frag_it = frag.begin(); frag_it != frag.end(); ++frag_it)      // Iterate through all segments in frag list container (it)                               
	{
		frag_it->findlength();												// Calculate the new length of the microvessel segments										
	}
	
	removeErrors();														// Remove bad segments
	
	//error_term = false;
	
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - displace_vessels
//		Use the displacement field from the FE solution to update microvessels into the current configuration
///////////////////////////////////////////////////////////////////////

void FEAngio::displace_vessels()
{
	if (disp_vess == true){											// If microvessel displacement is turned on
		int k = 0;													// Iterator for segment tips
		int elem_num = 0;											// Element number
		list<Segment>::iterator it;									// Iterator for the segment container FRAG
		double xpt = 0.; double ypt = 0.; double zpt = 0.;			// Position in Cartesian coordinates
		double xix = 0.; double xiy = 0.; double xiz = 0.;			// Position in the element's natural coordinates		
		double xold = 0.; double yold = 0.; double zold = 0.;		// Old position in Cartesian coordinates
		double shapeF[8] = {0.};									// Array containing the shape function values at the segment's position
		vect3 disp; vect3 weighted_disp;							// Displacement and weighted displacement vectors						
		Elem elem;													// Element containing the segment
    
		for (it = frag.begin(); it != frag.end(); ++it)             // Iterate through all segments in frag list container (it)                               
		{
			for (k=0; k<=1;++k)                                         // Iterate through both segment tips (k)
			{
				xpt = it->x[k];											// Set position to the current segment tip
				ypt = it->y[k];
				zpt = it->z[k];
		    
				int BC_face = 0;
				elem_num = it->tip_elem[k];								// Find the element that contains the segment tip
		    
				if (elem_num != -1)										// If the segment is inside the mesh and has a real element number...
				{
					elem = grid.ebin[elem_num];								// Obtain the element
		    
					grid.natcoord(xix, xiy, xiz, xpt, ypt, zpt, elem_num);	// Convert the position to natural coordinates
		    
					grid.shapefunctions(shapeF, xix, xiy, xiz);				// Obtain the values of the shape functions at this position
		    
					// Calculate the displacement vector by interpolating nodal displacement to the segment tip
					disp.x = shapeF[0]*(*elem.n1).u.x + shapeF[1]*(*elem.n2).u.x + shapeF[2]*(*elem.n3).u.x + shapeF[3]*(*elem.n4).u.x + shapeF[4]*(*elem.n5).u.x + shapeF[5]*(*elem.n6).u.x + shapeF[6]*(*elem.n7).u.x + shapeF[7]*(*elem.n8).u.x;
					disp.y = shapeF[0]*(*elem.n1).u.y + shapeF[1]*(*elem.n2).u.y + shapeF[2]*(*elem.n3).u.y + shapeF[3]*(*elem.n4).u.y + shapeF[4]*(*elem.n5).u.y + shapeF[5]*(*elem.n6).u.y + shapeF[6]*(*elem.n7).u.y + shapeF[7]*(*elem.n8).u.y;
					disp.z = shapeF[0]*(*elem.n1).u.z + shapeF[1]*(*elem.n2).u.z + shapeF[2]*(*elem.n3).u.z + shapeF[3]*(*elem.n4).u.z + shapeF[4]*(*elem.n5).u.z + shapeF[5]*(*elem.n6).u.z + shapeF[6]*(*elem.n7).u.z + shapeF[7]*(*elem.n8).u.z;
		    		    
					weighted_disp = disp*phi_stiff_factor;					// Calculate the weighted displacement vector
		
					xold = it->x[k]; yold = it->y[k]; zold = it->z[k];		// Store the old position
				
					// Update the segment tip position using the new weighted displacement vector
					it->x[k] = xold + weighted_disp.x;						
					it->y[k] = yold + weighted_disp.y;
					it->z[k] = zold + weighted_disp.z;
			
					if (grid.findelem(it->x[k],it->y[k],it->z[k]) == -1){	// If using the weighted displacement vector causes the segment to move outside the mesh...
						weighted_disp = disp;									// Update using the full displacement vector instead
						it->x[k] = xold + weighted_disp.x;
						it->y[k] = yold + weighted_disp.y;
						it->z[k] = zold + weighted_disp.z;}
				}
				else													// If the segment doesn't have a real element number...
				{
					it->mark_of_death = true;								// Kill the segment
					it->death_label = 9;
				}			
			}
		
			it->findunit();											// Recalculate the segment's unit vector based on it's new position
			//it->findphi();
		}
	}

    return;
}   


///////////////////////////////////////////////////////////////////////
// FEAngio - apply_sprout_forces
//		Apply sprout forces to the mesh for each active vessel tip
///////////////////////////////////////////////////////////////////////

void FEAngio::apply_sprout_forces(FEModel& fem, int load_curve, double scale)
{
	list<list<Segment>::iterator >::iterator tip_it;				// Iterator for the container that stores active growth tips
	Segment seg;													// Segment placeholder
	vec3d sprout_vect;												// Sprout vector

	double tipx = 0.; double tipy = 0.; double tipz = 0.;			// Position of the tip
	double magnitude = scale*sproutf;								// Scale the sprout magnitude

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (data.t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (data.t < 4.0)
		magnitude = (1.0/4.0)*data.t*scale;

	//#pragma omp parallel for
	for (tip_it = active_tips.begin(); tip_it != active_tips.end(); ++tip_it)		// For each active growth tip...
	{
		seg = (*(*tip_it));												// Obtain the growth tip

		if (seg.tip[0] == -1){											// If it's a -1 tip...
			tipx = seg.x[0];												// Obtain the position of the active tip
			tipy = seg.y[0];
			tipz = seg.z[0];
			
			sprout_vect.x = seg.x[0] - seg.x[1];							// Calculate the directional unit vector of the sprout
			sprout_vect.y = seg.y[0] - seg.y[1];
			sprout_vect.z = seg.z[0] - seg.z[1];
			sprout_vect = sprout_vect/sprout_vect.norm();

			(*tip_it)->bdyf_id[0] = create_body_force(sprout_vect, tipx, tipy, tipz, magnitude, b, fem, load_curve);}				// Create a new body force, set the tips body force ID
		
		if (seg.tip[1] == 1){											// If it's a +1 tip...
			tipx = seg.x[1];												// Obtain the position of the active tip
			tipy = seg.y[1];
			tipz = seg.z[1];
			
			sprout_vect.x = seg.x[1] - seg.x[0];							// Calculate the directional unit vector of the sprout
			sprout_vect.y = seg.y[1] - seg.y[0];
			sprout_vect.z = seg.z[1] - seg.z[0];
			sprout_vect = sprout_vect/sprout_vect.norm();

			(*tip_it)->bdyf_id[1] = create_body_force(sprout_vect, tipx, tipy, tipz, magnitude, b, fem, load_curve);}				// Create a new body force, set the tips body force ID
	}

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - create_body_force
//		Add a new body force entry into the body force field applyied to the mesh
///////////////////////////////////////////////////////////////////////

int FEAngio::create_body_force(vec3d sprout_vect, double xpt, double ypt, double zpt, double mag, double range, FEModel& fem, int load_curve)
{
	total_bdyf++;							// Iterate the total body force counter							


//--> SAM
	if (m_pmat)
	{
		m_pmat->AddSprout(vec3d(xpt, ypt, zpt), sprout_vect);
		return m_pmat->Sprouts() - 1;
	}
	else
	{
//<-- SAM
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(fem.GetBodyLoad(0));					// Obtain the body force class
		pbf->AddSprout(vec3d(xpt, ypt, zpt),sprout_vect);												// Add a new component to the body force for this active sprout tip
		return pbf->Sprouts() - 1;																		// Return the ID number for the body force
	}
}


///////////////////////////////////////////////////////////////////////
// FEAngio - update_body_forces
//		Update the sprout forces after a deformation
///////////////////////////////////////////////////////////////////////

void FEAngio::update_body_forces(FEModel& fem, double scale)
{
	list<Segment>::iterator frag_it;								// Iterator for the segment container FRAG
	vec3d sprout_vect;												// Sprout direction vector

	double tipx = 0.; double tipy = 0.; double tipz = 0.;			// Position of the tip
	double magnitude = scale*sproutf;								// Magnitude of the sprout force

	// Ramp up the sprout force magnitude up to time t = 4.0 days
	if (data.t == 0.0)
		magnitude = (1.0/4.0)*0.001*scale; 
	else if (data.t < 4.0)
		magnitude = (1.0/4.0)*data.t*scale;

	if (m_pmat)
	{
		//m_pmat->scale = magnitude;
		int NSP = m_pmat->Sprouts();
		for (int i=0; i<NSP; ++i)
		{
			FEAngioMaterial::SPROUT& sp = m_pmat->GetSprout(i);
			sp.bactive = false;
		}
	}
	else
	{
		FESproutBodyForce* pbf = dynamic_cast<FESproutBodyForce*>(fem.GetBodyLoad(0));			// Obtain the sprout body force field
		int NSP = pbf->Sprouts();										// Obtain the number of sprouts
		for (int i = 0; i < NSP; i++)									// Deactivate all sprout force components
		{
			FESproutBodyForce::SPROUT& sp = pbf->GetSprout(i);				// Obtain the sprout force component
			sp.active = false;												// Deactive the sprout force component
		}
		FEParameterList& pl = pbf->GetParameterList();										// Get the body force's parameter list
		FEParam* pa = pl.Find("a"); assert(pa);												// Get the sprout force magnitude parameter
		pa->value<double>() = magnitude*sproutf;													// Set the sprout force magnitude parameter
	}

	FEMesh& mesh = fem.GetMesh();									// Obtain the FE mesh

	//#pragma omp parallel for
	for (frag_it = frag.begin(); frag_it != frag.end(); ++frag_it)		// Iterate through each segment in the model...
	{
		const Segment& seg = (*frag_it);								// Obtain the segment, keep it constant to prevent changes

		if (((seg.tip[0] == -1) || (seg.tip_BC[0] == 1)) && (seg.bdyf_id[0] >= 0)){		  // Turn on the body force for any active -1 segment OR -1 segment that has encountered a boundary and stopped growing...
		//if ((seg.tip[0] == -1) && (seg.bdyf_id[0] >= 0)){									// Turn on the body force for any active -1 segment
			tipx = seg.x[0];																	// Obtain the tip position
			tipy = seg.y[0];
			tipz = seg.z[0];

			sprout_vect.x = seg.x[0] - seg.x[1];												// Calculate the sprout directional vector
			sprout_vect.y = seg.y[0] - seg.y[1];
			sprout_vect.z = seg.z[0] - seg.z[1];
			sprout_vect = sprout_vect/sprout_vect.norm();			

//--> SAM
			update_angio_sprout(seg.bdyf_id[0], true, vec3d(tipx, tipy, tipz), sprout_vect);
//<-- SAM
			}
		
		if (((seg.tip[1] == 1) || (seg.tip_BC[1] == 1)) && (seg.bdyf_id[1] >= 0)){		  // Turn on the body force for any active +1 segment OR +1 segment that has encountered a boundary and stopped growing...
		//if ((seg.tip[1] == 1) && (seg.bdyf_id[1] >= 0)){									// Turn on the body force for any active +1 segment
			tipx = seg.x[1];																	// Obtain the tip position
			tipy = seg.y[1];
			tipz = seg.z[1];
			
			sprout_vect.x = seg.x[1] - seg.x[0];												// Calculate the sprout directional vector
			sprout_vect.y = seg.y[1] - seg.y[0];
			sprout_vect.z = seg.z[1] - seg.z[0];
			sprout_vect = sprout_vect/sprout_vect.norm();

//--> SAM
			update_angio_sprout(seg.bdyf_id[1], true, vec3d(tipx, tipy, tipz), sprout_vect);
//<-- SAM
			}
	}
	
	return;
}

//--> SAM
void FEAngio::update_angio_sprout(int id, bool bactive, const vec3d& rc, const vec3d& sprout_vect)
{
	if (m_pmat)
	{
		FEAngioMaterial::SPROUT& sp = m_pmat->GetSprout(id);
		sp.bactive = true;
		sp.sprout = sprout_vect;
		m_pmat->UpdateSprout(sp, rc);
	}
	else
	{
		FESproutBodyForce::SPROUT& sp = m_pbf->GetSprout(id);		// Obtain the sprout component 
		sp.active = true;											// Set the sprout component to active
		sp.rc = rc;													// Set the tip position
		sp.sprout = sprout_vect;									// Set the sprout force directional vector
	}
}
//<-- SAM

///////////////////////////////////////////////////////////////////////
// FEAngio - create_branching_force
//		Create a new sprout force component for a newly formed branch		
///////////////////////////////////////////////////////////////////////

void FEAngio::create_branching_force(Segment& seg, FEModel& fem)
{
	double tipx = 0.; double tipy = 0.; double tipz = 0.;						// Position of the new sprout tip 			
	vec3d sprout_vect;															// Sprout for directional vector

	total_bdyf = 0;																// Obtain the total number of sprouts
	if (m_pmat) total_bdyf = m_pmat->Sprouts();
	else total_bdyf = m_pbf->Sprouts();
	
	if (seg.tip[0] == -1){														// If the new branch is a -1 segment... 														
		tipx = seg.x[0];															// Obtain the position of the new tip
		tipy = seg.y[0];
		tipz = seg.z[0];
		
		sprout_vect.x = seg.x[0] - seg.x[1];										// Calculate the sprout directional unit vector									
		sprout_vect.y = seg.y[0] - seg.y[1];
		sprout_vect.z = seg.z[0] - seg.z[1];
		sprout_vect = sprout_vect/sprout_vect.norm();
		
		seg.bdyf_id[0] = total_bdyf - 1;}											// Assign the body force ID
		
	if (seg.tip[1] == 1){														// If the new branch is a +1 segment...
		tipx = seg.x[1];															// Obtain the position of the new tip
		tipy = seg.y[1];
		tipz = seg.z[1];
		
		sprout_vect.x = seg.x[1] - seg.x[0];										// Calculate the sprout directional unit vector									
		sprout_vect.y = seg.y[1] - seg.y[0];
		sprout_vect.z = seg.z[1] - seg.z[0];
		sprout_vect = sprout_vect/sprout_vect.norm();
		
		seg.bdyf_id[1] = total_bdyf - 1;}											// Assign the body force ID

//--> SAM
	if (m_pmat)
	{
		m_pmat->AddSprout(vec3d(tipx, tipy, tipz), sprout_vect);
		total_bdyf = m_pmat->Sprouts();
	}
	else
	{
//!<-- SAM
		m_pbf->AddSprout(vec3d(tipx, tipy, tipz), sprout_vect);						// Add the new sprout component to the sprout force field
		total_bdyf = m_pbf->Sprouts();												// Update the total number of sprouts
	}
	
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - save_vessel_state
//		Save microvessel position at the current time point
///////////////////////////////////////////////////////////////////////

void FEAngio::save_vessel_state()
{
	list<Segment>::iterator it;													// Iterator for the segment container FRAG
		
	fprintf(stream,"%-5s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n","State","Time","X1","Y1","Z1","X2","Y2","Z2","Length");  // Write column labels to out_vess_state.ang
	
	for (it = frag.begin(); it != frag.end(); ++it)								// Iterate through all segments in frag list container (it)
	{
		fprintf(stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n",FE_state,it->TofBirth,it->x[0],it->y[0],it->z[0],it->x[1],it->y[1],it->z[1],it->length);  // Write to out_vess_state.ang
	}
	
	return;
}



///////////////////////////////////////////////////////////////////////
// FEAngio - save_bdy_forces
//		Save positions of the body forces at the current time step (This function needs to be re-written)
///////////////////////////////////////////////////////////////////////

void FEAngio::save_bdy_forces(FEModel& fem)
{
	int NBF = fem.BodyLoads();

	fprintf(bf_stream,"%-5s %-12s %-12s %-12s %-12s %-12s %-12s %-12s %-12s\n","State","Time","X1","Y1","Z1","X2","Y2","Z2","Length"); 

	for (int i = 0; i < NBF; ++i)
	{
		FEBodyForce* pbf = dynamic_cast<FEBodyForce*>(fem.GetBodyLoad(i));
		FEParameterList& pl = pbf->GetParameterList();
		FEParam* pa = pl.Find("a");
		FEParam* prc = pl.Find("rc");

		if (pa && prc)
		{
			if (pa->value<double>() != 0.0)
				fprintf(bf_stream,"%-5.2i %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f\n",FE_state, data.t, prc->value<vec3d>().x, prc->value<vec3d>().y, prc->value<vec3d>().z, prc->value<vec3d>().x + 1.0, prc->value<vec3d>().y + 1.0, prc->value<vec3d>().z + 1.0, 1.73); 
		}
	}

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - save_time
//		Save the current time information			
///////////////////////////////////////////////////////////////////////

void FEAngio::save_time()
{
	if (time_write_headers == true){												// If this is the first time writing to out_time.ang
		fprintf(time_stream,"%-5s %-12s %-12s\n","State","Vess Time","FE Time");		// Print the column labels
		time_write_headers = false;}													// Turn off the headers flag
	
	fprintf(time_stream,"%-5.2i %-12.7f %-12.7f\n",FE_state, data.t, ((double)FE_state - 1.0)*FE_time_step);	// Print out the FE state, the vessel growth model time, and the FE time
	
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - update_grid
//		Update the grid after a deformation using the FE mesh
///////////////////////////////////////////////////////////////////////

void FEAngio::update_grid(FEMesh &mesh)
{
	for (int i = 0; i < grid.Nn; ++i)								//For each node in the grid...
	{
		grid.nodes[i].u.x = mesh.Node(i).m_rt.x - grid.nodes[i].x;		// Calculate the displacement vector by finding the difference between the node in the deformed FE mesh and the undeformed grid
		grid.nodes[i].u.y = mesh.Node(i).m_rt.y - grid.nodes[i].y;
		grid.nodes[i].u.z = mesh.Node(i).m_rt.z - grid.nodes[i].z;
				
		grid.nodes[i].x = mesh.Node(i).m_rt.x;							// Update the grid node to the current position of the FE mesh
		grid.nodes[i].y = mesh.Node(i).m_rt.y;
		grid.nodes[i].z = mesh.Node(i).m_rt.z;
	}
	
	return;
}



///////////////////////////////////////////////////////////////////////
// FEAngio - update_ECM
//		Update the ECM field after a deformation
///////////////////////////////////////////////////////////////////////

void FEAngio::update_ECM()
{
	Elem elem;														// Element placeholder
	mat3 F;															// Deformation gradient tensor				
	double Jacob = 0.;												// Jacobian (i.e., determinant of F)
	
	vect3 coll_fib;													// Collagen fiber vector
	double ecm_den = 0.;											// Collagen density

	double e1 = 0.; double e2 = 0.; double e3 = 0.;					// Natural coordinates of each node within the element

	for (int i = 0; i < grid.Ne; ++i)								// For each element within the mesh...
	{
		elem = grid.ebin[i];											// Obtain the element
		
		for (int j = 1; j < 9; j++)										// For each node in the element...
		{
			if (j == 1){													// Lower front left node
				ecm_den = (*elem.n1).ecm_den0;									// Set the matrix density
				coll_fib.reassign((*elem.n1).collfib.x,(*elem.n1).collfib.y,(*elem.n1).collfib.z);	// Set the collagen fiber vector
				e1 = -1.;														// Set the natural coordinates of the node
				e2 = -1.;
				e3 = -1.;}
			
			if (j == 2){													// Lower front right node
				ecm_den = (*elem.n2).ecm_den0;									// Set the matrix density
				coll_fib.reassign((*elem.n2).collfib.x,(*elem.n2).collfib.y,(*elem.n2).collfib.z);	// Set the collagen fiber vector
				e1 = 1.;														// Set the natural coordinates of the node
				e2 = -1.;
				e3 = -1.;}

			if (j == 3){													// Lower back left node
				ecm_den = (*elem.n3).ecm_den0;									// Set the matrix density
				coll_fib.reassign((*elem.n3).collfib.x,(*elem.n3).collfib.y,(*elem.n3).collfib.z);	// Set the collagen fiber vector	
				e1 = -1.;														// Set the natural coordinates of the node
				e2 = 1.;
				e3 = -1.;}

			if (j == 4){													// Lower back right node
				ecm_den = (*elem.n4).ecm_den0;									// Set the matrix density
				coll_fib.reassign((*elem.n4).collfib.x,(*elem.n4).collfib.y,(*elem.n4).collfib.z);	// Set the collagen fiber vector
				e1 = 1.;														// Set the natural coordinates of the node
				e2 = 1.;
				e3 = -1.;}

			if (j == 5){													// Upper front left node
				ecm_den = (*elem.n5).ecm_den0;									// Set the matrix density
				coll_fib.reassign((*elem.n5).collfib.x,(*elem.n5).collfib.y,(*elem.n5).collfib.z);	// Set the collagen fiber vector
				e1 = -1.;														// Set the natural coordinates of the node
				e2 = -1.;
				e3 = 1.;}

			if (j == 6){													// Upper front right node
				ecm_den = (*elem.n6).ecm_den0;									// Set the matrix density
				coll_fib.reassign((*elem.n6).collfib.x,(*elem.n6).collfib.y,(*elem.n6).collfib.z);	// Set the collagen fiber vector
				e1 = 1.;														// Set the natural coordinates of the node
				e2 = -1.;
				e3 = 1.;}

			if (j == 7){													// Upper back left node
				ecm_den = (*elem.n7).ecm_den0;									// Set the matrix density
				coll_fib.reassign((*elem.n7).collfib.x,(*elem.n7).collfib.y,(*elem.n7).collfib.z);	// Set the collagen fiber vector
				e1 = -1.;														// Set the natural coordinates of the node
				e2 = 1.;
				e3 = 1.;}

			if (j == 8){													// Upper back right node
				ecm_den = (*elem.n8).ecm_den0;									// Set the matrix density
				coll_fib.reassign((*elem.n8).collfib.x,(*elem.n8).collfib.y,(*elem.n8).collfib.z);	// Set the collagen fiber vector
				e1 = 1.;														// Set the natural coordinates of the node
				e2 = 1.;
				e3 = 1.;}
			
			F = calculate_deform_tensor(elem, e1, e2, e3);					// Calculate the deformation gradient tensor
			Jacob = F.det();												// Calculate the Jacobian by taking the determinant of F
			
			coll_fib = F*coll_fib;											// Update the collagen fiber orientation vector into the current configuration using F		

			if (coll_fib.norm() != 0.)										// Normalize the collagen fiber vector to obtain the unit vector
				coll_fib = coll_fib/coll_fib.norm();						

			if (Jacob != 0.)												// Update matrix density using the Jacobian
				ecm_den = ecm_den/Jacob;		
			
			if (j == 1){													// Lower front left node
				if ((*elem.n1).updated == false){								// If the node hasn't been updated yet...
					(*elem.n1).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n1).collfib.y = coll_fib.y;
					(*elem.n1).collfib.z = coll_fib.z;
					(*elem.n1).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n1).updated = true;}										// Set the updated flag
				else if ((*elem.n1).updated == true){							// If the node has been updated...
					(*elem.n1).collfib.x = ((*elem.n1).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n1).collfib.y = ((*elem.n1).collfib.y + coll_fib.y)/2;
					(*elem.n1).collfib.z = ((*elem.n1).collfib.z + coll_fib.z)/2;				
					(*elem.n1).ecm_den = ((*elem.n1).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 2){													// Lower front right node
				if ((*elem.n2).updated == false){								// If the node hasn't been updated yet...
					(*elem.n2).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n2).collfib.y = coll_fib.y;
					(*elem.n2).collfib.z = coll_fib.z;					
					(*elem.n2).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n2).updated = true;}										// Set the updated flag
				else if ((*elem.n2).updated == true){							// If the node has been updated...
					(*elem.n2).collfib.x = ((*elem.n2).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n2).collfib.y = ((*elem.n2).collfib.y + coll_fib.y)/2;
					(*elem.n2).collfib.z = ((*elem.n2).collfib.z + coll_fib.z)/2;
					(*elem.n2).ecm_den = ((*elem.n2).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 3){													// Lower back left node
				if ((*elem.n3).updated == false){								// If the node hasn't been updated yet...
					(*elem.n3).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n3).collfib.y = coll_fib.y;
					(*elem.n3).collfib.z = coll_fib.z;						
					(*elem.n3).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n3).updated = true;}										// Set the updated flag
				else if ((*elem.n3).updated == true){							// If the node has been updated...
					(*elem.n3).collfib.x = ((*elem.n3).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n3).collfib.y = ((*elem.n3).collfib.y + coll_fib.y)/2;
					(*elem.n3).collfib.z = ((*elem.n3).collfib.z + coll_fib.z)/2;
					(*elem.n3).ecm_den = ((*elem.n3).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 4){													// Lower back right node
				if ((*elem.n4).updated == false){								// If the node hasn't been updated yet...
					(*elem.n4).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n4).collfib.y = coll_fib.y;
					(*elem.n4).collfib.z = coll_fib.z;						
					(*elem.n4).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n4).updated = true;}										// Set the updated flag
				else if ((*elem.n4).updated == true){							// If the node has been updated...
					(*elem.n4).collfib.x = ((*elem.n4).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n4).collfib.y = ((*elem.n4).collfib.y + coll_fib.y)/2;
					(*elem.n4).collfib.z = ((*elem.n4).collfib.z + coll_fib.z)/2;
					(*elem.n4).ecm_den = ((*elem.n4).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 5){													// Upper front left node
				if ((*elem.n5).updated == false){								// If the node hasn't been updated yet...
					(*elem.n5).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n5).collfib.y = coll_fib.y;
					(*elem.n5).collfib.z = coll_fib.z;
					(*elem.n5).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n5).updated = true;}										// Set the updated flag
				else if ((*elem.n5).updated == true){							// If the node has been updated...
					(*elem.n5).collfib.x = ((*elem.n5).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n5).collfib.y = ((*elem.n5).collfib.y + coll_fib.y)/2;
					(*elem.n5).collfib.z = ((*elem.n5).collfib.z + coll_fib.z)/2;
					(*elem.n5).ecm_den = ((*elem.n5).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
		
			if (j == 6){													// Upper front right node
				if ((*elem.n6).updated == false){								// If the node hasn't been updated yet...
					(*elem.n6).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n6).collfib.y = coll_fib.y;
					(*elem.n6).collfib.z = coll_fib.z;
					(*elem.n6).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n6).updated = true;}										// Set the updated flag
				else if ((*elem.n6).updated == true){							// If the node has been updated...
					(*elem.n6).collfib.x = ((*elem.n6).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n6).collfib.y = ((*elem.n6).collfib.y + coll_fib.y)/2;
					(*elem.n6).collfib.z = ((*elem.n6).collfib.z + coll_fib.z)/2;
					(*elem.n6).ecm_den = ((*elem.n6).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 7){													// Upper back left node
				if ((*elem.n7).updated == false){								// If the node hasn't been updated yet...
					(*elem.n7).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n7).collfib.y = coll_fib.y;
					(*elem.n7).collfib.z = coll_fib.z;
					(*elem.n7).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n7).updated = true;}										// Set the updated flag
				else if ((*elem.n7).updated == true){							// If the node has been updated...
					(*elem.n7).collfib.x = ((*elem.n7).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n7).collfib.y = ((*elem.n7).collfib.y + coll_fib.y)/2;
					(*elem.n7).collfib.z = ((*elem.n7).collfib.z + coll_fib.z)/2;
					(*elem.n7).ecm_den = ((*elem.n7).ecm_den + ecm_den)/2;}}		// Average together the new matrix density
			
			if (j == 8){													// Upper back right node
				if ((*elem.n8).updated == false){								// If the node hasn't been updated yet...
					(*elem.n8).collfib.x = coll_fib.x;								// Set the new collagen fiber orientation
					(*elem.n8).collfib.y = coll_fib.y;
					(*elem.n8).collfib.z = coll_fib.z;
					(*elem.n8).ecm_den = ecm_den;									// Set the new matrix density
					(*elem.n8).updated = true;}										// Set the updated flag
				else if ((*elem.n8).updated == true){							// If the node has been updated...
					(*elem.n8).collfib.x = ((*elem.n8).collfib.x + coll_fib.x)/2;	// Average together the new fiber orientation vector
					(*elem.n8).collfib.y = ((*elem.n8).collfib.y + coll_fib.y)/2;
					(*elem.n8).collfib.z = ((*elem.n8).collfib.z + coll_fib.z)/2;
					(*elem.n8).ecm_den = ((*elem.n8).ecm_den + ecm_den)/2;}}		// Average together the new matrix density

		}
	}

	
	for (int i = 0; i < grid.Nn; ++i){								// Turn off the updated flag for all nodes
		grid.nodes[i].updated = false;
		grid.nodes[i].ecm_den_store.push_back(grid.nodes[i].ecm_den);
		grid.nodes[i].ecm_fibril_store.push_back(grid.nodes[i].collfib);}
	
	//update_ecm_den_grad();										  // Update the ECM density gradient based on the solution from FEBio
		
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - calculate_deform_tensor
//		Calculate the deformation gradient tensor
///////////////////////////////////////////////////////////////////////

mat3 FEAngio::calculate_deform_tensor(Elem elem, const double ex, const double ey, const double ez)
{
	mat3 F;															// Deformation gradient tensor
	mat3 dXde;														// The tensor dX/de (derivative of reference position with respect to natural coordinates) 
	
	// Calculate the derviative of the shape functions evaluate at each node
	vect3 dN1; vect3 dN2; vect3 dN3; vect3 dN4; vect3 dN5; vect3 dN6; vect3 dN7; vect3 dN8;	
   	dN1 = shapefun_d1(ex, ey, ez, 1);								
    dN2 = shapefun_d1(ex, ey, ez, 2);
	dN3 = shapefun_d1(ex, ey, ez, 3);
	dN4 = shapefun_d1(ex, ey, ez, 4);
	dN5 = shapefun_d1(ex, ey, ez, 5);
	dN6 = shapefun_d1(ex, ey, ez, 6);
	dN7 = shapefun_d1(ex, ey, ez, 7);
	dN8 = shapefun_d1(ex, ey, ez, 8);

	// Calculate dX/de
	dXde.M11 = ((*elem.n1).x0)*dN1.x + ((*elem.n2).x0)*dN2.x + ((*elem.n3).x0)*dN3.x + ((*elem.n4).x0)*dN4.x + ((*elem.n5).x0)*dN5.x + ((*elem.n6).x0)*dN6.x + ((*elem.n7).x0)*dN7.x + ((*elem.n8).x0)*dN8.x;
	dXde.M12 = ((*elem.n1).x0)*dN1.y + ((*elem.n2).x0)*dN2.y + ((*elem.n3).x0)*dN3.y + ((*elem.n4).x0)*dN4.y + ((*elem.n5).x0)*dN5.y + ((*elem.n6).x0)*dN6.y + ((*elem.n7).x0)*dN7.y + ((*elem.n8).x0)*dN8.y;
	dXde.M13 = ((*elem.n1).x0)*dN1.z + ((*elem.n2).x0)*dN2.z + ((*elem.n3).x0)*dN3.z + ((*elem.n4).x0)*dN4.z + ((*elem.n5).x0)*dN5.z + ((*elem.n6).x0)*dN6.z + ((*elem.n7).x0)*dN7.z + ((*elem.n8).x0)*dN8.z;

	dXde.M21 = ((*elem.n1).y0)*dN1.x + ((*elem.n2).y0)*dN2.x + ((*elem.n3).y0)*dN3.x + ((*elem.n4).y0)*dN4.x + ((*elem.n5).y0)*dN5.x + ((*elem.n6).y0)*dN6.x + ((*elem.n7).y0)*dN7.x + ((*elem.n8).y0)*dN8.x;
	dXde.M22 = ((*elem.n1).y0)*dN1.y + ((*elem.n2).y0)*dN2.y + ((*elem.n3).y0)*dN3.y + ((*elem.n4).y0)*dN4.y + ((*elem.n5).y0)*dN5.y + ((*elem.n6).y0)*dN6.y + ((*elem.n7).y0)*dN7.y + ((*elem.n8).y0)*dN8.y;
	dXde.M23 = ((*elem.n1).y0)*dN1.z + ((*elem.n2).y0)*dN2.z + ((*elem.n3).y0)*dN3.z + ((*elem.n4).y0)*dN4.z + ((*elem.n5).y0)*dN5.z + ((*elem.n6).y0)*dN6.z + ((*elem.n7).y0)*dN7.z + ((*elem.n8).y0)*dN8.z;

	dXde.M31 = ((*elem.n1).z0)*dN1.x + ((*elem.n2).z0)*dN2.x + ((*elem.n3).z0)*dN3.x + ((*elem.n4).z0)*dN4.x + ((*elem.n5).z0)*dN5.x + ((*elem.n6).z0)*dN6.x + ((*elem.n7).z0)*dN7.x + ((*elem.n8).z0)*dN8.x;
	dXde.M32 = ((*elem.n1).z0)*dN1.y + ((*elem.n2).z0)*dN2.y + ((*elem.n3).z0)*dN3.y + ((*elem.n4).z0)*dN4.y + ((*elem.n5).z0)*dN5.y + ((*elem.n6).z0)*dN6.y + ((*elem.n7).z0)*dN7.y + ((*elem.n8).z0)*dN8.y;
	dXde.M33 = ((*elem.n1).z0)*dN1.z + ((*elem.n2).z0)*dN2.z + ((*elem.n3).z0)*dN3.z + ((*elem.n4).z0)*dN4.z + ((*elem.n5).z0)*dN5.z + ((*elem.n6).z0)*dN6.z + ((*elem.n7).z0)*dN7.z + ((*elem.n8).z0)*dN8.z;

	// Calculate the tensor dM, which is (dX/de)^-T * dN
	vect3 dM1; vect3 dM2; vect3 dM3; vect3 dM4; vect3 dM5; vect3 dM6; vect3 dM7; vect3 dM8;
	
	mat3 dXde_inv_trans;
	dXde_inv_trans = (dXde.invert()).transpose();

	dM1.x = dXde_inv_trans.M11*dN1.x + dXde_inv_trans.M12*dN1.y + dXde_inv_trans.M13*dN1.z;
	dM1.y = dXde_inv_trans.M21*dN1.x + dXde_inv_trans.M22*dN1.y + dXde_inv_trans.M23*dN1.z;
	dM1.z = dXde_inv_trans.M31*dN1.x + dXde_inv_trans.M32*dN1.y + dXde_inv_trans.M33*dN1.z;
	
	dM2.x = dXde_inv_trans.M11*dN2.x + dXde_inv_trans.M12*dN2.y + dXde_inv_trans.M13*dN2.z;
	dM2.y = dXde_inv_trans.M21*dN2.x + dXde_inv_trans.M22*dN2.y + dXde_inv_trans.M23*dN2.z;
	dM2.z = dXde_inv_trans.M31*dN2.x + dXde_inv_trans.M32*dN2.y + dXde_inv_trans.M33*dN2.z;

	dM3.x = dXde_inv_trans.M11*dN3.x + dXde_inv_trans.M12*dN3.y + dXde_inv_trans.M13*dN3.z;
	dM3.y = dXde_inv_trans.M21*dN3.x + dXde_inv_trans.M22*dN3.y + dXde_inv_trans.M23*dN3.z;
	dM3.z = dXde_inv_trans.M31*dN3.x + dXde_inv_trans.M32*dN3.y + dXde_inv_trans.M33*dN3.z;

	dM4.x = dXde_inv_trans.M11*dN4.x + dXde_inv_trans.M12*dN4.y + dXde_inv_trans.M13*dN4.z;
	dM4.y = dXde_inv_trans.M21*dN4.x + dXde_inv_trans.M22*dN4.y + dXde_inv_trans.M23*dN4.z;
	dM4.z = dXde_inv_trans.M31*dN4.x + dXde_inv_trans.M32*dN4.y + dXde_inv_trans.M33*dN4.z;

	dM5.x = dXde_inv_trans.M11*dN5.x + dXde_inv_trans.M12*dN5.y + dXde_inv_trans.M13*dN5.z;
	dM5.y = dXde_inv_trans.M21*dN5.x + dXde_inv_trans.M22*dN5.y + dXde_inv_trans.M23*dN5.z;
	dM5.z = dXde_inv_trans.M31*dN5.x + dXde_inv_trans.M32*dN5.y + dXde_inv_trans.M33*dN5.z;
	
	dM6.x = dXde_inv_trans.M11*dN6.x + dXde_inv_trans.M12*dN6.y + dXde_inv_trans.M13*dN6.z;
	dM6.y = dXde_inv_trans.M21*dN6.x + dXde_inv_trans.M22*dN6.y + dXde_inv_trans.M23*dN6.z;
	dM6.z = dXde_inv_trans.M31*dN6.x + dXde_inv_trans.M32*dN6.y + dXde_inv_trans.M33*dN6.z;

	dM7.x = dXde_inv_trans.M11*dN7.x + dXde_inv_trans.M12*dN7.y + dXde_inv_trans.M13*dN7.z;
	dM7.y = dXde_inv_trans.M21*dN7.x + dXde_inv_trans.M22*dN7.y + dXde_inv_trans.M23*dN7.z;
	dM7.z = dXde_inv_trans.M31*dN7.x + dXde_inv_trans.M32*dN7.y + dXde_inv_trans.M33*dN7.z;

	dM8.x = dXde_inv_trans.M11*dN8.x + dXde_inv_trans.M12*dN8.y + dXde_inv_trans.M13*dN8.z;
	dM8.y = dXde_inv_trans.M21*dN8.x + dXde_inv_trans.M22*dN8.y + dXde_inv_trans.M23*dN8.z;
	dM8.z = dXde_inv_trans.M31*dN8.x + dXde_inv_trans.M32*dN8.y + dXde_inv_trans.M33*dN8.z;

	// Calculate F 
	F.M11 = ((*elem.n1).x)*dM1.x + ((*elem.n2).x)*dM2.x + ((*elem.n3).x)*dM3.x + ((*elem.n4).x)*dM4.x + ((*elem.n5).x)*dM5.x + ((*elem.n6).x)*dM6.x + ((*elem.n7).x)*dM7.x + ((*elem.n8).x)*dM8.x;
	F.M12 = ((*elem.n1).x)*dM1.y + ((*elem.n2).x)*dM2.y + ((*elem.n3).x)*dM3.y + ((*elem.n4).x)*dM4.y + ((*elem.n5).x)*dM5.y + ((*elem.n6).x)*dM6.y + ((*elem.n7).x)*dM7.y + ((*elem.n8).x)*dM8.y;
	F.M13 = ((*elem.n1).x)*dM1.z + ((*elem.n2).x)*dM2.z + ((*elem.n3).x)*dM3.z + ((*elem.n4).x)*dM4.z + ((*elem.n5).x)*dM5.z + ((*elem.n6).x)*dM6.z + ((*elem.n7).x)*dM7.z + ((*elem.n8).x)*dM8.z;

	F.M21 = ((*elem.n1).y)*dM1.x + ((*elem.n2).y)*dM2.x + ((*elem.n3).y)*dM3.x + ((*elem.n4).y)*dM4.x + ((*elem.n5).y)*dM5.x + ((*elem.n6).y)*dM6.x + ((*elem.n7).y)*dM7.x + ((*elem.n8).y)*dM8.x;
	F.M22 = ((*elem.n1).y)*dM1.y + ((*elem.n2).y)*dM2.y + ((*elem.n3).y)*dM3.y + ((*elem.n4).y)*dM4.y + ((*elem.n5).y)*dM5.y + ((*elem.n6).y)*dM6.y + ((*elem.n7).y)*dM7.y + ((*elem.n8).y)*dM8.y;
	F.M23 = ((*elem.n1).y)*dM1.z + ((*elem.n2).y)*dM2.z + ((*elem.n3).y)*dM3.z + ((*elem.n4).y)*dM4.z + ((*elem.n5).y)*dM5.z + ((*elem.n6).y)*dM6.z + ((*elem.n7).y)*dM7.z + ((*elem.n8).y)*dM8.z;

	F.M31 = ((*elem.n1).z)*dM1.x + ((*elem.n2).z)*dM2.x + ((*elem.n3).z)*dM3.x + ((*elem.n4).z)*dM4.x + ((*elem.n5).z)*dM5.x + ((*elem.n6).z)*dM6.x + ((*elem.n7).z)*dM7.x + ((*elem.n8).z)*dM8.x;
	F.M32 = ((*elem.n1).z)*dM1.y + ((*elem.n2).z)*dM2.y + ((*elem.n3).z)*dM3.y + ((*elem.n4).z)*dM4.y + ((*elem.n5).z)*dM5.y + ((*elem.n6).z)*dM6.y + ((*elem.n7).z)*dM7.y + ((*elem.n8).z)*dM8.y;
	F.M33 = ((*elem.n1).z)*dM1.z + ((*elem.n2).z)*dM2.z + ((*elem.n3).z)*dM3.z + ((*elem.n4).z)*dM4.z + ((*elem.n5).z)*dM5.z + ((*elem.n6).z)*dM6.z + ((*elem.n7).z)*dM7.z + ((*elem.n8).z)*dM8.z;

	return F;
}



///////////////////////////////////////////////////////////////////////
// FEAngio - shapefun_d1
//		Evaluate the gradient of the shape functions
///////////////////////////////////////////////////////////////////////

vect3 FEAngio::shapefun_d1(const double xix, const double xiy, const double xiz, int node)
{
    vect3 out;														// Output vector
	
	if (node == 1)													// dN1/de
    {
        out.x = -(1 - xiy)*(1 - xiz)/8.;
        out.y = -(1 - xix)*(1 - xiz)/8.;
		out.z = -(1 - xix)*(1 - xiy)/8.;
    }
    
    if (node == 2)													// dN2/de
    {
        out.x = (1 - xiy)*(1 - xiz)/8.;
        out.y = -(1 + xix)*(1 - xiz)/8.;
        out.z = -(1 + xix)*(1 - xiy)/8.;
    }    
    
    if (node == 3)													// dN3/de
    {
        out.x = -(1 + xiy)*(1 - xiz)/8.;
        out.y = (1 - xix)*(1 - xiz)/8.;
        out.z = -(1 - xix)*(1 + xiy)/8.;
    }    
    
    if (node == 4)													// dN4/de
    {
        out.x = (1 + xiy)*(1 - xiz)/8.;
        out.y = (1 + xix)*(1 - xiz)/8.;
        out.z = -(1 + xix)*(1 + xiy)/8.;
    }    
    
    if (node == 5)													// dN5/de
    {
        out.x = -(1 - xiy)*(1 + xiz)/8.;
        out.y = -(1 - xix)*(1 + xiz)/8.;
        out.z = (1 - xix)*(1 - xiy)/8.;
    }    

    if (node == 6)													// dN6/de
    {
        out.x = (1 - xiy)*(1 + xiz)/8.;
        out.y = -(1 + xix)*(1 + xiz)/8.;
        out.z = (1 + xix)*(1 - xiy)/8.;
    }    
        
    if (node == 7)													// dN7/de
    {
        out.x = -(1 + xiy)*(1 + xiz)/8.;
        out.y = (1 - xix)*(1 + xiz)/8.;
        out.z = (1 - xix)*(1 + xiy)/8.;
    }   

    if (node == 8)													// dN8/de
    {
        out.x = (1 + xiy)*(1 + xiz)/8.;
        out.y = (1 + xix)*(1 + xiz)/8.;
        out.z = (1 + xix)*(1 + xiy)/8.;
    }   

    return out;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - adjust_mesh_stiffness
//		Adjust the stiffness of the mesh based on the microvessel population
///////////////////////////////////////////////////////////////////////

void FEAngio::adjust_mesh_stiffness(FEModel& fem)
{
	if (comp_mat == 0)													// If a composite consitutive model isn't being used, exit
		return;
		
	Elem elem;															// Element placeholder
	int elem_num = 0;													// Element number
	list<Segment>::iterator frag_it;									// Iterator for the segment container FRAG
	vect3 vess_vect;													// Vessel vector

	for (int i = 0; i < grid.Ne; ++i){									// For each element within the mesh...
		grid.ebin[i].alpha = 0.;											// Set the vessel volume fraction, alpha, to zero
		grid.ebin[i].fiber_orient.x = 0.;									// Set the vessel orientation vector to 0 (this is the element fiber_orient vector, which does not contain collagen fiber orientation information but rather the material fiber direction for a transversely isotropic material model)
		grid.ebin[i].fiber_orient.y = 0.;
		grid.ebin[i].fiber_orient.z = 0.;}
	
	int Nsub = 2;														// Number of subdivisions used to calculate vessel orientation
	double sub_scale = 1/(double)Nsub;									// Calculate the subdivision scale 
	double mid_x = 0.; double mid_y = 0.; double mid_z = 0.;			// Mid-point position of each subdivision
	double elem_volume = 0.;											// Element volume
	double subunit_volume = 0.;											// Subdivision volume
	double volume_fraction = 0.;										// Volume fraction

	for (frag_it = frag.begin(); frag_it != frag.end(); ++frag_it)		// For each segment...
	{
		Segment seg;													// Segment placeholder
		Segment subunit;												// Segment subdivision placeholder
	
		seg = (*frag_it);												// Obtain the segment
		
		for (int k = 1; k <= Nsub; k++)									// For each subdivision...
		{
			if (k == 1){													// If this is the first subdivision, find the origin of the segment
				if (seg.length > 0.){											// If it's a +1 segment...
					subunit.x[0] = seg.x[0];										
					subunit.y[0] = seg.y[0];
					subunit.z[0] = seg.z[0];}
				else{															// If it's a -1 segment...
					subunit.x[1] = seg.x[1];
					subunit.y[1] = seg.y[1];
					subunit.z[1] = seg.z[1];
					subunit.length = -1.;}}
			
			// Calculate the subdivision
			if (seg.length > 0.){										// If it's a +1 segment...			
				subunit.x[1] = subunit.x[0] + sub_scale*seg.length*seg.uvect.x;     
				subunit.y[1] = subunit.y[0] + sub_scale*seg.length*seg.uvect.y;     
				subunit.z[1] = subunit.z[0] + sub_scale*seg.length*seg.uvect.z;}
			else{														// If it's a -1 segment...
				subunit.x[0] = subunit.x[1] + sub_scale*seg.length*seg.uvect.x;     
				subunit.y[0] = subunit.y[1] + sub_scale*seg.length*seg.uvect.y;     
				subunit.z[0] = subunit.z[1] + sub_scale*seg.length*seg.uvect.z;}

			subunit.findlength();										// Find the length of the subdivision
			
			mid_x = (subunit.x[1] + subunit.x[0])/2;					// Find the midpoint of the subdivision
			mid_y = (subunit.y[1] + subunit.y[0])/2;
			mid_z = (subunit.z[1] + subunit.z[0])/2;

			elem_num = grid.findelem(mid_x, mid_y, mid_z);				// Find the element that the midpoint is within

			// Calculate the orientation of the subdivision
			if (seg.length > 0.){										// If it's a +1 segment...
				vess_vect.x = subunit.x[1] - subunit.x[0];
				vess_vect.y = subunit.y[1] - subunit.y[0];
				vess_vect.z = subunit.z[1] - subunit.z[0];}
			else{														// If it's a -1 segment...
				vess_vect.x = subunit.x[0] - subunit.x[1];
				vess_vect.y = subunit.y[0] - subunit.y[1];
				vess_vect.z = subunit.z[0] - subunit.z[1];}

			if (vess_vect.norm() != 0.)									// Normalize to find the unit vector
				vess_vect = vess_vect/vess_vect.norm();

			if (elem_num != -1)											// If the midpoint has a real element number...
			{
				elem_volume = grid.ebin[elem_num].volume;					// Calculate the volume of the element
				subunit_volume = pi*(data.vessel_width/2.)*(data.vessel_width/2.)*fabs(subunit.length);		// Find the volume of the subdivision
				volume_fraction = subunit_volume/elem_volume;				// Calculate the volume fraction

				grid.ebin[elem_num].alpha = grid.ebin[elem_num].alpha + volume_fraction;	// Add the volume fraction for each subdivision to alpha
			
				// Calculate the vessel orientation vector 
				if ((grid.ebin[elem_num].fiber_orient.x == 0) && (grid.ebin[elem_num].fiber_orient.y == 0) && (grid.ebin[elem_num].fiber_orient.z == 0)){	// If the vessel orientation vector hasn't been assigned yet...
					grid.ebin[elem_num].fiber_orient.x = vess_vect.x;			// Set the vessel orientation vector					
					grid.ebin[elem_num].fiber_orient.y = vess_vect.y;
					grid.ebin[elem_num].fiber_orient.z = vess_vect.z;}
				else{														// If it has been...	
					grid.ebin[elem_num].fiber_orient.x = (grid.ebin[elem_num].fiber_orient.x + vess_vect.x)/2;	// Average together the vessel orientation vector
					grid.ebin[elem_num].fiber_orient.y = (grid.ebin[elem_num].fiber_orient.y + vess_vect.y)/2;
					grid.ebin[elem_num].fiber_orient.z = (grid.ebin[elem_num].fiber_orient.z + vess_vect.z)/2;}
			}
			
			// Set the origin of the next subdivision to the end of the current one
			if (seg.length > 0.){
				subunit.x[0] = subunit.x[1];
				subunit.y[0] = subunit.y[1];
				subunit.z[0] = subunit.z[1];}
			else{
				subunit.x[1] = subunit.x[0];
				subunit.y[1] = subunit.y[0];
				subunit.z[1] = subunit.z[0];}
		}
	}
	
	double alpha = 0.;													// Volume fraction for the composite material model
	vect3 e1; vect3 e2; vect3 e3;										// Basis for the material coordinate system (e1 is the material direction, e2 and e3 are isotropic)
	
	FEMesh& mesh = fem.GetMesh();										// Get the FE mesh
	int J = mesh.Domains();												// Find the number of domains within the mesh
    int num_elem = 0;

	for (int k = 0; k < J; ++k){
		FEDomain& d = mesh.Domain(k);										// Obtain the domain
		
		for (int j = 0; j < d.Elements(); ++j)								// For each element within the domain...
		{
			FEElement& e = d.ElementRef(j);										// Obtain the element from the domain
			int nint = e.GaussPoints();											// Obtain the number of gauss points
		
			alpha = grid.ebin[num_elem].alpha;											// Obtain alpha from the grid element

			e1.x = grid.ebin[num_elem].fiber_orient.x;									// Set e1 to the vessel orientation vector
			e1.y = grid.ebin[num_elem].fiber_orient.y;
			e1.z = grid.ebin[num_elem].fiber_orient.z;

			if ((e1.x == 0) && (e1.y == 0) && (e1.z == 0)){						// If there is not vessels in the element, set the material basis to the global coordinate basis
				e1.reassign(1,0,0);
				e2.reassign(0,1,0);
				e3.reassign(0,0,1);}
			else{																// Else, set the other two directions to be orthogonal to the vessel orientation
				e2.y = 1;
				e2 = e1^e2;
				e3 = e1^e2;}
		
			for (int n = 0; n < nint; ++n)										// For each gauss point...
			{
				//FEMaterialPoint& mp = *e.m_State[n];								// Obtain the material point
				//FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();	// Obtain the mixture material point
				//pt.m_w[0] = alpha;													// Set the first weight factor to alpha
				//pt.m_w[1] = 1.0 - alpha;											// Set the second weight factor to 1 - alpha
			
				FEMaterialPoint& mp = *(e.GetMaterialPoint(n)->GetPointData(1));    // returns the second component of the mixture
				FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>(); // get the mixture material point
				pt.m_w[0] = alpha;
				pt.m_w[1] = 1.0 - alpha;

				if (comp_mat == 2){													// If the transversely isotropic material is being used...
					FEElasticMaterialPoint& pt2 = *mp.ExtractData<FEElasticMaterialPoint>();
					pt2.m_Q[0][0] = e1.x;												// Set the first column of Q to e1
					pt2.m_Q[1][0] = e1.y;
					pt2.m_Q[2][0] = e1.z;
					pt2.m_Q[0][1] = e2.x;												// Set the second column of Q to e2
					pt2.m_Q[1][1] = e2.y;
					pt2.m_Q[2][1] = e2.z;
					pt2.m_Q[0][2] = e3.x;												// Set the third column of Q to e3
					pt2.m_Q[1][2] = e3.y;
					pt2.m_Q[2][2] = e3.z;}
				}

			num_elem++;
		}
	}
	//FEDomain& d = mesh.Domain(0);										// Obtain the domain

 //   for (int j = 0; j < d.Elements(); ++j)								// For each element within the domain...
	//{
	//	FEElement& e = d.ElementRef(j);										// Obtain the element from the domain
	//	int nint = e.GaussPoints();											// Obtain the number of gauss points
	//	
	//	alpha = grid.ebin[j].alpha;											// Obtain alpha from the grid element

	//	e1.x = grid.ebin[j].fiber_orient.x;									// Set e1 to the vessel orientation vector
	//	e1.y = grid.ebin[j].fiber_orient.y;
	//	e1.z = grid.ebin[j].fiber_orient.z;

	//	if ((e1.x == 0) && (e1.y == 0) && (e1.z == 0)){						// If there is not vessels in the element, set the material basis to the global coordinate basis
	//		e1.reassign(1,0,0);
	//		e2.reassign(0,1,0);
	//		e3.reassign(0,0,1);}
	//	else{																// Else, set the other two directions to be orthogonal to the vessel orientation
	//		e2.y = 1;
	//		e2 = e1^e2;
	//		e3 = e1^e2;}
	//	
	//	for (int n = 0; n < nint; ++n)										// For each gauss point...
	//	{
	//		//FEMaterialPoint& mp = *e.m_State[n];								// Obtain the material point
	//		//FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>();	// Obtain the mixture material point
	//		//pt.m_w[0] = alpha;													// Set the first weight factor to alpha
	//		//pt.m_w[1] = 1.0 - alpha;											// Set the second weight factor to 1 - alpha
	//		
	//		FEMaterialPoint& mp = *(e.GetMaterialPoint(n)->GetPointData(1));    // returns the second component of the mixture
	//		FEElasticMixtureMaterialPoint& pt = *mp.ExtractData<FEElasticMixtureMaterialPoint>(); // get the mixture material point
	//		pt.m_w[0] = alpha;
	//		pt.m_w[1] = 1.0 - alpha;

	//		if (comp_mat == 2){													// If the transversely isotropic material is being used...
	//			FEElasticMaterialPoint& pt2 = *mp.ExtractData<FEElasticMaterialPoint>();
	//			pt2.m_Q[0][0] = e1.x;												// Set the first column of Q to e1
	//			pt2.m_Q[1][0] = e1.y;
	//			pt2.m_Q[2][0] = e1.z;
	//			pt2.m_Q[0][1] = e2.x;												// Set the second column of Q to e2
	//			pt2.m_Q[1][1] = e2.y;
	//			pt2.m_Q[2][1] = e2.z;
	//			pt2.m_Q[0][2] = e3.x;												// Set the third column of Q to e3
	//			pt2.m_Q[1][2] = e3.y;
	//			pt2.m_Q[2][2] = e3.z;}
	//		}
	//}
		
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - initialize_grid_volume
//		Set the initial grid volume
///////////////////////////////////////////////////////////////////////

void FEAngio::initialize_grid_volume()
{
	mat3 Jacob_mat;														// Jacobian matrix
	Elem elem;															// Element placeholder
	double ex = 0.; double ey = 0.; double ez = 0.;						// Position of the centroid in natural coordinates
	vect3 dN1; vect3 dN2; vect3 dN3; vect3 dN4; vect3 dN5; vect3 dN6; vect3 dN7; vect3 dN8;	// Arrays containing the derivatives of the shape functions
	double volume = 0.;													// Element volume

	// Calculate the value of the derviatives of the shape functions evaulated at the centroid
	dN1 = shapefun_d1(ex, ey, ez, 1);
    dN2 = shapefun_d1(ex, ey, ez, 2);
	dN3 = shapefun_d1(ex, ey, ez, 3);
	dN4 = shapefun_d1(ex, ey, ez, 4);
	dN5 = shapefun_d1(ex, ey, ez, 5);
	dN6 = shapefun_d1(ex, ey, ez, 6);
	dN7 = shapefun_d1(ex, ey, ez, 7);
	dN8 = shapefun_d1(ex, ey, ez, 8);
		
	for (int i = 0; i < grid.Ne; ++i){									// For each element within the mesh...
		elem = grid.ebin[i];												// Obtain the element
		
		// Construct the Jacobian matrix
		Jacob_mat.M11 = ((*elem.n1).x)*dN1.x + ((*elem.n2).x)*dN2.x + ((*elem.n3).x)*dN3.x + ((*elem.n4).x)*dN4.x + ((*elem.n5).x)*dN5.x + ((*elem.n6).x)*dN6.x + ((*elem.n7).x)*dN7.x + ((*elem.n8).x)*dN8.x;
		Jacob_mat.M21 = ((*elem.n1).x)*dN1.y + ((*elem.n2).x)*dN2.y + ((*elem.n3).x)*dN3.y + ((*elem.n4).x)*dN4.y + ((*elem.n5).x)*dN5.y + ((*elem.n6).x)*dN6.y + ((*elem.n7).x)*dN7.y + ((*elem.n8).x)*dN8.y;
		Jacob_mat.M31 = ((*elem.n1).x)*dN1.z + ((*elem.n2).x)*dN2.z + ((*elem.n3).x)*dN3.z + ((*elem.n4).x)*dN4.z + ((*elem.n5).x)*dN5.z + ((*elem.n6).x)*dN6.z + ((*elem.n7).x)*dN7.z + ((*elem.n8).x)*dN8.z;
		Jacob_mat.M12 = ((*elem.n1).y)*dN1.x + ((*elem.n2).y)*dN2.x + ((*elem.n3).y)*dN3.x + ((*elem.n4).y)*dN4.x + ((*elem.n5).y)*dN5.x + ((*elem.n6).y)*dN6.x + ((*elem.n7).y)*dN7.x + ((*elem.n8).y)*dN8.x;
		Jacob_mat.M22 = ((*elem.n1).y)*dN1.y + ((*elem.n2).y)*dN2.y + ((*elem.n3).y)*dN3.y + ((*elem.n4).y)*dN4.y + ((*elem.n5).y)*dN5.y + ((*elem.n6).y)*dN6.y + ((*elem.n7).y)*dN7.y + ((*elem.n8).y)*dN8.y;
		Jacob_mat.M32 = ((*elem.n1).y)*dN1.z + ((*elem.n2).y)*dN2.z + ((*elem.n3).y)*dN3.z + ((*elem.n4).y)*dN4.z + ((*elem.n5).y)*dN5.z + ((*elem.n6).y)*dN6.z + ((*elem.n7).y)*dN7.z + ((*elem.n8).y)*dN8.z;
		Jacob_mat.M13 = ((*elem.n1).z)*dN1.x + ((*elem.n2).z)*dN2.x + ((*elem.n3).z)*dN3.x + ((*elem.n4).z)*dN4.x + ((*elem.n5).z)*dN5.x + ((*elem.n6).z)*dN6.x + ((*elem.n7).z)*dN7.x + ((*elem.n8).z)*dN8.x;
		Jacob_mat.M23 = ((*elem.n1).z)*dN1.y + ((*elem.n2).z)*dN2.y + ((*elem.n3).z)*dN3.y + ((*elem.n4).z)*dN4.y + ((*elem.n5).z)*dN5.y + ((*elem.n6).z)*dN6.y + ((*elem.n7).z)*dN7.y + ((*elem.n8).z)*dN8.y;
		Jacob_mat.M33 = ((*elem.n1).z)*dN1.z + ((*elem.n2).z)*dN2.z + ((*elem.n3).z)*dN3.z + ((*elem.n4).z)*dN4.z + ((*elem.n5).z)*dN5.z + ((*elem.n6).z)*dN6.z + ((*elem.n7).z)*dN7.z + ((*elem.n8).z)*dN8.z;

		if (Jacob_mat.det() != 0.)											// Calculate the volume from the Jacobian matrix using the determinant
			volume = 8*Jacob_mat.det();

		grid.ebin[i].volume = volume;										// Set the element volume
		grid.ebin[i].volume0 = volume;}										// Set the element initial volume
	
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - update_grid_volume
//		Update grid volume for each element after a deformation
///////////////////////////////////////////////////////////////////////

void FEAngio::update_grid_volume()
{
	Elem elem;															// Element placeholder
	mat3 F;																// Deformation gradient tensor
	double Jacob = 0.;													// Jacobian (i.e., determinant of F)
	double ex = 0.; double ey = 0.; double ez = 0.;						// Position of the centroid in natural coordinates
	double new_volume = 0.;												// New element volume

	for (int i = 0; i < grid.Ne; ++i){									// For each element within the mesh...
		elem = grid.ebin[i];												// Obtain the element
		F = calculate_deform_tensor(elem, ex, ey, ez);						// Calculate the deformation gradient tensor
		Jacob = F.det();													// Calculate the Jacobian by taking the determinant of F

		new_volume = Jacob*elem.volume0;									// Calculate the new element volume using the Jacobian
		grid.ebin[i].volume = new_volume;}									// Store the new element volume
	
	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - assemble_sprout_nodes
//		Assemble all vessel segment nodes into a container, being careful not to repeat and track how many segments a node shares
///////////////////////////////////////////////////////////////////////

void FEAngio::assemble_sprout_nodes()
{
    list<Segment>::iterator it;											// Iterator for the segment container frag
    double xpt = 0.; double ypt = 0.; double zpt = 0.;					// Position of the sprout node
    //int pos = 1;
	double count = 1.0;													// Set the total segment node counter to 1
	bool exist = false;													// Set the exist flag to false

    for (it = frag.begin(); it != frag.end(); ++it)						// Iterate through all segments in frag list container (it)                               
	{
		for (int k=0; k<=1; ++k)											// Iterate through both segment tips (k)
		{
			exist = false;													// Set exist flag to false
			xpt = it->x[k];													// Set the position of the segment tip
			ypt = it->y[k];
			zpt = it->z[k];

			vector<double> snode;											// Create the container for the segment nodes
			snode.resize(5);												// Resize the container

			if (count == 1.0){												// If this is the first time the segment is being counted...
				snode[0] = count;												// Set the node number
				snode[1] = xpt;													// Set the x-position
				snode[2] = ypt;													// Set the y-position
				snode[3] = zpt;													// Set the z-position
				snode[4] = 1.0;													// Set the number of times this segment has been counted
				sprout_nodes.push_back(snode);									// Add the segment node to the container
				count = count + 1.0;}											// Add to the total segment node count
			else{															// If this isn't the first time the segment has been counted...
				for (int i = 0; i < int(count - 1.0); ++i)						// For each time the segment has been counted...
				{
					if ((xpt == sprout_nodes[i][1]) && (ypt == sprout_nodes[i][2]) && (zpt == sprout_nodes[i][3])){	// If we have the same segment...
						sprout_nodes[i][4] = sprout_nodes[i][4] + 1.0;				// Increase the number of times this segment has been counted by 1 
						exist = true;}												// Set exist flag to true
				}

				if (exist == false){										// If the exist flag is false
					snode[0] = count;											// Set the node number
					snode[1] = xpt;												// Set the x-position
					snode[2] = ypt;												// Set the y-position
					snode[3] = zpt;												// Set the z-position
					snode[4] = 1.0;												// Set the number of times this segment has been counted
					sprout_nodes.push_back(snode);								// Add the segment node to the container
					count = count + 1.0;}										// Add to the total segment node count
			}
		}
	}

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - output_params
//		Output parameter values after the simulation ends
///////////////////////////////////////////////////////////////////////

void FEAngio::output_params()
{
	FILE *param_stream;													// Parameter output file stream                                                                                                                     
	param_stream = fopen("out_params.ang","wt");						// Output the parameter output file
	
	fprintf(param_stream,"a = %5.5f \n",sproutf);						// Print the sprout force magnitude
	fprintf(param_stream,"tip range = %5.5f \n",tip_range);				// Print the sprout force range
	fprintf(param_stream,"phi_stiff_factor = %5.5f \n",phi_stiff_factor);	// Print the displacement stiffness factor
	fprintf(param_stream,"total_body_force = %10.5i \n",total_bdyf);		// Print the total number of body forces
	fclose(param_stream);												// Close the parameter output file

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - save_cropped_vessels
//		Output parameter values after the simulation ends
///////////////////////////////////////////////////////////////////////

void FEAngio::save_cropped_vessels()
{
	FILE *data_stream;                                                          // Open stream to 'grid.ang' (stream2)                                                                   
	data_stream = fopen("out_cropped_data.ang","wt");
	
	FILE *seg_conn_stream;                                                           // Open stream to 'out_seg_conn.ang' (stream)
	seg_conn_stream = fopen("out_cropped_seg_conn.ang","wt");

	list<Segment>::iterator it;
		
	double xmin = 0.; double xmax = 0.; double ymin = 0.; double ymax = 0.; double zmin = 0.; double zmax = 0.;
	
	xmin = ((grid.xrange[1] + grid.xrange[0])/2) - 2548/2;
	xmax = ((grid.xrange[1] + grid.xrange[0])/2) + 2548/2;
	ymin = ((grid.yrange[1] + grid.yrange[0])/2) - 2548/2;
	ymax = ((grid.yrange[1] + grid.yrange[0])/2) + 2548/2;
	zmin = ((grid.zrange[1] + grid.zrange[0])/2) - 1500/2;
	zmax = ((grid.zrange[1] + grid.zrange[0])/2) + 1500/2;

	for (it = frag.begin(); it != frag.end(); ++it)                         // Iterate through all segments in frag list container (it)
	{
		if (((it->x[0] <= xmax) && (it->x[0] >= xmin)) && ((it->x[1] <= xmax) && (it->x[1] >= xmin))){
			if (((it->y[0] <= ymax) && (it->y[0] >= ymin)) && ((it->y[1] <= ymax) && (it->y[1] >= ymin))){
				if (((it->z[0] <= zmax) && (it->z[0] >= zmin)) && ((it->z[1] <= zmax) && (it->z[1] >= zmin))){
					fprintf(data_stream,"%-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-12.7f %-5.2i %-5.2i\n",it->TofBirth,it->x[0],it->y[0],it->z[0],it->x[1],it->y[1],it->z[1],it->length,it->seg_num,it->label);  // Write to data.ang		
					fprintf(seg_conn_stream,"%-5.2i %-5.2i %-5.2i %-5.2i %-5.2i\n",it->seg_num,it->seg_conn[0][0],it->seg_conn[0][1],it->seg_conn[1][0],it->seg_conn[1][1]);  // Write to seg_conn.ang
				}
			}
		}
	}
	
	fclose(data_stream);
	fclose(seg_conn_stream);

	return;

}


///////////////////////////////////////////////////////////////////////
// FEAngio - setup_sprout_verification
//		Set up the sprout verification simulation
///////////////////////////////////////////////////////////////////////

void FEAngio::setup_sprout_verification()
{
	list<Segment>::iterator it;
	it = frag.begin();

	double xpos = 0;
	double ypos = 50;
	double zpos = 0;

	it->x[1] = xpos; it->x0[1] = it->x[1];
	it->y[1] = ypos; it->y0[1] = it->y[1];
	it->z[1] = zpos; it->z0[1] = it->z[1];

	it->uvect.x = 0;
	it->uvect.y = 1;
	it->uvect.z = 0;
	
	it->tip[1] = 1;
	it->tip[0] = 0;

	it->x[0] = xpos; it->x0[0] = it->x[0];
	it->y[0] = grid.yrange[0];
	it->z[0] = xpos; it->z0[0] = it->z[0];

	it->length = it->y[1] - it->y[0];
	it->mark_of_death = false;

	return;
}


///////////////////////////////////////////////////////////////////////
// FEAngio - update_ecm_den_grad
//		Calculate the density gradient for each element (this function may not work properly)
///////////////////////////////////////////////////////////////////////

void FEAngio::update_ecm_den_grad()
{
	Elem elem;
	vect3 elem_den_grad;
	double ex = 0.; double ey = 0.; double ez = 0.;	
	vect3 dN1; vect3 dN2; vect3 dN3; vect3 dN4; vect3 dN5; vect3 dN6; vect3 dN7; vect3 dN8;
	
	dN1 = shapefun_d1(ex, ey, ez, 1);
    dN2 = shapefun_d1(ex, ey, ez, 2);
	dN3 = shapefun_d1(ex, ey, ez, 3);
	dN4 = shapefun_d1(ex, ey, ez, 4);
	dN5 = shapefun_d1(ex, ey, ez, 5);
	dN6 = shapefun_d1(ex, ey, ez, 6);
	dN7 = shapefun_d1(ex, ey, ez, 7);
	dN8 = shapefun_d1(ex, ey, ez, 8);
	
	for (int j = 0; j < grid.Nn; ++j){
		grid.nodes[j].updated = false;}


	for (int i = 0; i < grid.Ne; ++i){
		elem = grid.ebin[i];
		
		elem_den_grad.x = ((*elem.n1).ecm_den)*dN1.x + ((*elem.n2).ecm_den)*dN2.x + ((*elem.n3).ecm_den)*dN3.x + ((*elem.n4).ecm_den)*dN4.x + ((*elem.n5).ecm_den)*dN5.x + ((*elem.n6).ecm_den)*dN6.x + ((*elem.n7).ecm_den)*dN7.x + ((*elem.n8).ecm_den)*dN8.x;
		elem_den_grad.y = ((*elem.n1).ecm_den)*dN1.y + ((*elem.n2).ecm_den)*dN2.y + ((*elem.n3).ecm_den)*dN3.y + ((*elem.n4).ecm_den)*dN4.y + ((*elem.n5).ecm_den)*dN5.y + ((*elem.n6).ecm_den)*dN6.y + ((*elem.n7).ecm_den)*dN7.y + ((*elem.n8).ecm_den)*dN8.y;
		elem_den_grad.z = ((*elem.n1).ecm_den)*dN1.z + ((*elem.n2).ecm_den)*dN2.z + ((*elem.n3).ecm_den)*dN3.z + ((*elem.n4).ecm_den)*dN4.z + ((*elem.n5).ecm_den)*dN5.z + ((*elem.n6).ecm_den)*dN6.z + ((*elem.n7).ecm_den)*dN7.z + ((*elem.n8).ecm_den)*dN8.z;
	
		if ((*elem.n1).updated == false){
			(*elem.n1).ecm_den_grad = elem_den_grad;
			(*elem.n1).updated = true;}
		else if ((*elem.n1).updated == true){
			(*elem.n1).ecm_den_grad = (elem_den_grad + (*elem.n1).ecm_den_grad)*0.5;}

		if ((*elem.n2).updated == false){
			(*elem.n2).ecm_den_grad = elem_den_grad;
			(*elem.n2).updated = true;}
		else if ((*elem.n2).updated == true){
			(*elem.n2).ecm_den_grad = (elem_den_grad + (*elem.n2).ecm_den_grad)*0.5;}
		
		if ((*elem.n3).updated == false){
			(*elem.n3).ecm_den_grad = elem_den_grad;
			(*elem.n3).updated = true;}
		else if ((*elem.n3).updated == true){
			(*elem.n3).ecm_den_grad = (elem_den_grad + (*elem.n3).ecm_den_grad)*0.5;}

		if ((*elem.n4).updated == false){
			(*elem.n4).ecm_den_grad = elem_den_grad;
			(*elem.n4).updated = true;}
		else if ((*elem.n4).updated == true){
			(*elem.n4).ecm_den_grad = (elem_den_grad + (*elem.n4).ecm_den_grad)*0.5;}

		if ((*elem.n5).updated == false){
			(*elem.n5).ecm_den_grad = elem_den_grad;
			(*elem.n5).updated = true;}
		else if ((*elem.n5).updated == true){
			(*elem.n5).ecm_den_grad = (elem_den_grad + (*elem.n5).ecm_den_grad)*0.5;}

		if ((*elem.n6).updated == false){
			(*elem.n6).ecm_den_grad = elem_den_grad;
			(*elem.n6).updated = true;}
		else if ((*elem.n6).updated == true){
			(*elem.n6).ecm_den_grad = (elem_den_grad + (*elem.n6).ecm_den_grad)*0.5;}

		if ((*elem.n7).updated == false){
			(*elem.n7).ecm_den_grad = elem_den_grad;
			(*elem.n7).updated = true;}
		else if ((*elem.n7).updated == true){
			(*elem.n7).ecm_den_grad = (elem_den_grad + (*elem.n7).ecm_den_grad)*0.5;}

		if ((*elem.n8).updated == false){
			(*elem.n8).ecm_den_grad = elem_den_grad;
			(*elem.n8).updated = true;}
		else if ((*elem.n8).updated == true){
			(*elem.n8).ecm_den_grad = (elem_den_grad + (*elem.n8).ecm_den_grad)*0.5;}
	}
	
	return;
}


void FEAngio::update_sprout_stress_scaling()
{
	double y0 = -0.004; double x0 = 3.0; double b = 0.5436; double a = 1.0081;

	if (m_pmat)
		m_pmat->scale = y0 + a/(1 + exp(-(data.t - x0)/b));
	
	return;
}

void FEAngio::circ_gel()
{
	Elem elem;
	double xmax = grid.xrange[1];
	double ymax = grid.yrange[1];

	for (int i = 0; i < grid.Ne; ++i){
		elem = grid.ebin[i];

		// Check face 1 for right boundary
		if ((((*elem.n1).x == xmax) && ((*elem.n2).x == xmax) && ((*elem.n5).x == xmax) && ((*elem.n6).x == xmax)) || (((*elem.n1).y == ymax) && ((*elem.n2).y == ymax) && ((*elem.n5).y == ymax) && ((*elem.n6).y == ymax))){
			grid.ebin[i].f1.BC = true;
			grid.ebin[i].f1.bc_type = 'w';
		}

		// Check face 2 for right boundary
		if ((((*elem.n2).x == xmax) && ((*elem.n4).x == xmax) && ((*elem.n6).x == xmax) && ((*elem.n8).x == xmax)) || (((*elem.n2).y == ymax) && ((*elem.n4).y == ymax) && ((*elem.n6).y == ymax) && ((*elem.n8).y == ymax))){
			grid.ebin[i].f2.BC = true;
			grid.ebin[i].f2.bc_type = 'w';
		}

		// Check face 3 for right boundary
		if ((((*elem.n3).x == xmax) && ((*elem.n4).x == xmax) && ((*elem.n7).x == xmax) && ((*elem.n8).x == xmax)) || (((*elem.n3).y == ymax) && ((*elem.n4).y == ymax) && ((*elem.n7).y == ymax) && ((*elem.n8).y == ymax))){
			grid.ebin[i].f3.BC = true;
			grid.ebin[i].f3.bc_type = 'w';
		}

		// Check face 4 for right boundary
		if ((((*elem.n1).x == xmax) && ((*elem.n3).x == xmax) && ((*elem.n5).x == xmax) && ((*elem.n7).x == xmax)) || (((*elem.n1).y == ymax) && ((*elem.n3).y == ymax) && ((*elem.n5).y == ymax) && ((*elem.n7).y == ymax))){
			grid.ebin[i].f4.BC = true;
			grid.ebin[i].f4.bc_type = 'w';
		}

		double xpt; double ypt; 	
		vect3 r;
		double d;
		double ep = 200.;

		//// Check face 1 for circumfrential boundary
		//xpt = ((*elem.n1).x + (*elem.n2).x + (*elem.n5).x + (*elem.n6).x)/4;
		//ypt = ((*elem.n1).y + (*elem.n2).y + (*elem.n5).y + (*elem.n6).y)/4;

		//r.x = xpt - xmax; r.y = ypt - ymax;
		//d = r.norm();

		//if (abs(d - xmax) <= ep){
		//	grid.ebin[i].f1.BC = true;
		//	grid.ebin[i].f1.bc_type = 'i';}
	
		// Check face 2 for circumfrential boundary
		xpt = ((*elem.n2).x + (*elem.n4).x + (*elem.n6).x + (*elem.n8).x)/4;
		ypt = ((*elem.n2).y + (*elem.n4).y + (*elem.n6).y + (*elem.n8).y)/4;

		r.x = xpt - xmax; r.y = ypt - ymax;
		d = r.norm();

		if (fabs(d - xmax) <= ep){
			grid.ebin[i].f2.BC = true;
			grid.ebin[i].f2.bc_type = 'i';}

		//// Check face 3 for circumfrential boundary
		//xpt = ((*elem.n3).x + (*elem.n4).x + (*elem.n7).x + (*elem.n8).x)/4;
		//ypt = ((*elem.n3).y + (*elem.n4).y + (*elem.n7).y + (*elem.n8).y)/4;

		//r.x = xpt - xmax; r.y = ypt - ymax;
		//d = r.norm();

		//if (abs(d - xmax) <= ep){
		//	grid.ebin[i].f3.BC = true;
		//	grid.ebin[i].f3.bc_type = 'i';}

		//// Check face 4 for circumfrential boundary
		//xpt = ((*elem.n1).x + (*elem.n3).x + (*elem.n5).x + (*elem.n7).x)/4;
		//ypt = ((*elem.n1).y + (*elem.n3).y + (*elem.n5).y + (*elem.n7).y)/4;

		//r.x = xpt - xmax; r.y = ypt - ymax;
		//d = r.norm();

		//if (abs(d - xmax) <= ep){
		//	grid.ebin[i].f4.BC = true;
		//	grid.ebin[i].f4.bc_type = 'i';}
	}

	return;
}
