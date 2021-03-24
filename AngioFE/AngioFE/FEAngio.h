///////////////////////////////////////////////////////////////////////
// FEAngio.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The ANGIO class contains all the functionality of the growth model.
///////////////////////////////////////////////////////////////////////

#pragma once

#include "Angio.h"
#include "FECore/FEModel.h"
#include "Fileout.h"
#include "FEAngInput.h"
#include "Profiler.h"
#include "FEAngioMaterial.h"
#include "FESproutBodyForce.h"

class FEAngio : public Angio
{
// Public functions:
public:
	FEAngio(FEAngInput &input);
	~FEAngio();

	void Growth(FEModel& fem);
	void Branch(list<Segment>::iterator it, FEModel& fem);
	void Subgrowth(int sub_steps, FEModel& fem);
	void displace_vessels();
	void apply_sprout_forces(FEModel& fem, int load_curve, double scale);
	int create_body_force(vec3d sprout_vect, double xpt, double ypt, double zpt, double mag, double range, FEModel& fem, int load_curve);
	void update_body_forces(FEModel& fem, double scale);
	void create_branching_force(Segment& seg, FEModel& fem);
	void save_vessel_state();
	void save_bdy_forces(FEModel& fem);
	void save_time();
	void update_grid(FEMesh &mesh);
	void update_ECM();
	mat3 calculate_deform_tensor(Elem elem, double e1, double e2, double e3);
	vect3 shapefun_d1(const double xix, const double xiy, const double xiz, int node);
	void adjust_mesh_stiffness(FEModel& fem);
	void initialize_grid_volume();
	void update_grid_volume();
	void update_ecm_den_grad();
	void assemble_sprout_nodes();
	void output_params();
	void save_cropped_vessels();
	void setup_sprout_verification();
	void update_sprout_stress_scaling();
	void circ_gel();

//--> SAM
	void update_angio_sprout(int i, bool bactive, const vec3d& rc, const vec3d& sprout_vect);
//<-- SAM

// Public fields:
public:
	double sproutf;
	double tip_range;
	double b;
	
	double phi_stiff_factor;

	int total_bdyf;
	double FE_time_step;
	int FE_state;
	
	bool yes_branching;
	bool yes_anast;

	Fileout fileout;
	
	FILE *stream;                                                           
	FILE *bf_stream;

	FILE *time_stream;
	bool time_write_headers;

	int comp_mat;

	vector<vector<double> > sprout_nodes;

	bool feconv;													// Flag indicating convergence of the FE model

	int sub_cycles;

	bool sprout_verify;

	bool disp_vess;

//--> SAM
public:
	FESproutBodyForce*	m_pbf;	//!< sprout body-force
	FEAngioMaterial*	m_pmat;	//!< the angio-material pointer
//<-- SAM

	//Profiler profiler;
};

