///////////////////////////////////////////////////////////////////////
// Angio.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// The ANGIO class handles all the microvessel growth within the model.
///////////////////////////////////////////////////////////////////////



#pragma once

#include <stdlib.h>
#include <list>
#include "time.h"

#include "Grid.h"
#include "Data.h"
#include "Input.h"
#include "Culture.h"
#include "vect3.h"

class Segment;

using namespace std;

class Angio
{
public:
	Angio(Input &input);
	virtual ~Angio();
	
	void seedFrag();
	void initBranch();
	void updateTime();
	void updateTotalLength();
	virtual void Growth();
	void updateLength();
	virtual void Branch(list<Segment>::iterator it);
	virtual void Fuse();
    void check4anast(list<Segment>::iterator it, int k);
	//void search_nearby_neighbors(list<Segment>::iterator it, int k, int elem_num);
	//void scan_4_segs(list<Segment>::iterator it, int k, int elem_num);
    void anastomose(double dist0, double dist1, int k, list<Segment>::iterator it, list<Segment>::iterator it2);
	list<Segment> returnFragList();
	void removeErrors();
	virtual void displace_vessels();
	void find_active_tips();
	void find_newborns();
	void update_elem_tags();
	void kill_dead_segs();
	void initialize_seg_positions();
	void enforce_global_BC();
	
public:
    bool first_adj;
    double half_cell_l;           // Half the length of one grid element in the x direction, use in determining variable time step
	double p;
	double q;
    
	list<Segment> frag;
    
	//list<vect3> seg_nodes;

	list<list<Segment>::iterator > active_tips;
	
	list<list<Segment>::iterator > newborns;

    Grid grid;
	Data data;
	Culture cult;  

	time_t ranseed;

	FILE *killed_segs_stream; 

	bool kill_off;
};


