///////////////////////////////////////////////////////////////////////
// Fileout.h
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
// The FILEOUT class writes the output files that contain the results
// of the simulation.
///////////////////////////////////////////////////////////////////////

#pragma once

#include <list>
#include <time.h>
#include <fstream>
#include <vector>

class Grid;
class Data;
class Segment;
class Angio;
class FEAngio;

using namespace std;

class Fileout
{
public:
	Fileout();
	virtual ~Fileout();
	
	void timestart();
	void writeTracking(Data &data);
	void closeTracking();
	void printStatus(Data &data);
	void dataout(FEAngio &feangio);
    void writeData(FEAngio &feangio);
	void writeGrid(Data &data, Grid &grid);
	void writeNodes(Data &data, Grid &grid);
	void writeEconn(Data &data, Grid &grid);
	void writeGrad(Data &data, Grid &grid);
	void writeAngle(list<Segment> &frag);
	void writeCollFib(Grid &grid, bool initial);
	void writeECMDen(Grid &grid);
	void writeECMDenGrad(Grid &grid);
	void writeBC(Grid &grid);
	void printtime();
	void printrandseed(int randseed);
	void printsproutnodes(vector<vector<double> > sprout_nodes);
	void writeSegConn(list<Segment> &frag);
	void writeECMDenStore(Grid &grid);
	void writeECMFibrilStore(Grid &grid);

public:
    time_t start, stop;                                         // Declare time-keeping variables start, stop, and t_seconds
	double t_seconds;
	ofstream logstream;
    FILE* stream3;      
	FILE *stream;                                                           // Open stream to 'data.ang' (stream)

	int num_fe_timesteps;
};