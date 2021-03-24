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

class Grid;
class Data;
class Segment;
class Angio;

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
	void dataout(Angio &angio, Data &data, Grid &grid);
    void writeData(list<Segment> &frag);
	void writeGrid(Data &data, Grid &grid);
	void writeNodes(Data &data, Grid &grid);
	void writeEconn(Data &data, Grid &grid);
	void writeGrad(Data &data, Grid &grid);
	void writeAngle(list<Segment> &frag);
	void writeBC(Grid &grid);
	void printtime();
	void printrandseed(int randseed);
	void writeSegConn(list<Segment> &frag);

public:
    time_t start, stop;                                         // Declare time-keeping variables start, stop, and t_seconds
	double t_seconds;
	ofstream logstream;
    FILE* stream3;                                              // Open stream to 'tracking.ang' (stream3)
};