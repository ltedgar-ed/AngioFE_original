///////////////////////////////////////////////////////////////////////
// Profiler.h
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// Description
///////////////////////////////////////////////////////////////////////



#pragma once

#include <time.h>
#include <fstream>

using namespace std;

class Profiler
{
public:
	Profiler();
	virtual ~Profiler();
	
	void time_in(int profiler);
	void time_out();

public:
	ofstream profstream;

	time_t start, stop;                                         // Declare time-keeping variables start, stop, and t_seconds
	
	double total_time_1, total_time_2, total_time_3;
	int current_profiler;
};