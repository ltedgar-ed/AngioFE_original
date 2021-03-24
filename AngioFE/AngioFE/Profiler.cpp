///////////////////////////////////////////////////////////////////////
// Profiler.cpp
///////////////////////////////////////////////////////////////////////

#include "StdAfx.h"
#include "Profiler.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Profiler::Profiler()
{
    profstream.open("out_profiler.ang");
	total_time_1 = 0;
	total_time_2 = 0;
	total_time_3 = 0;
	current_profiler = 0;
}

Profiler::~Profiler()
{
    profstream << "Time in fem.Solve: " << total_time_1 << endl;
	profstream << "Time in FESproutBodyForce.force: " << total_time_2 << endl;
	profstream << "Time in FESproutBodyForce.stiff: " << total_time_3 << endl;
	profstream.close();
}



///////////////////////////////////////////////////////////////////////
// Member Functions
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
// time_in
///////////////////////////////////////////////////////////////////////

void Profiler::time_in(int profiler)
{
	current_profiler = profiler;
	
	time(&start);
	
	return;
}



///////////////////////////////////////////////////////////////////////
// time_out
///////////////////////////////////////////////////////////////////////

void Profiler::time_out()
{
	time(&stop);					                                 // Stop the timer
	
	if (current_profiler == 1)
		total_time_1 = total_time_1 + difftime(stop, start);                 // Calculate the simulation time in seconds
	
	if (current_profiler == 2)
		total_time_2 = total_time_2 + difftime(stop, start); 

	if (current_profiler == 3)
		total_time_3 = total_time_3 + difftime(stop, start); 
	
	return;
}
