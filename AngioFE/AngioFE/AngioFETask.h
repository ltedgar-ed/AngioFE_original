///////////////////////////////////////////////////////////////////////
// AngioFETask.h
///////////////////////////////////////////////////////////////////////

#pragma once

#include "FECore/FECoreKernel.h"
#include "FECore/FECoreTask.h"

class AngioFETask : public FECoreTask
{
public:
	AngioFETask(FEModel* pfem);
	~AngioFETask(void);

	bool Init(const char* szfile);

	bool Run();
};
