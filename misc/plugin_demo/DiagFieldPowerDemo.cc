/*
 * Demo plugin for field diagnostics.
 *
 * Computes power for every slice and integration step.
 *
 * C. Lechner, EuXFEL
 */

#include <iostream>
#include <math.h>

#include "DiagFieldPowerDemo.h"

using namespace std;

extern "C" DiagFieldHookedBase *factory(void)
{
	DiagFieldHookedBase *p = new DiagFieldHookedDemo;

	return(p);
}

extern "C" void destroy(DiagFieldHookedBase *p)
{
	delete p;
}

DiagFieldHookedDemo::DiagFieldHookedDemo() {}
DiagFieldHookedDemo::~DiagFieldHookedDemo() {}

// After calling this function, a copy of returned obj_names should be generated by caller (it might change/become invalid at a later time)
bool DiagFieldHookedDemo::get_infos(DiagFieldHookInfos *pi)
{
	my_obj_names.clear();
	my_obj_names.push_back("my_power");
	my_obj_names.push_back("abc");
	pi->obj_names = &my_obj_names;
	pi->info_txt = "Demo for plugin diagnostics";

	return(true);
}

void DiagFieldHookedDemo::doit(DiagFieldHookData *pd)
{
	/* ALWAYS AS FIRST STEP: check version of data structure (and if the data layout is what we expect). This crashes the program if there is a mismatch. */
	verify_datastructure(pd);

	if((pd->verbose) && (pd->mpi_rank==0)) {
		cout << "DiagFieldHookedDemo: "
		     << "iz=" << pd->iz
		     << ", is=" << pd->is
		     << ", harm=" << pd->harm
		     << ", ngrid=" << pd->ngrid << endl;
	}

	int iy, ix;
	complex<double> loc = 0;
	double wei = 0;
	double sum = 0;
	for(iy=0; iy<pd->ngrid; iy++) {
		for(ix=0; ix<pd->ngrid; ix++) {
			int idx = iy*pd->ngrid+ix; // see for instance src/Core/Diagnostic.cpp, function 'DiagField::getTags' (commit id efcc090)

			loc = pd->datain->at(idx);
			wei = loc.real()*loc.real() + loc.imag()*loc.imag();
			sum += wei;
		}
	}

	double power = scale_to_power(pd, sum);

	/* data organization as in the obj_names provided by calling 'get_infos' */
	pd->dataout->at(0) = power; // this goes to 'my_power'
	pd->dataout->at(1) = 42.;   // this goes to 'abc'
}
