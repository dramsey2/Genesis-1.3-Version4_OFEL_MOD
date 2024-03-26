#ifndef __GENESIS_TRACKBEAM__
#define __GENESIS_TRACKBEAM__

#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <map>
#include <stdlib.h>

#include "EFieldSolver.h"

class Beam;
#include "Undulator.h"
#include "Particle.h"


using namespace std;

class TrackBeam{
 public:
   TrackBeam();
   virtual ~TrackBeam();
   void track(double, Beam *, Undulator *,  EFieldSolver , bool);
   //void initEField(double rmax, int ngrid, int nz, int nphi, double lambda, bool longr);
   void (TrackBeam::* ApplyX)  (double, double, double, double*, double*, double*, double*, double, double);
   void (TrackBeam::* ApplyY)   (double, double, double, double*, double*, double*, double*, double, double);
   void applyDrift(double, double, double, double*, double*, double*, double*, double, double);
   void applyFQuad(double, double, double, double*, double*, double*, double*, double, double);
   void applyDQuad(double, double, double, double*, double*, double*, double*, double, double);
   void applyCorrector(Beam *, double, double);
   void applyChicane(Beam *, double, double, double, double,double);
   void applyR56(Beam *, Undulator *, double);

 private:

   double er{};
   //EFieldSolver efield;
   void matmul(double a[][4], double b[][4]);

};

//inline void TrackBeam::initEField(double rmax, int ngrid, int nz, int nphi, double lambda, bool longr) {
//    efield.init(rmax, ngrid, nz, nphi, lambda, longr);
//}

#endif
