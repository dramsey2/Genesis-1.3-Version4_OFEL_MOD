#include "TrackBeam.h"
#include "Field.h"
#include "Beam.h"


TrackBeam::TrackBeam(){}
TrackBeam::~TrackBeam(){}


void TrackBeam::track(double delz, Beam *beam,Undulator *und, EFieldSolver efield, bool lastStep=true)
{

  // get undulator parameter for the given step
  double aw,dax,day,ku,kx,ky;
  double qf,dqx,dqy;
  double cx,cy;
  double angle,lb,ld,lt;

  double gamma0=und->getGammaRef();
  und->getUndulatorParameters( &aw,&dax,&day,&ku,&kx,&ky);
  und->getQuadrupoleParameters(&qf,&dqx,&dqy);
  und->getCorrectorParameters(&cx,&cy);
  und->getChicaneParameters(&angle,&lb,&ld,&lt);


  double betpar0=sqrt(1-(1+aw*aw)/gamma0/gamma0);
  
  // effective focusing in x and y with the correct energy dependence
  double qquad=qf*gamma0;
  double qnatx=kx*aw*aw/gamma0/betpar0;  // kx has already the scaling with ku^2
  double qnaty=ky*aw*aw/gamma0/betpar0;  // same with ky

  double qx= qquad+qnatx;
  double qy=-qquad+qnaty;

  double xoff= qquad*dqx+qnatx*dax;
  double yoff=-qquad*dqy+qnaty*day;
  // cout << "qnaty: " << qnaty/gamma0 << " Gamma0: " << gamma0 << endl;

  if (lastStep){
    if ((cx!=0) || (cy!=0)) { this->applyCorrector(beam,cx*gamma0,cy*gamma0); }
  } else {
    if (angle!=0) { this->applyChicane(beam,angle,lb,ld,lt,gamma0); }
  }
  // handle the different cases (drift, focusing and defocusing) with function pointers to member functions

  if (qx==0){ 
    this->ApplyX=&TrackBeam::applyDrift; 
  }else{
    xoff=xoff/qx;
    if (qx>0){
      this->ApplyX=&TrackBeam::applyFQuad;
    } else {
      this->ApplyX=&TrackBeam::applyDQuad;
    }
  }

  if (qy==0){ 
    this->ApplyY=&TrackBeam::applyDrift; 
  }else{
    yoff=yoff/qy;
    if (qy>0){
      this->ApplyY=&TrackBeam::applyFQuad;
    } else {
      this->ApplyY=&TrackBeam::applyDQuad;
    }
  }

  for (int i = 0; i < beam->beam.size(); i++)
  {
    efield.shortRange(&beam->beam.at(i), beam->current.at(i), sqrt(1/(1-betpar0*betpar0)), i);

    for (int j = 0; j < beam->beam.at(i).size(); j++)
    {
      Particle *p = &beam->beam.at(i).at(j);
      er = efield.getERField(j);

      double gammaz = sqrt(p->gamma * p->gamma - 1 - aw * aw - p->px * p->px - p->py * p->py); // = gamma*betaz=gamma*(1-(1+aw*aw)/gamma^2);
#ifdef G4_DBGDIAG
// G4_DBGDIAG: add test against negative radicand? Note that the particles probably already made lots of noise elsewhere.
#endif
      (this->*ApplyX)(delz,qx, kx, &(p->x), &(p->y), &(p->px), &(p->py), gammaz, xoff);
	    // Do not apply Y. It will be updated in ApplyX through the RK4
      //(this->*ApplyY)(delz,qy, ky, &(p->y), x_temp, &(p->py),gammaz,yoff);
    }
  }

  return;
} 


void TrackBeam::applyDrift(double delz, double qf, double kx, double* x, double* y, double* px, double* py, double gammaz, double dx)
{
  *x+=(*px)*delz/gammaz;
  return;
}


void TrackBeam::applyFQuad(double delz, double qf, double kx, double* x, double* y, double* px, double* py, double gammaz, double dx)
{
  double foc=sqrt(qf/gammaz);
  double omg=foc*delz;
  double a1=cos(omg);
  double a2=sin(omg)/foc;
  double a3=-a2*foc*foc;
  double xtmp=*x-dx;
  *x =a1*xtmp+a2*(*px)/gammaz+dx;
  *px=a3*xtmp*gammaz+a1*(*px);
  return;
}


/*
void TrackBeam::applyDQuad(double delz, double qf, double kx, double* x, double* px, double gammaz, double dx)
{
  double foc=sqrt(-qf/gammaz);
  double omg=foc*delz;
  double a1=cosh(omg);
  double a2=sinh(omg)/foc;
  double a3=a2*foc*foc;
  double xtmp=*x-dx;
  *x =a1*xtmp+a2*(*px)/gammaz+dx;
  *px=a3*xtmp*gammaz+a1*(*px);
  return;
}
*/


void TrackBeam::applyDQuad(double delz, double qf, double kx, double* x, double* y, double* px, double* py, double gammaz, double dx)
{
  double betpar0 = sqrt(1-1/(gammaz*gammaz));

	double focsq = -1*qf / gammaz;
	double dxgdz = (*px);
	double dygdz = (*py);
	double xtmp = (*x);
	double ytmp = (*y);


	// first step
	double K2dxgdz = 0;
	double K2dygdz = 0;
	double K2x = 0;
	double K2y = 0;

	double kxrsq = kx * (xtmp * xtmp + ytmp * ytmp);
  double angle = atan2(ytmp,xtmp);
	K2dxgdz += gammaz*focsq * kxrsq * (kxrsq - 2) * xtmp * exp(kxrsq) / 4/betpar0 - (1-betpar0*betpar0)*er*cos(angle)/betpar0;
	K2dygdz += gammaz*focsq * kxrsq * (kxrsq - 2) * ytmp * exp(kxrsq) / 4/betpar0 - (1-betpar0*betpar0)*er*sin(angle)/betpar0;
	K2x += dxgdz/gammaz;
	K2y += dygdz/gammaz;


	// second step
	dxgdz += 0.5 * delz * K2dxgdz;
	dygdz += 0.5 * delz * K2dygdz;
	xtmp += 0.5 * delz * K2x;
	ytmp += 0.5 * delz * K2y;

	double K3dxgdz = K2dxgdz;
	double K3dygdz = K2dygdz;
	double K3x = K2x;
	double K3y = K2y;

	K2dxgdz = 0;
	K2dygdz = 0;
	K2x = 0;
	K2y = 0;

	kxrsq = kx * (xtmp * xtmp + ytmp * ytmp);
  angle = atan2(ytmp,xtmp);
	K2dxgdz += gammaz*focsq * kxrsq * (kxrsq - 2) * xtmp * exp(kxrsq) / 4/betpar0 - (1-betpar0*betpar0)*er*cos(angle)/betpar0;
	K2dygdz += gammaz*focsq * kxrsq * (kxrsq - 2) * ytmp * exp(kxrsq) / 4/betpar0 - (1-betpar0*betpar0)*er*sin(angle)/betpar0;
	K2x += dxgdz/gammaz;
	K2y += dygdz/gammaz;

	// third step
	dxgdz += 0.5 * delz * (K2dxgdz - K3dxgdz);
	dygdz += 0.5 * delz * (K2dygdz - K3dygdz);
	xtmp += 0.5 * delz * (K2x - K3x);
	ytmp += 0.5 * delz * (K2y - K3y);

	K3dxgdz /= 6;
	K3dygdz /= 6;
	K3x /= 6;
	K3y /= 6;

	K2dxgdz *= -0.5;
	K2dygdz *= -0.5;
	K2x *= -0.5;
	K2y *= -0.5;

	kxrsq = kx * (xtmp * xtmp + ytmp * ytmp);
  angle = atan2(ytmp,xtmp);
	K2dxgdz += gammaz*focsq * kxrsq * (kxrsq - 2) * xtmp * exp(kxrsq) / 4/betpar0 - (1-betpar0*betpar0)*er*cos(angle)/betpar0;
	K2dygdz += gammaz*focsq * kxrsq * (kxrsq - 2) * ytmp * exp(kxrsq) / 4/betpar0 - (1-betpar0*betpar0)*er*sin(angle)/betpar0;
	K2x += dxgdz/gammaz;
	K2y += dygdz/gammaz;

	// fourth step
	dxgdz += delz * K2dxgdz;
	dygdz += delz * K2dygdz;
	xtmp += delz * K2x;
	ytmp += delz * K2y;

	K3dxgdz -= K2dxgdz;
	K3dygdz -= K2dygdz;
	K3x -= K2x;
	K3y -= K2y;

	K2dxgdz *= 2;
	K2dygdz *= 2;
	K2x *= 2;
	K2y *= 2;

	kxrsq = kx * (xtmp * xtmp + ytmp * ytmp);
  angle = atan2(ytmp,xtmp);
	K2dxgdz += gammaz*focsq * kxrsq * (kxrsq - 2) * xtmp * exp(kxrsq) / 4/betpar0 - (1-betpar0*betpar0)*er*cos(angle)/betpar0;
	K2dygdz += gammaz*focsq * kxrsq * (kxrsq - 2) * ytmp * exp(kxrsq) / 4/betpar0 - (1-betpar0*betpar0)*er*sin(angle)/betpar0;
	K2x += dxgdz/gammaz;
	K2y += dygdz/gammaz;

	dxgdz += delz * (K3dxgdz + K2dxgdz / 6.0);
	dygdz += delz * (K3dygdz + K2dygdz / 6.0);
	xtmp += delz * (K3x + K2x / 6.0);
	ytmp += delz * (K3y + K2y / 6.0);


	// pass out
	*px = dxgdz;
	*py = dygdz;
	*x = xtmp;
	*y = ytmp;


  return;
}

void TrackBeam::applyCorrector(Beam* beam, double cx, double cy)
{

	for (int i = 0; i < beam->beam.size(); i++) {
		for (int j = 0; j < beam->beam.at(i).size(); j++) {
			beam->beam.at(i).at(j).px += cx;
			beam->beam.at(i).at(j).py += cy;
		}
	}
	return;
}


void TrackBeam::applyChicane(Beam *beam, double angle, double lb, double ld, double lt, double gamma0)
{ 
  // the tracking is done my applying the transfer matrix for the chicane and  backtracking for a drift over the length of the chicane
  // the effect of the R56 is applied here to the particle phase.  Then the normal tracking should do the momentum dependent change in the 
  // longitudinal position

  // the transfer matrix order is
  //  m -> bp -> ep -> d1 -> en -> bn -> d2 -> bn -> en-> d1 -> ep-> bp ->d3
  

  double m[4][4];
  double d1[4][4];
  double d2[4][4];
  double d3[4][4];
  double bp[4][4];
  double bn[4][4];
  double ep[4][4];
  double en[4][4];

  // construct the transfer matrix
  for (int i=0; i<4;i++){
    for (int j=0; j<4; j++){
      m[i][j]=0;
      d1[i][j]=0;
      d2[i][j]=0;
      d3[i][j]=0;
      bp[i][j]=0;
      bn[i][j]=0;
      ep[i][j]=0;
      en[i][j]=0;
    }
    m[i][i]=1;
    d1[i][i]=1;
    d2[i][i]=1;
    d3[i][i]=1;
    bp[i][i]=1;
    bn[i][i]=1;
    ep[i][i]=1;
    en[i][i]=1;
  }
  d1[0][1]=ld/cos(angle);  // drift between dipoles
  d1[2][3]=ld/cos(angle);   
  d2[0][1]=lt-4*lb-2*ld;   // drift in the middle
  d2[2][3]=lt-4*lb-2*ld;   
  d3[0][1]=-lt;            // negative drift over total chicane to get a zero element
  d3[2][3]=-lt;

  double R=lb/sin(angle);
  double Lpath=R*angle; 
  bp[2][3]=Lpath;  // positive deflection angle
  bp[0][0]=cos(angle);
  bp[0][1]=R*sin(angle);
  bp[1][0]=-sin(angle)/R;
  bp[1][1]=cos(angle);
 
  bn[2][3]=Lpath; // negative deflection angle
  bn[0][0]=cos(-angle);
  bn[0][1]=R*sin(-angle)*-1;
  bn[1][0]=-sin(-angle)/R*-1;
  bn[1][1]=cos(-angle);

  double efoc=tan(angle)/R;
  ep[1][0]=efoc;
  ep[3][2]=-efoc;
  en[1][0]=-efoc*-1;
  en[3][2]=efoc*-1;


  this->matmul(m,bp);
  this->matmul(m,ep);  
  this->matmul(m,d1);
  this->matmul(m,en);
  this->matmul(m,bn);
  this->matmul(m,d2);
  this->matmul(m,bn);
  this->matmul(m,en);
  this->matmul(m,d1);
  this->matmul(m,ep);
  this->matmul(m,bp);

  // transport matrix has been cross checked with Madx.  
  /*
  cout << "lt = " << lt << " angle = " << angle << " lb = " << lb << " ld = " << ld <<  endl;
  cout << m[0][0] << " " << m[0][1] << endl;
  cout << m[1][0] << " " << m[1][1] << endl;
  cout << m[2][2] << " " << m[2][3] << endl;
  cout << m[3][2] << " " << m[3][3] << endl;
  */

  this->matmul(m,d3);  // transport backwards because the main tracking still has to do the drift
  

  for (int i=0; i<beam->beam.size();i++){
    for (int j=0; j<beam->beam.at(i).size();j++){
      Particle *p=&beam->beam.at(i).at(j);
      double gammaz=sqrt(p->gamma*p->gamma-1- p->px*p->px - p->py*p->py); // = gamma*betaz=gamma*(1-(1+aw*aw)/gamma^2);

      double tmp=p->x;
      p->x =m[0][0]*tmp        +m[0][1]*p->px/gammaz;
      p->px=m[1][0]*tmp*gammaz +m[1][1]*p->px;
      tmp=p->y;
      p->y =m[2][2]*tmp        +m[2][3]*p->py/gammaz;
      p->py=m[3][2]*tmp*gammaz +m[3][3]*p->py;
      
    }
  }

  return;
}

void TrackBeam::matmul(double m[][4], double e[][4])
{
  double t[4][4];

  for (int i=0;i<4;i++){
    for (int j=0; j<4;j++){
      t[i][j]=0;
      for (int k=0;k<4;k++){
	t[i][j]+=e[i][k]*m[k][j];
      }
    }
  }
  for (int i=0;i<4;i++){
    for (int j=0; j<4;j++){
      m[i][j]=t[i][j];
    }
  }
}

void TrackBeam::applyR56(Beam *beam, Undulator *und, double lambda0)
{
  double angle,lb,ld,lt;

  double gamma0=und->getGammaRef();
  und->getChicaneParameters(&angle,&lb,&ld,&lt);
  if (angle==0) { return;}
  double R56=(4*lb/sin(angle)*(1-angle/tan(angle))+2*ld*tan(angle)/cos(angle))*angle;
  //    cout << "R56: " << R56 << endl;
  R56=R56*4*asin(1)/lambda0/gamma0;
  for (int i=0; i<beam->beam.size();i++){
    for (int j=0; j<beam->beam.at(i).size();j++){
      beam->beam.at(i).at(j).theta+=R56*(beam->beam.at(i).at(j).gamma-gamma0);
    }
  }
  return;

}
