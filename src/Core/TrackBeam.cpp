#include "TrackBeam.h"
#include "Beam.h"

TrackBeam::TrackBeam() {}
TrackBeam::~TrackBeam() {}

void TrackBeam::track(double delz, Beam *beam, Undulator *und, EFieldSolver efield, bool lastStep = true)
{

  // get undulator parameter for the given step
  double aw, dax, day, ku, kx, ky;
  double qf, dqx, dqy;
  double cx, cy;
  double angle, lb, ld, lt;

  double gamma0 = und->getGammaRef();
  und->getUndulatorParameters(&aw, &dax, &day, &ku, &kx, &ky);
  und->getQuadrupoleParameters(&qf, &dqx, &dqy);
  und->getCorrectorParameters(&cx, &cy);
  und->getChicaneParameters(&angle, &lb, &ld, &lt);
  double z = und->getz();

  // First pass to pick which to add to apply

  betpar0 = sqrt(1 - (1 + aw * aw) / gamma0 / gamma0);

  // effective focusing in x and y with the correct energy dependence
  double qquad = qf * gamma0;
  double qnatx = kx * aw * aw / gamma0 / betpar0; // kx has already the scaling with ku^2
  double qnaty = ky * aw * aw / gamma0 / betpar0; // same with ky

  double qx = qquad + qnatx;
  double qy = -qquad + qnaty;

  double xoff = qquad * dqx + qnatx * dax;
  double yoff = -qquad * dqy + qnaty * day;
  // cout << "qnaty: " << qnaty/gamma0 << " Gamma0: " << gamma0 << endl;

  if (lastStep)
  {
    if ((cx != 0) || (cy != 0))
    {
      this->applyCorrector(beam, cx * gamma0, cy * gamma0);
    }
  }
  else
  {
    if (angle != 0)
    {
      this->applyChicane(beam, angle, lb, ld, lt, gamma0);
    }
  }
  // handle the different cases (drift, focusing and defocusing) with function pointers to member functions

  if (qx == 0)
  {
    this->ApplyX = &TrackBeam::applyDrift;
  }
  else
  {
    xoff = xoff / qx;
    if (qx > 0)
    {
      this->ApplyX = &TrackBeam::applyFQuad;
    }
    else
    {
      this->ApplyX = &TrackBeam::applyDQuad;
    }
  }

  if (qy == 0)
  {
    this->ApplyY = &TrackBeam::applyDrift;
  }
  else
  {
    yoff = yoff / qy;
    if (qy > 0)
    {
      this->ApplyY = &TrackBeam::applyFQuad;
    }
    else
    {
      this->ApplyY = &TrackBeam::applyDQuad;
    }
  }

  for (int i = 0; i < beam->beam.size(); i++)
  {

    efield.shortRange(&beam->beam.at(i), beam->current.at(i), sqrt(1 / (1 - betpar0 * betpar0)), i);

    for (int j = 0; j < beam->beam.at(i).size(); j++)
    {
      Particle *p = &beam->beam.at(i).at(j);
      er = efield.getERField(j);
      double awsq_at_p =  (und->faw(p->x, p->y))*(und->faw(p->x, p->y))*aw*aw; // should use the value of aw where the electron is
      //double gammaz = sqrt(p->gamma * p->gamma - 1 - aw * aw - p->px * p->px - p->py * p->py); // = gamma*betaz=gamma*(1-(1+aw*aw)/gamma^2);
      double gammaz = sqrt(p->gamma * p->gamma - 1 - awsq_at_p - p->px * p->px - p->py * p->py); // = gamma*betaz=gamma*(1-(1+aw*aw)/gamma^2);
      betpar0 = gammaz/(p->gamma);

      // use particle numbers // Note aw still has 1/(1+z)^(1/2) but is fixed in the push
      qx = kx * aw * aw / (p->gamma)/ betpar0; // kx has already the scaling with ku^2

#ifdef G4_DBGDIAG
// G4_DBGDIAG: add test against negative radicand? Note that the particles probably already made lots of noise elsewhere.
#endif
      //double x_temp = p->x;
      //double y_temp = p->y;
      (this->*ApplyX)(delz, qx, kx, ku, z, &(p->x), &(p->y), &(p->px), &(p->py), gammaz, xoff);
      // ApplyX will update Y for me
      // this->*ApplyY)(delz,qy, ky, &(p->y), x_temp, &(p->py), gammaz, yoff);
    }
  }

  return;
}

void TrackBeam::applyDrift(double delz, double qf, double kx, double ku, double z, double *x, double *y, double *px, double *py, double gammaz, double dx)
{
  *x += (*px) * delz / gammaz;
  return;
}

void TrackBeam::applyFQuad(double delz, double qf, double kx, double ku, double z, double *x, double *y, double *px, double *py, double gammaz, double dx)
{
  double foc = sqrt(qf / gammaz);
  double omg = foc * delz;
  double a1 = cos(omg);
  double a2 = sin(omg) / foc;
  double a3 = -a2 * foc * foc;
  double xtmp = *x - dx;
  *x = a1 * xtmp + a2 * (*px) / gammaz + dx;
  *px = a3 * xtmp * gammaz + a1 * (*px);
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

/*
void TrackBeam::applyDQuad(double delz, double qf, double kx, double* x, double y, double* px, double gammaz, double dx)
{

    // qf has the same sign as kx but it needs to be negated below

    double xtmp = *x - dx;
    double foc;
    double omg;

    if ((1 + 2 * kx * xtmp * xtmp) >= 0)
    {

        foc = sqrt(-qf / gammaz) * sqrt(1 + 2 * kx *( xtmp * xtmp)) * exp(kx * (xtmp * xtmp + y * y) / 2);
        omg = foc * delz;

        double a1x = 2 * kx * xtmp * xtmp * xtmp / (1 + 2 * kx * xtmp * xtmp);
        double a2x = cosh(omg) / (1 + 2 * kx * xtmp * xtmp);
        double a3x = sinh(omg) / foc;

        double a1p = cosh(omg);
        double a2p = foc * sinh(omg) / (1 + 2 * kx * xtmp * xtmp);

        *x = a1x + a2x * xtmp + a3x * (*px) / gammaz + dx;
        *px = a1p * (*px) + a2p * xtmp * gammaz;

    }
    else
    {
        foc = sqrt(-qf / gammaz) * sqrt(-1 - 2 * kx * xtmp * xtmp) * exp(kx * (xtmp * xtmp + y * y) / 2);
        omg = foc * delz;

        double a1x = 2 * kx * xtmp * xtmp * xtmp / (1 + 2 * kx * xtmp * xtmp);
        double a2x = cos(omg) / (1 + 2 * kx * xtmp * xtmp);
        double a3x = sin(omg) / foc;

        double a1p = cos(omg);
        double a2p = -1 * foc * sin(omg) / (1 + 2 * kx * xtmp * xtmp);

        *x = a1x + a2x * xtmp + a3x * (*px) / gammaz + dx;
        *px = a1p * (*px) + a2p * xtmp * gammaz;
    }
    return;
}

*/

void TrackBeam::applyDQuad(double delz, double qf, double kx, double ku, double z, double *x, double *y, double *px, double *py, double gammaz, double dx)
{
  // ensure q is negated
  // *The er terms need to be divided by w0* //
  // Remember gammaz is betaz*gamma
  //double betpar0 = sqrt(1 - 1 / (gammaz * gammaz));
  double w0 = sqrt(2) * sqrt(-1 / kx);
  double ZR = -ku / (2 * kx);
  double z_norm = z / ZR;
  double factor_on_aw = 1 / sqrt(1 + z_norm * z_norm);
  double qf_fixed = qf / factor_on_aw / factor_on_aw; // Now aw in focsq will be the max value of aw

  double focsq = -1 * qf_fixed / gammaz;

  double dxdz = (*px) / gammaz;
  double dydz = (*py) / gammaz;
  double xtmp = (*x);
  double ytmp = (*y);

  double xtmp_norm = xtmp / w0;
  double ytmp_norm = ytmp / w0;
  double dxdz_norm = dxdz / w0;
  double dydz_norm = dydz / w0;

  // first step
  double K2dxdz = 0;
  double K2dydz = 0;
  double K2x = 0;
  double K2y = 0;

  double angle = atan2(ytmp_norm, xtmp_norm);
  double exp_and_Arg = exp(-2 * (xtmp_norm * xtmp_norm + ytmp_norm * ytmp_norm) / (1 + z_norm * z_norm));
  double com = focsq * exp_and_Arg/ (1 + z_norm * z_norm)/ (1 + z_norm * z_norm);
  K2dxdz += xtmp_norm * com - (1 - betpar0 * betpar0) * er * cos(angle) / gammaz/w0;
  K2dydz += ytmp_norm * com - (1 - betpar0 * betpar0) * er * sin(angle) / gammaz/w0;
  K2x += dxdz_norm;
  K2y += dydz_norm;

  // second step
  dxdz_norm += 0.5 * delz * K2dxdz;
  dydz_norm += 0.5 * delz * K2dydz;
  xtmp_norm += 0.5 * delz * K2x;
  ytmp_norm += 0.5 * delz * K2y;

  double K3dxdz = K2dxdz;
  double K3dydz = K2dydz;
  double K3x = K2x;
  double K3y = K2y;

  K2dxdz = 0;
  K2dydz = 0;
  K2x = 0;
  K2y = 0;

  angle = atan2(ytmp_norm, xtmp_norm);
  exp_and_Arg = exp(-2 * (xtmp_norm * xtmp_norm + ytmp_norm * ytmp_norm)) / (1 + z_norm * z_norm);
  com = focsq * exp_and_Arg/ (1 + z_norm * z_norm)/ (1 + z_norm * z_norm);
  K2dxdz += xtmp_norm * com - (1 - betpar0 * betpar0) * er * cos(angle) / gammaz/w0;
  K2dydz += ytmp_norm * com - (1 - betpar0 * betpar0) * er * sin(angle) / gammaz/w0;
  K2x += dxdz_norm;
  K2y += dydz_norm;

  // third step
  dxdz_norm += 0.5 * delz * (K2dxdz - K3dxdz);
  dydz_norm += 0.5 * delz * (K2dydz - K3dydz);
  xtmp_norm += 0.5 * delz * (K2x - K3x);
  ytmp_norm += 0.5 * delz * (K2y - K3y);

  K3dxdz /= 6;
  K3dydz /= 6;
  K3x /= 6;
  K3y /= 6;

  K2dxdz *= -0.5;
  K2dydz *= -0.5;
  K2x *= -0.5;
  K2y *= -0.5;

  angle = atan2(ytmp_norm, xtmp_norm);
  exp_and_Arg = exp(-2 * (xtmp_norm * xtmp_norm + ytmp_norm * ytmp_norm)) / (1 + z_norm * z_norm);
  com = focsq * exp_and_Arg/ (1 + z_norm * z_norm)/ (1 + z_norm * z_norm);
  K2dxdz += xtmp_norm * com - (1 - betpar0 * betpar0) * er * cos(angle) / gammaz/w0;
  K2dydz += ytmp_norm * com - (1 - betpar0 * betpar0) * er * sin(angle) / gammaz/w0;
  K2x += dxdz_norm;
  K2y += dydz_norm;

  // fourth step
  dxdz_norm += delz * K2dxdz;
  dydz_norm += delz * K2dydz;
  xtmp_norm += delz * K2x;
  ytmp_norm += delz * K2y;

  K3dxdz -= K2dxdz;
  K3dydz -= K2dydz;
  K3x -= K2x;
  K3y -= K2y;

  K2dxdz *= 2;
  K2dydz *= 2;
  K2x *= 2;
  K2y *= 2;

  angle = atan2(ytmp_norm, xtmp_norm);
  exp_and_Arg = exp(-2 * (xtmp_norm * xtmp_norm + ytmp_norm * ytmp_norm)) / (1 + z_norm * z_norm);
  com = focsq * exp_and_Arg/ (1 + z_norm * z_norm)/ (1 + z_norm * z_norm);
  K2dxdz += xtmp_norm * com - (1 - betpar0 * betpar0) * er * cos(angle) / gammaz/w0;
  K2dydz += ytmp_norm * com - (1 - betpar0 * betpar0) * er * sin(angle) / gammaz/w0;
  K2x += dxdz_norm;
  K2y += dydz_norm;

  dxdz_norm += delz * (K3dxdz + K2dxdz / 6);
  dydz_norm += delz * (K3dydz + K2dydz / 6);
  xtmp_norm += delz * (K3x + K2x / 6);
  ytmp_norm += delz * (K3y + K2y / 6);

  // pass out
  *px = gammaz * dxdz_norm * w0;
  *py = gammaz * dydz_norm * w0;
  *x = xtmp_norm * w0;
  *y = ytmp_norm * w0;

  return;
}

void TrackBeam::applyCorrector(Beam *beam, double cx, double cy)
{

  for (int i = 0; i < beam->beam.size(); i++)
  {
    for (int j = 0; j < beam->beam.at(i).size(); j++)
    {
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
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      m[i][j] = 0;
      d1[i][j] = 0;
      d2[i][j] = 0;
      d3[i][j] = 0;
      bp[i][j] = 0;
      bn[i][j] = 0;
      ep[i][j] = 0;
      en[i][j] = 0;
    }
    m[i][i] = 1;
    d1[i][i] = 1;
    d2[i][i] = 1;
    d3[i][i] = 1;
    bp[i][i] = 1;
    bn[i][i] = 1;
    ep[i][i] = 1;
    en[i][i] = 1;
  }
  d1[0][1] = ld / cos(angle); // drift between dipoles
  d1[2][3] = ld / cos(angle);
  d2[0][1] = lt - 4 * lb - 2 * ld; // drift in the middle
  d2[2][3] = lt - 4 * lb - 2 * ld;
  d3[0][1] = -lt; // negative drift over total chicane to get a zero element
  d3[2][3] = -lt;

  double R = lb / sin(angle);
  double Lpath = R * angle;
  bp[2][3] = Lpath; // positive deflection angle
  bp[0][0] = cos(angle);
  bp[0][1] = R * sin(angle);
  bp[1][0] = -sin(angle) / R;
  bp[1][1] = cos(angle);

  bn[2][3] = Lpath; // negative deflection angle
  bn[0][0] = cos(-angle);
  bn[0][1] = R * sin(-angle) * -1;
  bn[1][0] = -sin(-angle) / R * -1;
  bn[1][1] = cos(-angle);

  double efoc = tan(angle) / R;
  ep[1][0] = efoc;
  ep[3][2] = -efoc;
  en[1][0] = -efoc * -1;
  en[3][2] = efoc * -1;

  this->matmul(m, bp);
  this->matmul(m, ep);
  this->matmul(m, d1);
  this->matmul(m, en);
  this->matmul(m, bn);
  this->matmul(m, d2);
  this->matmul(m, bn);
  this->matmul(m, en);
  this->matmul(m, d1);
  this->matmul(m, ep);
  this->matmul(m, bp);

  // transport matrix has been cross checked with Madx.
  /*
  cout << "lt = " << lt << " angle = " << angle << " lb = " << lb << " ld = " << ld <<  endl;
  cout << m[0][0] << " " << m[0][1] << endl;
  cout << m[1][0] << " " << m[1][1] << endl;
  cout << m[2][2] << " " << m[2][3] << endl;
  cout << m[3][2] << " " << m[3][3] << endl;
  */

  this->matmul(m, d3); // transport backwards because the main tracking still has to do the drift

  for (int i = 0; i < beam->beam.size(); i++)
  {
    for (int j = 0; j < beam->beam.at(i).size(); j++)
    {
      Particle *p = &beam->beam.at(i).at(j);
      double gammaz = sqrt(p->gamma * p->gamma - 1 - p->px * p->px - p->py * p->py); // = gamma*betaz=gamma*(1-(1+aw*aw)/gamma^2);

      double tmp = p->x;
      p->x = m[0][0] * tmp + m[0][1] * p->px / gammaz;
      p->px = m[1][0] * tmp * gammaz + m[1][1] * p->px;
      tmp = p->y;
      p->y = m[2][2] * tmp + m[2][3] * p->py / gammaz;
      p->py = m[3][2] * tmp * gammaz + m[3][3] * p->py;
    }
  }

  return;
}

void TrackBeam::matmul(double m[][4], double e[][4])
{
  double t[4][4];

  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      t[i][j] = 0;
      for (int k = 0; k < 4; k++)
      {
        t[i][j] += e[i][k] * m[k][j];
      }
    }
  }
  for (int i = 0; i < 4; i++)
  {
    for (int j = 0; j < 4; j++)
    {
      m[i][j] = t[i][j];
    }
  }
}

void TrackBeam::applyR56(Beam *beam, Undulator *und, double lambda0)
{
  double angle, lb, ld, lt;

  double gamma0 = und->getGammaRef();
  und->getChicaneParameters(&angle, &lb, &ld, &lt);
  if (angle == 0)
  {
    return;
  }
  double R56 = (4 * lb / sin(angle) * (1 - angle / tan(angle)) + 2 * ld * tan(angle) / cos(angle)) * angle;
  //    cout << "R56: " << R56 << endl;
  R56 = R56 * 4 * asin(1) / lambda0 / gamma0;
  for (int i = 0; i < beam->beam.size(); i++)
  {
    for (int j = 0; j < beam->beam.at(i).size(); j++)
    {
      beam->beam.at(i).at(j).theta += R56 * (beam->beam.at(i).at(j).gamma - gamma0);
    }
  }
  return;
}
