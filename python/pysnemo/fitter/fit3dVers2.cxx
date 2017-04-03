#include <TGraph2DErrors.h>
#include <TVirtualFitter.h>
#include <Math/Vector3D.h>
#include <TMath.h>

TGraph2DErrors* tg;
TVirtualFitter* line_fit(TGraph2DErrors* data)
{
  void myline(int &, double *, double &sum, double* par, int);

  tg = new TGraph2DErrors(data->GetN(),data->GetX(),data->GetY(),data->GetZ(),data->GetEX(),data->GetEY(),data->GetEZ());

  // using the xy,xz start value estimates from data
  double* tx=data->GetX();
  double* ty=data->GetY();
  double* tz=data->GetZ();
  double dx = tx[1] - tx[0];
  double dy = ty[1] - ty[0];
  double dz = tz[1] - tz[0];
  double xysl = dy/dx; // slopes on projections
  double xzsl = dz/dx;
  double xyi = ty[0] - tx[0]*xysl; // intercepts
  double xzi = tz[0] - tx[0]*xzsl;
  double inipar[4] = {xyi,xysl,xzi,xzsl};
  
//  TVirtualFitter::SetDefaultFitter("Minuit2"); 
  TVirtualFitter* fitter = TVirtualFitter::Fitter(0, 4);
  // make it quiet
  double pp=-1.0;
  fitter->ExecuteCommand("SET PRINTOUT", &pp, 1);
  fitter->ExecuteCommand("SET NOWARN", &pp, 1);
  // set the function
  fitter->SetFCN(myline);
  
  fitter->SetParameter(0, "intercept x-y", inipar[0],0.01,0.0,0.0);
  fitter->SetParameter(1, "slope x-y", inipar[1],0.01,0.0,0.0);
  fitter->SetParameter(2, "intercept x-z", inipar[2],0.01,0.0,0.0);
  fitter->SetParameter(3, "slope x-z", inipar[3],0.01,0.0,0.0);
  // minimize, robust first, then fine
  //  fitter->ExecuteCommand("SIMPLEX", 0, 0);
  int ierr = fitter->ExecuteCommand("MIGRAD", 0, 0);

  delete tg;

  if (ierr>0)
    return 0;
  else
    return fitter;
}

TVirtualFitter* helix_fit(TGraph2DErrors* data, double bfield)
{
  void myhelix(int &, double *, double &sum, double* par, int);

  tg = new TGraph2DErrors(data->GetN(),data->GetX(),data->GetY(),data->GetZ(),data->GetEX(),data->GetEY(),data->GetEZ());

  // start values, 6 parameter
  double* tx=data->GetX();
  double* ty=data->GetY();
  double* tz=data->GetZ();
  double dx = tx[1] - tx[0];
  double dy = ty[1] - ty[0];
  double dz = tz[1] - tz[0];
  double diag = dx*dx + dy*dy;
  double tanlambda = dz/(TMath::Sqrt(diag));
  double inipar[6] = {tx[0],ty[0],tz[0],1.0,tanlambda,bfield};
  
  TVirtualFitter::SetDefaultFitter("Minuit2"); 
  TVirtualFitter* fitter1 = TVirtualFitter::Fitter(0, 6);
  // make it quiet
  double pp=-1.0;
  fitter1->ExecuteCommand("SET PRINTOUT", &pp, 1);
  fitter1->ExecuteCommand("SET NOWARN", &pp, 1);
  // set the function
  fitter1->SetFCN(myhelix);
  
  fitter1->SetParameter(0, "xref", inipar[0],0.01,0.0,0.0);
  fitter1->SetParameter(1, "yref", inipar[1],0.01,0.0,0.0);
  fitter1->SetParameter(2, "zref", inipar[2],0.01,0.0,0.0);
  fitter1->SetParameter(3, "invrad", inipar[3],0.01,0.0,0.0);
  fitter1->SetParameter(4, "slope", inipar[4],0.01,0.0,0.0);
  fitter1->SetParameter(5, "B-field", bfield,0.01,0.0,0.0);
  fitter1->FixParameter(0);
  fitter1->FixParameter(1);
  fitter1->FixParameter(2);
  fitter1->FixParameter(5);
  // minimize, robust first, then fine
  int ierr;
  //  fitter1->ExecuteCommand("SIMPLEX", 0, 0);
  ierr = fitter1->ExecuteCommand("MIGRAD", 0, 0);
  
  double start1 = fitter1->GetParameter(3);
  double start2 = fitter1->GetParameter(4);
  delete fitter1;
  
  TVirtualFitter* fitter2 = TVirtualFitter::Fitter(0, 6);
  fitter2->ExecuteCommand("SET PRINTOUT", &pp, 1);
  fitter2->ExecuteCommand("SET NOWARN", &pp, 1);
  fitter2->SetFCN(myhelix);
  fitter2->SetParameter(0, "xref", inipar[0],0.01,0.0,0.0);
  fitter2->SetParameter(1, "yref", inipar[1],0.01,0.0,0.0);
  fitter2->SetParameter(2, "zref", inipar[2],0.01,0.0,0.0);
  fitter2->SetParameter(3, "invrad", start1,0.01,0.0,0.0);
  fitter2->SetParameter(4, "slope", start2,0.01,0.0,0.0);
  fitter2->SetParameter(5, "B-field", bfield,0.01,0.0,0.0);
  fitter2->FixParameter(5);
  ierr = fitter2->ExecuteCommand("MIGRAD", 0, 0);
  
  delete tg;
  
  if (ierr>0) {
    cout << "Error number: " << ierr << endl;
    return 0;
  }
  else
    return fitter2;
}

double distance(double x,double y,double z,double ex,double ey,double ez,double *p)
{
  // Distance point x,y,z to helix
  ROOT::Math::XYZVector weight(ex,ey,ez); 
  double w2 = weight.Mag2();
  
  double xref = p[0];
  double yref = p[1];
  double z0 = p[2];
  double omega = p[3];
  double tanlambda = p[4];
  double bfield = p[5];
  
  // from helix constructor
  if (omega==0.0) return 0.0;
  double charge = omega / TMath::Abs(omega);
  double radius = 1.0 / TMath::Abs(omega);
  double phi0 = TMath::ATan(-xref/yref);
  double xcentre = xref + radius*TMath::Cos(phi0 - TMath::PiOver2()*charge);
  double ycentre = yref + radius*TMath::Sin(phi0 - TMath::PiOver2()*charge);
  
  // from distance calc
  double sign = bfield/TMath::Abs(bfield);
  double newphi0 = TMath::ATan2(yref - ycentre, xref - xcentre);
  double phi;
  double offset;
  if (charge*sign>0.0) {
    phi = TMath::ATan2(y - ycentre, x - xcentre);
    offset = (TMath::Pi() + newphi0)/TMath::TwoPi();
  }
  else {
    phi = TMath::ATan2(y - ycentre, x - xcentre) + TMath::Pi();
    offset = (newphi0)/TMath::TwoPi() - 0.5;
  }    
  int nCircles;
  double xCirc = offset + (sign*charge*(z+z0))/(TMath::TwoPi()*tanlambda*radius);
  if (xCirc==1.0)
    nCircles = (int)(xCirc-1.0);
  else
    nCircles = (int)xCirc;

  double dphi;  
  if (charge*sign>0.0) 
    dphi = TMath::TwoPi()*(Double_t)nCircles + phi - newphi0;
  else
    dphi = TMath::TwoPi()*(Double_t)nCircles + phi - newphi0 - TMath::Pi();
    
  double zonhelix = z0 + sign*charge*radius*tanlambda*dphi;
  double distz =  TMath::Abs(zonhelix - z);
  double distxy = TMath::Sqrt((xcentre-x)*(xcentre-x) + (ycentre-y)*(ycentre-y));
  distxy = TMath::Abs(distxy - radius);  
  double d2 = distxy*distxy + distz*distz;
  return d2/w2; 
}


double distance2(double x,double y,double z,double ex,double ey,double ez,double *p)
{
  // assume never parallel to z,y plane = foil plane
  ROOT::Math::XYZVector xp(x,y,z); 
  ROOT::Math::XYZVector weight(ex,ey,ez); 
  ROOT::Math::XYZVector x0(0. , p[0], p[2]); 
  ROOT::Math::XYZVector x1(1. , p[0] + p[1], p[2] + p[3]); 
  ROOT::Math::XYZVector u = (x1-x0).Unit(); 
  double d2 = ((xp-x0).Cross(u)).Mag2(); 
  double w2 = weight.Mag2();
  return d2/w2; 
}

// FCN interface function, fixed by MIGRAD, only 'sum' and 'par' needed
void myline(int &, double *, double &sum, double* par, int)
{
  double distance2(double x,double y,double z,double ex,double ey,double ez,double *p);
  // tg is the TGraph2DErrors object; global here
  double * x = tg->GetX();
  double * y = tg->GetY();
  double * z = tg->GetZ();
  double * ex = tg->GetEX();
  double * ey = tg->GetEY();
  double * ez = tg->GetEZ();
  int npoints = tg->GetN();
  sum = 0.0;
  for (int i  = 0; i < npoints; ++i) { 
    sum += distance2(x[i],y[i],z[i],ex[i],ey[i],ez[i],par); ;
  }
}


// FCN interface function, fixed by MIGRAD, only 'sum' and 'par' needed
void myhelix(int &, double *, double &sum, double* par, int)
{
  double distance(double x,double y,double z,double ex,double ey,double ez,double *p);
  // tg is the TGraph2DErrors object; global here
  double * x = tg->GetX();
  double * y = tg->GetY();
  double * z = tg->GetZ();
  double * ex = tg->GetEX();
  double * ey = tg->GetEY();
  double * ez = tg->GetEZ();
  int npoints = tg->GetN();
  sum = 0.0;
  for (int i = 0; i < npoints; ++i) { 
    sum += distance(x[i],y[i],z[i],ex[i],ey[i],ez[i],par); ;
  }
}
