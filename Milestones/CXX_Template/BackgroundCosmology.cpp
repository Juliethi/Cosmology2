#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaLambda,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaLambda(OmegaLambda),
  Neff(Neff), 
  TCMB(TCMB)
{
  H0 = Constants.H0_over_h*h; //SI units
  OmegaR = 2*M_PI*M_PI/30 * pow(Constants.k_b*TCMB, 4)/(pow(Constants.hbar,3)*pow(Constants.c,5)) * 8*M_PI*Constants.G/(3*H0*H0);
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
  int x_start = -20;
  int x_end = 5;
  int n = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, n);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){
    detadx[0] = Constants.c/Hp_of_x(x);
    return GSL_SUCCESS;
  };

  Vector eta_ic = {0.0};
  ODESolver ode;
  ode.solve(detadx, x_array, eta_ic);
  auto all_data = ode.get_data();
  Vector eta = Vector(n);
  
  for(int i = 0; i < all_data.size(); i++){
    eta[i] = all_data[i][0];
  }

  eta_of_x_spline.create(x_array,eta);

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  return H0 * sqrt((OmegaB + OmegaCDM)*exp(-3*x) + OmegaR*exp(-4*x) + OmegaLambda);
}

double BackgroundCosmology::Hp_of_x(double x) const{
  return H_of_x(x) * exp(x);
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
  double A = (OmegaB + OmegaCDM)*exp(-x) + 2*OmegaR*exp(-2*x) - 2*OmegaLambda*exp(2*x);
  double B = (OmegaB + OmegaCDM)*exp(-x) + OmegaR*exp(-2*x) + OmegaLambda*exp(2*x);
  return -H0/2 * A/sqrt(B);
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;
  double H_over_H0 = H_of_x(x)/H0;
  return OmegaB/(exp(3*x)*H_over_H0*H_over_H0);
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;
  double H_over_H0 = H_of_x(x)/H0;
  return OmegaR/(exp(4*x)*H_over_H0*H_over_H0);
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  return 0.0;
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;
  double H_over_H0 = H_of_x(x)/H0;
  return OmegaCDM/(exp(3*x)*H_over_H0*H_over_H0);
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;
  double H_over_H0 = H_of_x(x)/H0;
  return OmegaLambda/(H_over_H0*H_over_H0);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
  //return 0;
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB() const{ 
  return TCMB; 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -16.0;
  const double x_max =  2.0;
  const int    n_pts =  1000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

