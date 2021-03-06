#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================
  /*
  Vector k_array;
  Vector log_k_array = log(k_array);
  */
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);

  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();

  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  
  //=========================================================================
  // TODO: Do the same for polarization...
  //=========================================================================
  // ...
  // ...
  // ...
  // ...
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  int npts = 1e4;
  Vector log_z_array = Utils::linspace(log(1e-8), log(4e4), npts+1); 
  Vector z_array = exp(log_z_array);
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    Vector j_ell_array(npts+1);
    for(int j=0; j<npts+1; j++){
      double z = z_array[j];
      j_ell_array[j] = Utils::j_ell(ell, z);
    }

    Spline j_ell_spline_i;
    j_ell_spline_i.create(log_z_array,j_ell_array, "j_ell_spline_i");
    j_ell_splines[i] = j_ell_spline_i;

    // Make the j_ell_splines[i] spline
  }

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  
  int npts = 2e3;
  Vector x_array = Utils::linspace(x_start, x_end, npts+1);

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t ik = 0; ik < k_array.size(); ik++){
    double k = k_array[ik];

    std::cout << "ik= " << ik << std::endl;

    for(size_t iell = 0; iell < ells.size(); iell++){
      double integral = 0;
      int ell = ells[iell];
      Spline j_ell_spline = j_ell_splines[iell];
      ODEFunction dThetadx = [&](double x, const double *Theta_ell, double *dThetadx){
        double source_function = pert->get_Source_T(x, k);
        double eta = cosmo->eta_of_x(x);
        double eta_0 = cosmo->eta_of_x(0);
        double bessel_term = k*(eta_0 - eta);

        if(bessel_term < 1e-8){
          bessel_term = 1e-8;
        }
        double j_ell = j_ell_spline(log(bessel_term));
        dThetadx[0] = source_function*j_ell;
      return GSL_SUCCESS;
      };

      Vector Theta_ic = {0.0};
      ODESolver ode;
      double hstart = 1e-3, abserr = 1e-6, relerr = 1e-6;
      ode.set_accuracy(hstart, abserr, relerr);
      ode.solve(dThetadx, x_array, Theta_ic);
      result[iell][ik] = ode.get_data_by_component(0)[npts];
    }
    

    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    // ...
    // ...
    // ...

    // Store the result for Source_ell(k) in results[ell][ik]
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);


  for(size_t iell = 0; iell < ells.size(); iell++){
    Spline theta_iell_spline;
    theta_iell_spline.create(k_array, thetaT_ell_of_k[iell], "theta_iell_spline");
    thetaT_ell_of_k_spline[iell] = theta_iell_spline;
  }

  // Spline the result and store it in thetaT_ell_of_k_spline
  // ...
  // ...
  // ...
  // ...

  //============================================================================
  // TODO: Solve for ThetaE_ell(k) and spline
  //============================================================================

}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();

  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  Vector result(nells);
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_t iell = 0; iell < nells; iell++){

    ODEFunction dCelldlogk = [&](double log_k, const double *C_ell, double *dCelldlogk){
      double k = exp(log_k);
      //double k_pivot = kpivot_mpc/Constants.Mpc;
      double f_ell = f_ell_spline[iell](k);
      double g_ell = g_ell_spline[iell](k);
      double P = primordial_power_spectrum(k);
    
      dCelldlogk[0] = 4*M_PI*P*f_ell*g_ell;
      return GSL_SUCCESS;
      };
    
    Vector C_ell_ini{0.0};
    ODESolver ode;
    double hstart = 1e-3, abserr = 1e-6, relerr = 1e-6;
    ode.set_accuracy(hstart, abserr, relerr);
    ode.solve(dCelldlogk, log_k_array, C_ell_ini);
    result[iell] = ode.get_data_by_component(0)[n_k-1];
  }
  

  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double k = k_mpc;
  double c = Constants.c;
  double Psi = pert->get_Psi(x,k_mpc);
  double Omega_M = cosmo->get_OmegaB(0) + cosmo->get_OmegaCDM(0);
  double H0 = cosmo->get_H0();
  double C = 2*M_PI*M_PI/pow(k,3);
  double delta_M = (2.*c*c*k*k*Psi*exp(x))/(3.*Omega_M*H0*H0);
  double pofk = abs(delta_M*delta_M)*primordial_power_spectrum(k)*C;

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}

//====================================================
// Output functions
//====================================================

void PowerSpectrum::output_integrand_theta_ell(std::string filename) const{
  std::ofstream fp(filename.c_str());
  double npts = 1e4;
  Vector x_array = Utils::linspace(x_start, x_end, npts+1);
  for(int ix = 0; ix < npts; ix++){
      double x = x_array[ix];
      double k = 340*cosmo->get_H0()/Constants.c;
      double source_function = pert->get_Source_T(x, k);
      double eta = cosmo->eta_of_x(x);
      double eta_0 = cosmo->eta_of_x(0);
      double j_ell = j_ell_splines[19](log(k*(eta_0 - eta))); //iell = 19 -> ell=100
      double integrand = j_ell*source_function;
      fp << x << " ";
      fp << integrand << " ";
      fp << "\n"; 
  }
}

void PowerSpectrum::output_matter_pk(std::string filename) const{
  std::ofstream fp(filename.c_str());
  
  double x = 0;
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);
  //double k_peak = cosmo->get_H0()*(cosmo->get_OmegaB(0) + cosmo->get_OmegaCDM(0))*sqrt(2/cosmo->get_OmegaR(0));
  double x_eq = -8.57; //from milestone 1
  double k_peak = cosmo->Hp_of_x(x_eq)/Constants.c;
  fp << k_peak << " ";
  fp << 0 << " "; //only to make it easier to read the file in python afterwards
  fp << "\n";

  for(int ik=0;ik < n_k; ik++){
    double k = k_array[ik]; ///Constants.Mpc;
    double pofk = get_matter_power_spectrum(x, k);
    fp << k << " ";
    fp << pofk << " ";
    fp <<"\n";
  }
}

void PowerSpectrum::output_theta_l(const int iell, std::string filename) const{
  std::ofstream fp(filename.c_str());
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);
  Spline theta_spline = thetaT_ell_of_k_spline[iell];
  for(int ik = 0; ik < n_k; ik++){
    double k = k_array[ik];
    double theta = theta_spline(k);
    fp << k << " ";
    fp << theta << " ";
    fp <<"\n";
  }
  


}

//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";
    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
}

