#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp,
    double h,
    double TCMB) :
  cosmo(cosmo),
  Yp(Yp),
  h(h),
  TCMB(TCMB)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  //int x_start = -20;
  //int x_end = 0;
  //int n = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_array(npts_rec_arrays);
  Vector ne_array(npts_rec_arrays);

  //Vector x_array;
  //Vector Xe_arr;
  //Vector ne_arr;

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      //std::cout << "Xe current=" << Xe_current << std::endl;
      Xe_array[i] = Xe_current;
      ne_array[i] = ne_current;

    } else {

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      double x_start_peebles = x_array[i];
      Vector X_peebles = Utils::linspace(x_start_peebles, x_end, npts_rec_arrays-i);

      //Not sure if I want to use linspace or for loop here, both seem to do the same but... linspace feels dangerous, too python-y. It looks good tho:)

      //Vector X_peebles(npts_rec_arrays-i);
      //for(int j=0; j<npts_rec_arrays-i; j++){
      //  X_peebles[j] = x_array[j+i];
      //}

      //double hm = X_peebles.size();
      //std::cout << "size of x_peebles=" << hm << std::endl;
      Vector peebles_ic ={Xe_current};

      ODESolver ode;
      ode.solve(dXedx, X_peebles, peebles_ic);
      auto Xe_array_peebles = ode.get_data_by_component(0);
      
      const double OmegaB = cosmo->get_OmegaB();
      const double H0_over_h   = Constants.H0_over_h;
      const double H0 = H0_over_h*h;
      double rho_c = 3*H0*H0/(8*M_PI*Constants.G);

      for(int j = 0; j < npts_rec_arrays-i; j++){
        Xe_array[j+i] = Xe_array_peebles[j];
        //std::cout << Xe_array_peebles[j] << std::endl;

        double n_b = OmegaB*rho_c/Constants.m_H*exp(-3*X_peebles[j]);
        ne_array[j+i] = Xe_array_peebles[j]*n_b;
      }

      break;
    
    }
  }
  //Creating splines of log(Xe) and log(ne). The functions Xe_of_x and ne_of_x will... unlog this further down.



  Vector logXe = log(Xe_array);
  Vector logne = log(ne_array);
  log_Xe_of_x_spline.create(x_array, logXe, "log_Xe_spline");
  log_ne_of_x_spline.create(x_array, logne);



  //Solving Xe for ONLY Saha to compare this with full solution 
  Vector Xe_array_saha_only(npts_rec_arrays);
  for(int i = 0; i < npts_rec_arrays; i++){

    auto Xe_ne_data_Saha = electron_fraction_from_saha_equation(x_array[i]);
    const double Xe_current_Saha = Xe_ne_data_Saha.first;
    Xe_array_saha_only[i] = Xe_current_Saha;
    //std::cout << Xe_current_Saha << std::endl;
  }
  Vector logXe_Saha_only = log(Xe_array_saha_only);
  log_Xe_of_x_Saha_only_Spline.create(x_array,logXe_Saha_only);



  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;
  const double H0 = H0_over_h*h; //SI UNITS

  // Fetch cosmological parameters
  //const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaB = cosmo->get_OmegaB(); //todays value of omega_B
  double rho_c = 3*H0*H0/(8*M_PI*G);
  

  // Electron fraction and number density
  //double Xe = 0.0;
  //double ne = 0.0;
  double n_b = OmegaB*rho_c/(m_H*pow(a,3)); //also n_h
  double T_b = TCMB/a;
  double C = 1/(n_b*pow(hbar,3))*pow(m_e*T_b*k_b/(2*M_PI), 3.0/2.0)*exp(-epsilon_0/(k_b*T_b));
  //std::cout << "C=" << C << std::endl;
  double Xe;

  //End points. If C is very large(it is hot), a is very small? => super duper early universe => X_e = 1
  if(C > 1e6){
    //std::cout << "henlo i am big" << std::endl;
    Xe = 1;
  } 

  //If C is very small, the Saha equation approaches 0, to mitigate instabilites we set it to something small. we could set it to 0 but that would create problems later when taking the log
  if(C<1e-20){
    //std::cout << "henlo i am smol" << std::endl;
    Xe = 1e-50;
  }

  else{
    //std::cout << "henlo i am just right" << std::endl;
    Xe = -C/2.0 + sqrt(C*(C+4))/2.0;
  }
  
  //std::cout << "Xe=" << Xe << std::endl;

  double ne = Xe*n_b;

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){
  
  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;
  //const double H0 = H0_over_h*h;

  const double OmegaB = cosmo->get_OmegaB(); //todays value of omega_B
  const double H = cosmo -> H_of_x(x);
  const double H0 = cosmo->get_H0();
  
  double rho_c = 3*H0*H0/(8*M_PI*G);
  
  double n_b = OmegaB*rho_c/(m_H*pow(a,3)); //also n_h
  double T_b = TCMB/a;

  double n_1s = (1-X_e)*n_b;
  double Lambda_alpha = H*pow(3*epsilon_0,3)/(8*8*M_PI*M_PI*n_1s)*pow(hbar*c,-3);
  double Lambda_2s1s = 8.227;
  double phi2 = 0.448*log(epsilon_0/(k_b*T_b));
  double alpha2 = 8./sqrt(3*M_PI)*c*sigma_T*sqrt(epsilon_0/(k_b*T_b))*phi2;
  double beta = alpha2/pow(hbar,3)*pow(m_e*k_b*T_b/(2*M_PI),3./2.)*exp(-epsilon_0/(k_b*T_b));
  double beta2 = alpha2/pow(hbar,3)*pow(m_e*k_b*T_b/(2*M_PI),3./2.)*exp(-epsilon_0/(4.*k_b*T_b));

  double Cr = (Lambda_2s1s + Lambda_alpha)/(Lambda_2s1s+Lambda_alpha + beta2);
  
  dXedx[0] = Cr/H*(beta*(1-X_e) - n_b*alpha2*X_e*X_e);
  

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = npts_rec_arrays;
  //double x_end = 2;
  double x_start = -12; 
  //double where_x0 = 0.0; 
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){
    double n_e = ne_of_x(x);
    double H = cosmo->H_of_x(x);
    dtaudx[0] = - Constants.c*n_e*Constants.sigma_T/H;
    return GSL_SUCCESS;
  };


  //we dont know what the intial conditions of tau is, we only know that it should be 0 in the last point, so we set it to some arbitrary value now and scale it later
  Vector tau_ic = {1e3}; 

  ODESolver ode;
  ode.solve(dtaudx, x_array, tau_ic);
  auto tau_unscaled = ode.get_data_by_component(0);
  auto dtaudx_array = ode.get_derivative_data_by_component(0);
  
  Vector tau_scaled(npts_rec_arrays);
  
  //for(int i = 0; i < npts_rec_arrays; i++){
  //  std::cout << tau_unscaled[i] << std::endl;
  //}

  //temporary
  //tau_of_x_spline.create(x_array,x_array); 
  


  double scale_difference = tau_unscaled[npts_rec_arrays-1]; //the value of tau unscaled when x=0
  //std::cout << scale_difference << " " << x_array[npts_rec_arrays-1] << std::endl;

  //now we need to rescale: 
  for(int i = 0; i < npts_rec_arrays; i++){
    tau_scaled[i] = tau_unscaled[i] - scale_difference;
    //std::cout << tau_scaled[i] << std::endl;
  }

  tau_of_x_spline.create(x_array,tau_scaled);
  dtaudx_of_x_spline.create(x_array,dtaudx_array, "dtaudx");
  //ddtaudx is created in the function? 


  Vector g_tilde_array(npts_rec_arrays);
  for(int i; i < npts_rec_arrays; i++){
    g_tilde_array[i] = -dtaudx_array[i]*exp(-tau_scaled[i]);
    //std::cout << g_tilde_array[i] << std::endl;
  }

  g_tilde_of_x_spline.create(x_array,g_tilde_array);


  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  return dtaudx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

double RecombinationHistory::Xe_of_x_Saha_only(double x) const{
  return exp(log_Xe_of_x_Saha_only_Spline(x));
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = -12;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << Xe_of_x_Saha_only(x) << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

