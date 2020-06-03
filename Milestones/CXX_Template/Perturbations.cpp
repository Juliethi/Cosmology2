#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================

  //Using the quadratic spacing found in Callin
  //Vector k_array(n_k);
  /*
  for(int i = 0; i < n_k; i++){
    double c = i/100.;
    k_array[i] = k_min + (k_max - k_min)*c*c;
  }
  */

  
  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);
  
  
  Vector x_array = Utils::linspace(x_start,x_end,n_x);

  Vector2D kx_array(n_x, Vector(n_k));
  Vector3D all_solutions(Constants.n_ell_tot_full, kx_array);

  Vector2D Psi_vector(n_x, Vector(n_k));




  //pluss1 because we also want to include Psi
  //all_solutions[y_i][k][x]


  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){
    //std::cout << "ik=" << ik << std::endl;
    // Progress bar...
    //std::cout << "ik" << ik << std::endl;

    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];
    //std::cout << "sug pikk 1" << std::endl;

    // Find value to integrate to
    double x_end_tight_temp = get_tight_coupling_time(k);


    //finding what index that corresponds the closest to x_end_tight to ensure that everything matches with the full x-array 
    int index_x_end_tight = (int) round((x_end_tight_temp-x_start)/(x_end-x_start)*(n_x-1));
    double x_end_tight = x_array[index_x_end_tight];
    int n_x_tc = index_x_end_tight+1;
  
    Vector x_array_tc = Utils::linspace(x_start, x_end_tight, n_x_tc);
    Vector x_array_after_tc = Utils::linspace(x_end_tight, x_end, n_x - index_x_end_tight);




    //===================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);


    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };


   
    ODESolver ode_tc;
    ode_tc.solve(dydx_tight_coupling, x_array_tc, y_tight_coupling_ini);


    //Looping over all components, adding them to the full-solutions-vector   
    for(int i = 0; i < Constants.n_ell_tot_tc; i++){
      auto solution_i = ode_tc.get_data_by_component(i);
      for(int j = 0; j < n_x_tc; j++){
        all_solutions[i][j][ik] = solution_i[j];
      }
    }

    //====i===============================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    // auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    //SETTING INITIAL CONDITIONS AFTER TIGHT COUPLING
    Vector y_tc_last_value(Constants.n_ell_tot_full); 
    
    for(int i = 0; i<Constants.n_ell_tot_tc; i++){                   
      y_tc_last_value[i] = all_solutions[i][index_x_end_tight][ik];
    }

    const double Hp = cosmo->Hp_of_x(x_end_tight);
    const double ck_over_Hp = Constants.c*k/Hp;
    const double dtaudx = rec->dtaudx_of_x(x_end_tight);

    double Theta1 = y_tc_last_value[Constants.ind_start_theta_tc+1];
    double Theta2 = -20./(45*dtaudx)*ck_over_Hp*Theta1;
    y_tc_last_value[Constants.ind_start_theta_tc+2] = Theta2;


    for(int l = 3; l<Constants.n_ell_theta; l++){ 
      double theta_l = -l/(2.*l+1)*ck_over_Hp*y_tc_last_value[l-1];
      y_tc_last_value[Constants.ind_start_theta + l]= theta_l; 
    }


    // The full ODE system
    
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    ODESolver ode_full;
    ode_full.solve(dydx_full, x_array_after_tc, y_tc_last_value);
    
    for(int i = 0; i < Constants.n_ell_tot_full; i++){
      auto solution_i = ode_full.get_data_by_component(i);
      for(int j = n_x_tc; j < n_x; j++){
        //std::cout << j << std::endl;
        all_solutions[i][j][ik] = solution_i[j-index_x_end_tight-1];
      }
    }

    //Calculating Psi
    for(int ix = 0; ix<n_x; ix++){
      
      double x = x_array[ix];
      double Hp = cosmo->Hp_of_x(x);
      double H0 = cosmo->get_H0();
      double OmegaR = cosmo->get_OmegaR(0);
      double dtaudx = rec->dtaudx_of_x(x);

      const double ck_over_Hp = Constants.c*k/Hp;
      //theta0 index = ind_start_theta, theta_l index = ind_start_theta + l
      double Theta1 = all_solutions[Constants.ind_start_theta+1][ix][ik];
      //double Theta2 = -20./(45*dtaudx)*ck_over_Hp*Theta1;
      double Theta2 = all_solutions[Constants.ind_start_theta+2][ix][ik];
      double Phi = all_solutions[Constants.ind_Phi][ix][ik];
      //Solving psi 
      Psi_vector[ix][ik] = -Phi - 12.*H0*H0/(Constants.c*Constants.c*k*k*exp(2*x))*OmegaR*Theta2;
      //Psi_vector[ix][ik] = -Phi;
      
    }


    //===================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    //===================================================================
    //...
    //...

  }
  Utils::StartTiming("integrateperturbation");


 Vector2D all_solutions_flattened(Constants.n_ell_tot_full, Vector(n_x*n_k));
 Vector Psi_flattened(n_x*n_k);

 for(int ix = 0; ix<n_x; ix++){
   for(int ik = 0; ik<n_k; ik++){
     Psi_flattened[ix+n_x*ik] = Psi_vector[ix][ik];
   }
 }



  for(int i=0; i<Constants.n_ell_tot_full; i++){
    for(int ix=0; ix<n_x; ix++){
      for(int ik=0; ik < n_k; ik++){
        all_solutions_flattened[i][ix+n_x*ik] = all_solutions[i][ix][ik];
      }
    }  
  }
  

  //I feel like this shoild be possible to do in a loop what am i missing here grrrr

  
  
  delta_cdm_spline.create(x_array, k_array, all_solutions_flattened[Constants.ind_deltacdm], "delta_cdm_spline"); //0
  delta_b_spline.create(x_array, k_array, all_solutions_flattened[Constants.ind_deltab],"delta_b_spline" ); //1
  v_cdm_spline.create(x_array, k_array, all_solutions_flattened[Constants.ind_vcdm], "v_cdm_spline"); //2
  v_b_spline.create(x_array, k_array,all_solutions_flattened[Constants.ind_vb], "v_b_spline"); //3
  Phi_spline.create(x_array,k_array, all_solutions_flattened[Constants.ind_Phi], "Phi_spline"); //4
  Psi_spline.create(x_array, k_array, Psi_flattened, "Psi_spline" ); 

  Theta_0_spline.create(x_array, k_array, all_solutions_flattened[Constants.ind_start_theta_tc], "Theta_0_spline"); //5
  Theta_1_spline.create(x_array, k_array, all_solutions_flattened[Constants.ind_start_theta_tc+1], "Theta_1_spline"); //6
  Theta_2_spline.create(x_array, k_array, all_solutions_flattened[Constants.ind_start_theta_tc+2], "Theta_2_spline"); //6

  
  
  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  // ...
  // ...
  // ...
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;


  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  //double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  
  const double Hp = cosmo->Hp_of_x(x);
  const double Psi_ic = -2.0/3.0;
  const double ck = Constants.c*k;

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  
  Phi = -Psi_ic;

  delta_cdm = -3.0/2.0*Psi_ic;
  delta_b = delta_cdm;
  v_cdm = -ck/(2.0*Hp)*Psi_ic;

  v_b = v_cdm;
  
  //std::cout << v_cdm << " " << k*Constants.Mpc << std::endl;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = -0.5*Psi_ic;
  Theta[1] = ck/(6.*Hp)*Psi_ic;
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  //const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  //double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  //double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  // ...
  // ...
  // ...

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  // ...
  // ...

  // SET: Photon temperature perturbations (Theta_ell)
  // ...
  // ...



  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;
  //getting spline objects from recombination class to utilize the binary search for value
  Spline log_Xe_spline = rec->log_Xe_of_x_spline;
  Spline dtaudx_spline = rec->dtaudx_of_x_spline;

  //Condition one: tight coupling needs to end before recombination. For simplicity, I say recombintaion starts when Xe=0.99
  double x1 = Utils::binary_search_for_value(log_Xe_spline, log(0.99), {x_start, x_end},1e-5);
  
  //Condition two: Since dtaudx is strictly negative, and im not sure how to take the absoulte value of a spline object but would very much
  //like to keep using the binary search for value function, i check when dtaudx=-10 instead, because that should be the same as abs(dtaudx)=10 in this case 
  
  double x2 = Utils::binary_search_for_value(dtaudx_spline, -10, {x_start,x_end}, 1e-5);

  //I assume the third condition is close to the first two
  double x_high = std::min(x1,x2)+2;
  double x_low = std::max(x1,x2)-2;
  double npts = 1000;
  Vector x_search = Utils::linspace(x_low,x_end,npts);
  double x3;

  for(int i=0; i<npts; i++){
    double x = x_search[i];
    double const Hp = cosmo->Hp_of_x(x);
    double const dtau = rec->dtaudx_of_x(x);
    double abs_dtau = abs(dtau);    
    double ck_over_Hp = Constants.c*k/Hp;
    if(abs_dtau < 10*ck_over_Hp){
      //std::cout << "X3 FOUND TO BE=" << x << std::endl;
      x3 = x;
      break;
      }
    if(i == npts-1){
      std::cout << "No x3 found in given interval" << std::endl;
      x3 = 1; 
    }
  }

  double x_min1 = std::min(x1,x2);
  double x_min2 = std::min(x3,x_min1);
  x_tight_coupling_end = x_min2;
  //std::cout << "x tight coulping end = " << x_tight_coupling_end << std::endl;


  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // TODO: Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // ...
  // ...
  //Vector k_array;
  //Vector x_array;


  Vector log_k_array = Utils::linspace(log(k_min), log(k_max), n_k);
  Vector k_array = exp(log_k_array);
  Vector x_array = Utils::linspace(x_start,x_end,n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    //for(auto ix = 0; ix < 2; ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      //for(auto ik = 0; ik < 2; ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      const double Hp = cosmo->Hp_of_x(x);
      const double dhpdx = cosmo->dHpdx_of_x(x);
      const double ddhpdx = cosmo->ddHpddx_of_x(x);

      const double g = rec->g_tilde_of_x(x);
      const double dgdx = rec->dgdx_tilde_of_x(x);
      const double ddgdx = rec->ddgddx_tilde_of_x(x);
      const double tau = rec->tau_of_x(x);

      const double theta0 = get_Theta(x, k, 0);
      const double theta2 = get_Theta(x, k, 2);
      const double dtheta2dx = Theta_2_spline.deriv_x(x,k);
      const double ddtheta2dx = Theta_2_spline.deriv_xx(x,k);

      const double Psi = get_Psi(x,k);
      const double Phi = get_Phi(x,k);

      const double dPsidx = Psi_spline.deriv_x(x,k);
      const double dPhidx = Phi_spline.deriv_x(x,k);

      const double vb = get_v_b(x, k);
      const double dvbdx = v_b_spline.deriv_x(x,k);

      const double c = Constants.c;
      double ck = c*k;

      double A = g*(theta0 + Psi + 1./4.*theta2);
      double B = exp(-tau)*(dPsidx - dPhidx);
      double C = 1./ck*(dhpdx*g*vb + Hp*dgdx*vb + Hp*g*dvbdx);
      double D = 3./(4.*ck*ck) * (dhpdx*dhpdx + Hp*ddhpdx)*g*theta2 + 3.*Hp*dhpdx*(dgdx*theta2 + g*dtheta2dx) + Hp*Hp*(ddgdx*theta2 + 2.*dgdx*dtheta2dx + g*ddtheta2dx);
      // Temperatur source
      ST_array[index] = A + B - C + D;
      /*
      std::cout << "A= " << A << std::endl;
      std::cout << "B= " << B << std::endl;
      std::cout << "C= " << C << std::endl;
      std::cout << "D= " << D << std::endl;
      */
      
      //std::cout << ST_array[0] << std::endl;
      
      

    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");


  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  //const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  //double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================
  
  const double dtaudx = rec-> dtaudx_of_x(x);
  const double ddtaudx = rec-> ddtauddx_of_x(x);
  const double Hp = cosmo->Hp_of_x(x);
  const double dHpdx = cosmo->dHpdx_of_x(x);
  const double Omega_r = cosmo->get_OmegaR(0);
  const double Omega_b = cosmo->get_OmegaB(0);
  const double Omega_cdm = cosmo->get_OmegaCDM(0);
  const double H0 = cosmo->get_H0();
  const double R = 4.*Omega_r/(3*Omega_b*exp(x));
  const double ck_over_Hp = Constants.c*k/Hp;
  const double dHpdx_over_Hp = dHpdx/Hp;

  const double Theta2 = -20./(45*dtaudx)*ck_over_Hp*Theta[1];
  const double Psi = -Phi - 12.*H0*H0/(Constants.c*Constants.c*k*k*exp(2*x))*Omega_r*Theta2;

  dPhidx = Psi - ck_over_Hp*ck_over_Hp*1./3.*Phi + H0*H0/(2.*Hp*Hp)*(Omega_cdm*exp(-x)*delta_cdm + Omega_b*exp(-x)*delta_b + 4.*Omega_r*exp(-2*x)*Theta[0]);
  dThetadx[0] = -ck_over_Hp*Theta[1]- dPhidx;  
  
  double q = (-((1-R)*dtaudx + (1+R)*ddtaudx)*(3.*Theta[1]+v_b) - ck_over_Hp*Psi + (1-dHpdx_over_Hp)*ck_over_Hp*(-Theta[0] + 2.*Theta2) - ck_over_Hp*dThetadx[0])/((1+R)*dtaudx + dHpdx_over_Hp - 1);
  
  ddelta_cdmdx = ck_over_Hp*v_cdm - 3.*dPhidx;
  ddelta_bdx = ck_over_Hp*v_b - 3.*dPhidx;
  dv_cdmdx = -v_cdm - ck_over_Hp*Psi;
  dv_bdx = 1./(1+R)*(-v_b - ck_over_Hp*Psi + R*(q + ck_over_Hp*(-Theta[0] + 2*Theta2)-ck_over_Hp*Psi));


  dThetadx[1] = 1./3.*(q-dv_bdx);

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  //const double *Theta_p         = &y[Constants.ind_start_thetap];
  //const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  //double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  //double *dNudx           = &dydx[Constants.ind_start_nu];

  
  const double dtaudx = rec-> dtaudx_of_x(x);
  const double ddtaudx = rec-> ddtauddx_of_x(x);

  const double Hp = cosmo->Hp_of_x(x);
  const double dHpdx = cosmo->dHpdx_of_x(x);
  const double eta = cosmo->eta_of_x(x);

  const double Omega_r = cosmo->get_OmegaR(0);
  const double Omega_b = cosmo->get_OmegaB(0);
  const double Omega_cdm = cosmo->get_OmegaCDM(0);
  const double H0 = cosmo->get_H0();
 
  const double R = 4*Omega_r/(3*Omega_b*exp(x));
  const double ck_over_Hp = Constants.c*k/Hp;
  const double dHpdx_over_Hp = dHpdx/Hp;

  double Psi = -Phi - 12.0*H0*H0/(Constants.c*Constants.c*k*k*exp(2*x))*Omega_r*Theta[2];

  dThetadx[0] = -ck_over_Hp*Theta[1]- dPhidx;
  dThetadx[1] = 1./3.*ck_over_Hp*Theta[0] - 2./3.*ck_over_Hp*Theta[2] + 1./3.*ck_over_Hp*Psi + dtaudx*(Theta[1] + 1./3*v_b);  
  
  for(int l=2; l<Constants.n_ell_theta-1; l++){
    double kd = (l==2) ? 1 : 0;
    dThetadx[l] = l/(2.*(l+1))*ck_over_Hp*Theta[l-1] - (l+1)/(2.*l+1)*ck_over_Hp*Theta[l]+dtaudx*(Theta[l] - 1/10.*Theta[2]*kd);
  }

  int l_max = Constants.n_ell_theta-1;
  dThetadx[l_max] = ck_over_Hp*Theta[l_max-1] - Constants.c*(l_max+1.)/(Hp*eta)*Theta[l_max]+dtaudx*Theta[l_max];
  
  dPhidx = Psi - ck_over_Hp*ck_over_Hp*1./3.*Phi + H0*H0/(2.*Hp*Hp)*(Omega_cdm*exp(-x)*delta_cdm + Omega_b*exp(-x)*delta_b + 4.*Omega_r*exp(-2*x)*Theta[0]);

  ddelta_cdmdx = ck_over_Hp*v_cdm - 3.*dPhidx;
  ddelta_bdx = ck_over_Hp*v_b - 3.*dPhidx;
  dv_cdmdx = -v_cdm - ck_over_Hp*Psi;
  dv_bdx = -v_b - ck_over_Hp*Psi +dtaudx*R*(3.*Theta[1] + v_b);


  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  if(ell==0)
    return Theta_0_spline(x,k);
  if(ell==1)
    return Theta_1_spline(x,k);
  if(ell==2)
    return Theta_2_spline(x,k);
  if(ell>2)
    std::cout << "no valid theta_ell value, abort mission" << std::endl;
    return -1.0;
  //return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * Constants.c * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    fp << get_v_b(x,k)       << " ";
    fp << get_v_cdm(x,k)     << " ";
    fp << get_delta_b(x,k)   << " ";
    fp << get_delta_cdm(x,k)   << " ";
    //fp << get_Pi(x,k)        << " ";
    fp << get_Source_T(x,k)  << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

