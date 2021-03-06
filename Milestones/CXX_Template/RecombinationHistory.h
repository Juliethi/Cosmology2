#ifndef _RECOMBINATION_HISTORY_HEADER
#define _RECOMBINATION_HISTORY_HEADER
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Utils.h"
#include "BackgroundCosmology.h"

using Vector = std::vector<double>;

class RecombinationHistory{
  private:

    // The cosmology we use
    BackgroundCosmology *cosmo = nullptr;
    
    // Helium fraction
    double Yp;
    double h;
    double TCMB;
 
    // The start and end points for recombination arrays (can be modified)
    const double x_start  = Constants.x_start;
    const double x_end    = Constants.x_end;
    
    // Numbers of points of Xe,ne array (modify as you see fit)
    const int npts_rec_arrays = 4000;
  
    // Xe for when to switch between Saha and Peebles
    const double Xe_saha_limit = 0.99;

    //===============================================================
    // [1] Computation of Xe (Saha and Peebles equation)
    //===============================================================
 
    // Compute Xe from the Saha equation
    std::pair<double,double> electron_fraction_from_saha_equation(double x) const;
    
    // Right hand side of the dXedx Peebles equation
    int rhs_peebles_ode(double x, const double *y, double *dydx);
    
    // Solve for Xe 
    void solve_number_density_electrons();
    
    //===============================================================
    // [2] Compute tau and visibility functions
    //===============================================================

    // The two things we need to solve: Xe/ne and tau
    void solve_for_optical_depth_tau();

    // Splines contained in this class
    //Spline log_Xe_of_x_spline{"Xe"};
    Spline log_Xe_of_x_Saha_only_Spline{"Xe_Saha_only"};
    Spline log_ne_of_x_spline{"ne"};
    Spline tau_of_x_spline{"tau"}; 
    //Spline dtaudx_of_x_spline{"dtaudx"};
    Spline ddtaudx_of_x_spine{"ddtaudx"};
    Spline g_tilde_of_x_spline{"g"}; 
    Spline dgdx_of_x_spline{"dgdx"};
    Spline ddgdx_of_x_spine{"ddgx"}; 

  public:

    // Construtors
    RecombinationHistory() = delete;
    RecombinationHistory(
        BackgroundCosmology *cosmo, 
        double Yp,
        double h,
        double TCMB
        );

    // Do all the solving
    void solve();
    
    // Print some useful info about the class
    void info() const;

    // Output some data to file
    void output(const std::string filename) const;

    //Making some splines public so I can get them in the perturbations code, since the binary search for values only take splines?
    Spline log_Xe_of_x_spline{"Xe"};
    Spline dtaudx_of_x_spline{"dtaudx"};

    // Get functions that we must implement
    double tau_of_x(double x) const;
    double dtaudx_of_x(double x) const;
    double ddtauddx_of_x(double x) const;
    double g_tilde_of_x(double x) const;
    double dgdx_tilde_of_x(double x) const;
    double ddgddx_tilde_of_x(double x) const;
    double Xe_of_x(double x) const;
    double ne_of_x(double x) const;
    double get_Yp() const;
    double Xe_of_x_Saha_only(double x) const;
};

#endif
