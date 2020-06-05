#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");
  
  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.7;
  double OmegaB      = 0.046;
  double OmegaCDM    = 0.224;
  double OmegaLambda = 0.72995;
  double Neff        = 3.046;
  double TCMB        = 2.725;

  // Recombination parameters
  double Yp          = 0;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaLambda, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("cosmology.txt");

  // Remove when module is completed
  //return 0;

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp, h, TCMB);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  
  // Remove when module is completed
  //return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  //double kvalue = 0.01/Constants.Mpc;
  double kvalue1 = 0.1/Constants.Mpc;
  pert.output(kvalue1, "perturbations_k0.1.txt");

  double kvalue2 = 0.01/Constants.Mpc;
  pert.output(kvalue2, "perturbations_k0.01.txt");

  double kvalue3 = 0.001/Constants.Mpc;
  pert.output(kvalue3, "perturbations_k0.001.txt");

  double kvalue4 = Constants.k_min; //0.0005/Mpc
  pert.output(kvalue4, "perturbations_k0.0005.txt");

  double kvalue5 = Constants.k_max; //0.3/Mpc
  pert.output(kvalue5, "perturbations_k0.3.txt");

  //pert.output_source_function_terms(kvalue1, "sourcetermsk01.txt");
  //pert.output_source_function_terms(kvalue2, "sourcetermsk001.txt");
  //pert.output_source_function_terms(kvalue3, "sourcetermsk0001.txt");

  // Remove when module is completed
  //return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert);
  power.solve();
  power.output("cells.txt");
  //power.output("cells_SW.txt");
  //power.output("cells_ISW.txt");
  //power.output("cells_Doppler.txt");
  //power.output("cells_lastterm.txt");
  //power.output_integrand_theta_ell("integrand.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
