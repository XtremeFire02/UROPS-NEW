#include <cmath>
//#include <iostream>

#include "ship.h"
#include "wave.h"

using namespace std;

double ship::bioFoul_Cf(const double& V, const double& Cf, const double& Rn, const int& tempCel) const
{
    //const double deltaCf_Granville = bioFoul_GranvilleCorrection(V, Cf, Rn, tempCel); //needs iter hence slow

    // cout << "Townsin" << deltaCf_Townsin << endl;
    // cout << "Granville" << deltaCf_Granville << endl;
    return bioFoul_townsin(Rn);
}

// Reference: [6] Hakim, M., Suastika, I., & Utama, I. (2023). A practical empirical formula for the calculation of ship added friction-resistance due to (bio)fouling. Ocean Engineering, 271, 113744. 
double ship::bioFoul_townsin(const double& Rn) const
{
    const double ks_m{ks * 1e-6};   // Convert from micrometers to meters
    return max((44 * (pow(ks_m / Lpp, wave::oneThird) - 10 * pow(Rn, -wave::oneThird)) + 0.125) * 1e-3, 0.0);    // Eq 1
}

/*
// Needs more details to work out the actual algorithm, current instructions too unclear
// Reference: [6] Hakim, M., Suastika, I., & Utama, I. (2023). A practical empirical formula for the calculation of ship added friction-resistance due to (bio)fouling. Ocean Engineering, 271, 113744.
double ship::bioFoul_granville(const double& V, const double& Cf, const double& Rn, const int& tempCel) const
{
    const double ks_m{ks * 1e-6};   // Convert from micrometers to meters
    constexpr double kappa{0.40};              // von Kármán constant
    constexpr int maxIter{100};

    //const double eta{wave::waterAbsViscosity(tempCel)};     // Absolute Viscosity
    const double nu{wave::waterAbsViscosity(tempCel) / wave::waterDensity(tempCel)};     // Kinematic Viscosity

    // Reference: [7] Schultz, M. P., & Myers, A. (2003). Comparison of three roughness function determination methods. Experiments in Fluids, 35(4), 372–379. https://doi.org/10.1007/s00348-003-0686-x
    // Obtain roughness function
    auto f = [&](double k) {return 1 / kappa * log(1.0 + k); };  // Eq 11 of [7]
    
    const double _U_tau = V * sqrt(Cf * 0.5);       // Eq 8 of [6]
    const double _ksPlus = ks_m * _U_tau / nu;      // Eq 7 of [6] 

    const double Cfs{Cf};
    double deltaUplus{0}, Cfr{0}, ksPlus{_ksPlus};

    for (int iter = 0; iter < maxIter; iter++)
    {
        // Calc deltaUprime
        double deltaUprime = f(ksPlus);             // Eq 7 of [6] 

        // Calc Cfr                             // Eq 6 of [6]
        const double correction = deltaUprime * kappa / log(10.0);
        const double denom = log10(Rn - correction) - 1.729;
        Cfr = 0.0795 / (denom * denom);

        // Calc U_tau
        const double U_tauNew = V * sqrt(Cfr* 0.5);
        ksPlus = ks_m * U_tauNew / nu;
        deltaUprime = f(ksPlus);

        // Calc deltaUplus                       // Eq 2 of [6] 
        const double sqrt_Cfs = sqrt(Cfs * 0.5);
        const double sqrt_Cfr = sqrt(Cfr * 0.5);

        const double deltaUplusNew = (1.0 / sqrt_Cfs) - (1.0 / sqrt_Cfr)
            - 19.7 * (sqrt_Cfs - sqrt_Cfr)
            - (1.0 / kappa) * deltaUprime * sqrt_Cfr;

        // Calc ksPlus                           // Eq 3 of [6] 
        double bracket = 1.0 - (1.0 / kappa) * sqrt_Cfr + (1.0 / kappa) * (1.5 / kappa - deltaUprime) * Cfr * 0.5;
        double ksPlusNew = (ks_m / Lpp) * (Rn * Cfr * 0.5) * 1.0 / sqrt_Cfr * bracket;

        // Check convergence                    
        if (abs(deltaUplusNew - deltaUplus) < wave::ep && abs(ksPlusNew - ksPlus) < wave::ep)
        {
            deltaUplus = deltaUplusNew;
            ksPlus = ksPlusNew;
            break;
        }

        deltaUplus = deltaUplusNew;
        ksPlus = ksPlusNew;
    }
    // cout << Cfr - Cfs << endl;
    return Cfr - Cfs; //will this be neg
}
*/
