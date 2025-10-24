#include <cmath>
#include <iostream>

#include "prop.h"
#include "wave.h"

using namespace std;

//these are used for byTable
constinit const double prop::Rfrac{0.1}; //we iterate prop over range of resistance: calmR * (1 + i), i = 0...size - 1
constinit const int prop::Rsize{31}; //21 - 1 = 30 bins of 0.1 = factor of 3
constinit const int prop::assumedTemp{16};

prop::prop(ship* shipP, const int& knots, const int& tempCel)
{
    if (shipP == nullptr)
    {
        std::cout << "shipPtr is null" << endl;
        return;
    }

    hullEff = (1 - shipP->prop_t) / (1 - shipP->prop_w);   // Eq 14 of [15]

    //The ratio of PE to PT is called the hull efficiency and for most ships is a little greater than unity. This is because the propeller gains from the energy already imparted to the water by the hull

    if (hullEff > 1)
        hullEff = 1;

    const double Va_Dprop{wave::knotsToMeterPerSec(knots) * (1 - shipP->prop_w) * shipP->prop_diameter};   // Eq J-17 of ITTC 2021 [13] * Dprop

    Cth_fac = 8 / (numbers::pi * (1 - shipP->prop_t) * wave::waterDensity(tempCel) * Va_Dprop * Va_Dprop);    // Eq 13 of [15]
}

// Reference: [4] Kristensen, Lutzen. (2013). Prediction of Resistance and Propulsion Power of Ships
// Page 12, Propulsive Efficiencies
double prop::diameter(const char& shipTyChar, const double& shipTmax)
{
    switch (shipTyChar)
    {
        case 't': return 0.395 * shipTmax + 1.30;    break;
        case 'g': return 0.395 * shipTmax + 1.30;    break;          // approximation, value not provided, using tanker value
        case 'b':  return 0.395 * shipTmax + 1.30;   break;
        case 'a':  return 0.623 * shipTmax - 0.16;   break;          // approximation, value not provided, using container value
        case 'c':  return 0.623 * shipTmax - 0.16;   break;
        case 'r': return 0.713 * shipTmax - 0.08;    break;
        default:
            std::cout << "Error: Unknown ship type '" << endl;
            return -1;
    }
}

double prop::Bseries_gFn(const double& R) const
{
    const double Cth{R * Cth_fac};

    const double idealEff{2 / (1 + sqrt(Cth + 1))};

    const double g = Cth >= 7 ? 0.85 : 0.59 + 0.177 * Cth - 0.0462 * pow(Cth, 2) + 0.00518 * pow(Cth, 3) - 0.000205 * pow(Cth, 4); //Pg 15 of [4]

    return idealEff * g * prop::rotationEff * prop::shaftEff;
}

// Reference: [15] Kim, Y., Steen, S., Kramel, D., Muri, H., & Strømman, A.H. (2023c).Modelling of ship resistance and power consumption for the global fleet : The MariTEAM model.Ocean Engineering, 281, 114758. https ://doi.org/10.1016/j.oceaneng.2023.114758
// Equations 11 - 14
double prop::Bseries(const double& R) const
{
    const double Cth{R * Cth_fac};

    //double fCth = min(0.65, 0.81 - 0.014 * C_th);  // Eq 12 of [15], sct : should be max by pg 14 of [6]

    const double openWaterEff{2 / (1 + sqrt(Cth + 1)) * max(0.65, 0.81 - 0.014 * Cth)};    // Eq 11 of [15], Eq 12 of [15]
    return hullEff * openWaterEff * prop::rotationEff * prop::shaftEff; // Eq 10 of [15]
}

// Reference: [4] Kristensen, Lutzen. (2013). Prediction of Resistance and Propulsion Power of Ships
// Page 12, Propulsive Efficiencies
double prop::wakeFrac(ship* shipP, const double& Cb, const double& ldr, const double& Dprop)
{
    const double a{0.1 * shipP->B / shipP->Loa + 0.149};
    const double b{0.05 * shipP->B / shipP->Loa + 0.449};
    const double c{585 - (5027 * shipP->B / shipP->Loa) + 11700 * pow(shipP->B / shipP->Loa, 2)};
    const double w_1{a + b / (c * pow(0.98 - Cb, 3) + 1)};
    const double w_2{(0.025 * static_cast<double>(shipP->hullForm)) / (100 * pow(Cb - 0.7, 2) + 1)};
    const double w_3{max(- 0.18 + 0.00756 / ((Dprop / shipP->Loa) + 0.002), 0.1)};

    const double w = (shipP->shipTyChar == 't' || shipP->shipTyChar == 'g' || shipP->shipTyChar == 'b') ?
        prop::makeFrac(0.7 * (w_1 + w_2 + w_3) - 0.45 + 0.08 * ldr) :
        prop::makeFrac(w_1 + w_2 + w_3);

    if (shipP->shipTyChar == 't' || shipP->shipTyChar == 'b')
        return 0.7 * w - 0.45 + 0.08 * ldr; //Appendix G of [4], ldr is same as M
    else
        return w;
}

// Reference: [4] Kristensen, Lutzen. (2013). Prediction of Resistance and Propulsion Power of Ships
// Page 12, Propulsive Efficiencies
double prop::thrustDeductionFac(ship* shipP, const double& Cb, const double& ldr, const double& Dprop)
{
    const double d{0.625 * shipP->B / shipP->Loa + 0.08};
    const double e{0.165 - (0.25 * shipP->B / shipP->Loa)};
    const double f{525 - (8060 * shipP->B / shipP->Loa) + 20300 * pow(shipP->B / shipP->Loa, 2)};
    const double t_1{d + e / (f * pow(0.98 - Cb, 3) + 1)};
    const double t_2{-0.01 * static_cast<double>(shipP->hullForm)};
    const double t_3{2 * (Dprop / shipP->Loa - 0.04)};

    const double t = (shipP->shipTyChar == 't' || shipP->shipTyChar == 'g' || shipP->shipTyChar == 'b') ?
        prop::makeFrac(t_1 + t_2 + t_3 - 0.26 + 0.04 * ldr) :
        prop::makeFrac(t_1 + t_2 + t_3);

    if (shipP->shipTyChar == 't' || shipP->shipTyChar == 'b')
        return t - 0.26 + 0.04 * ldr; //Appendix G of [4], ldr is same as M
    else
        return t;
}

// If have prop_expandedAreaRatio, prop_pitchRatio, prop_bladeNum, prop_rpm
// Use Oosterveld, M.W.C.; van Oossanen, P. Further Calcr-Analysed Data of the Wageningen B-screw Series
// Regression implementation to be done

//We can use https://github.com/mkergoat/bseries
//and find the max value, using Handymax Paper, A_E / A_O = 0.53, pitch ratio P/D = 0.74 on online website gives about 0.63 max, the below fn gives 0.59
//Surprising Behaviour of the Wageningen B-Screw Series Polynomials, in 2020 highlights some problems, so we ditch this method

double prop::getEff(const double& R) const
{
    const double eff_B{Bseries(R)};
    //const double eff_B_gFn{Bseries_gFn(R)};

    //std::cout << "propEff_B " << eff_B << " propEff_B_g " << eff_B_gFn << endl;

    if (isFrac(eff_B)) //priority
    {
        //std::cout << "propEff_B not frac " << propEff_B << endl;
        return eff_B;
    }/*
    else if (isFrac(eff_B_gFn))
    {
        //std::cout << "propEff_B_g not frac " << propEff_B_g << endl;
        return eff_B_gFn;
    }*/
    else
    {
        //std::cout << "propEff no value " << propEff_B << " " << propEff_B_g << endl;
        return 0.6;
    }
}

bool prop::isFrac(const double& v)
{
    return v > wave::ep && v < 1;
}

// To prevent div by 0, force (0, 1) excluding boundary
double prop::makeFrac(const double& v)
{
    if (v > 1)
        return 1 - wave::ep;
    else if (v < wave::ep)
        return wave::ep;
    else
        return v;
}
