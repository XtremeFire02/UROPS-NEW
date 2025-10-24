#include <cmath>
#include <iostream>
#include <map>

#include "wind.h"
#include "wave.h"
#include "ship.h"

using namespace std;

constinit const double wind::windMag_limU{24};
constinit const double wind::windMagDefaultStep{1};

constinit const int wind::angleDefaultSize{11}; //relative angle is [0, 180], 10 * 18 = 180, numBins = 10, table size = numBins + 1
constinit const int wind::angleDefaultStep{18};

// Reference: [10] Karlsson, J. (2012). Performace Modelling: Bunker Consumption as a Function of Vessel Speed. http://www.diva-portal.org/smash/record.jsf?pid=diva2:549800
// Table 2 of [10]
// A slightly different version is found in https://www.engineersedge.com/fluid_flow/wind_force_scales_15299.html
const d1 wind::bfTable{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
const d1 wind::windBfTable{0, 0.3, 1.6, 3.5, 5.5, 8.0, 10.8, 13.9, 17.2, 20.8, 24.5, 28.5, 32.6};
const d1 wind::swhBfTable{0, 0.001, 0.2, 0.5, 1.0, 2.0, 3.0, 4.0, 5.5, 7.5, 10.0, 12.5, 16.0}; //corrected for 12,11,10,9

const d1 wind::kdwt280TankerConvBowLadenCx{-0.96, -0.93, -0.85, -0.73, -0.62, -0.47, -0.34, -0.17, -0.06, 0.05, 0.14, 0.22, -0.29, 0.40, 0.53, 0.66, 0.75, 0.79, 0.77};
const d1 wind::lngPrismaticIntegratedCx{-1.02, -0.99, -0.93, -0.67, -0.48, -0.30, -0.14, -0.03, 0.05, 0.11, 0.24, 0.41, 0.56, 0.68, 0.79, 0.88, 0.92, 0.91};
const d1 wind::generalCargoCx{-0.60, -0.87, -1.00, -1.99, -0.88, -0.85, -0.65, -0.42, -0.27, -0.09, 0.09, 0.49, 0.84, 1.39, 1.47, 1.34, 0.92, 0.82};
const d1 wind::teu6800ContainerLadenCx{-1.01, -1.11, -1.11, -0.95, -0.74, -0.46, -0.13, 0.28, 0.77, 1.09, 1.18, 1.19, 0.94};

double wind::interpolBf(const d1& windTable, const d1& bfTable, const double& wind)
{
    //because bin interval are not equal, we cannot use ship::idxOfUniformTable

    if (wind <= windTable.front())
        return bfTable.front();

    if (wind >= windTable.back())
        return bfTable.back();

    int idx{static_cast<int>(windTable.size() - 2)};    //note -2

    for (int i = 0; i < static_cast<int>(windTable.size() - 2); i++)
    {
        if (wind > windTable[i] && wind <= windTable[i + 1])
        {
            idx = i;
            break;
        }
    }

    return wave::interpolate(windTable[idx], windTable[idx + 1], bfTable[idx], bfTable[idx + 1], wind);
}

double wind::mergeWindWaveBf(const double& windMag, const double& swh)
{
    //assert(wave::bfTable.size() == wave::windBfTable.size());
    //assert(wave::bfTable.size() == wave::swhBfTable.size());

    const double windBf{wind::interpolBf(wind::windBfTable, wind::bfTable, windMag)};
    const double waveBf{wind::interpolBf(wind::swhBfTable, wind::bfTable, swh)};

    return max(windBf, waveBf);
}

wind::wind(const ship* _shipP, const int& _knotsLb, const int& _knotsSize, const double& _windMag, const double& _windDegRelShip)
{
    init(_shipP, _knotsLb, _knotsSize, _windMag, _windDegRelShip);
}

// Main Reference:  [15] Kim, Y., Steen, S., Kramel, D., Muri, H., & Strømman, A. H. (2023c). Modelling of ship resistance and power consumption for the global fleet: The MariTEAM model. Ocean Engineering, 281, 114758. https://doi.org/10.1016/j.oceaneng.2023.114758
// Convention: Follows from convention, 0 indicates head winds
void wind::init(const ship* _shipP, const int& _knotsLb, const int& _knotsSize, const double& _windMag, const double& _windDegRelShip)
{
    if (_shipP == nullptr)
    {
        cout << "shipPtr is null" << endl;
        return;
    }

    if (_knotsLb <= 0)
    {
        cout << "knots <= 0" << endl;
        return;
    }

    double windMagTmp{0};

    if (_windMag < wave::smallMag)
        windMagTmp = 0;
    else if (_windMag >= wind::windMag_limU)
        windMagTmp = wind::windMag_limU;
    else
        windMagTmp = _windMag;

    windMagSq = windMagTmp * windMagTmp;

    knotsLb = _knotsLb;
    knotsSize = _knotsSize;

    switch (_shipP->shipTyChar)
    {
        case 't': Cx = wind::kdwt280TankerConvBowLadenCx; break; //we do a copy so that there is memory contention when we use omp
        case 'g': Cx = wind::lngPrismaticIntegratedCx; break;
        case 'a': Cx = wind::generalCargoCx; break;
        case 'c': Cx = wind::teu6800ContainerLadenCx; break;
        default:
            // Default to tanker for now
            Cx = kdwt280TankerConvBowLadenCx;
            break;
    }

    Cx_lb = 0; //angleOfAttack = [0, 180], so lb = 0
    Cx_size = static_cast<int>(Cx.size());
    Cx_step = 180.0 / (Cx_size - 1); //step = interval

    //calc here faster for table
    const double half_Axv{0.5 * wind::Axv(_shipP->B, _shipP->Loa, _shipP->shipTyChar)};
    Cda_0_half_Axv = Cda(0) * half_Axv;
    Cda_phi_half_Axv = Cda(_windDegRelShip) * half_Axv;

    cos_windRelShip_2Mag = cos(wave::degToRad(_windDegRelShip)) * 2 * windMagTmp;

    V = wave::knotsToMeterPerSec(_knotsLb);
    Vsq = V * V;
}

// Reference: [19] Kitamura, F., Ueno, M., Fujiwara, T., & Sogihara, N. (2017b). Estimation of above water structural parameters and wind loads on ships. Ships and Offshore Structures, 12(8), 1100–1108. https://doi.org/10.1080/17445302.2017.1316556
//sct: a has 9 numbers, some are unused cos our ships are limited? is there a idx 0?
double wind::Axv(const double& shipB, const double& Loa, const char& shipTyChar)
{
    constexpr array<double, 9> a{0.000e+00, 0.000e+00, 1.792e+01, 2.606e+01, 1.132e+00, 1.018e+00, 0.000e+00, -2.765e-02, 5.127e-01};
    constexpr array<double, 9> b{-3.211e-05, -3.303e-05, 1.140e+00, -2.447e+00, -6.409e-02, -5.437e-02, 8.992e-02, 5.182e-03, 8.065e-02};
    constexpr array<double, 9> c{2.571e-02, 2.094e-02, -7.515e+01, 2.183e+02, 4.221e+00, 4.203e+00, 8.500e+00, 8.259e-01, -6.399e-01};

    constexpr array<int, 9> form{4, 4, 1, 1, 3, 3, 3, 6, 3};

    //static_assert(a.size() == b.size(), "a, b not equal size");
    //static_assert(a.size() == c.size(), "a, c not equal size");
    //static_assert(a.size() == form.size(), "a, form not equal size");

    const map<char, int> shipCharToIdx{
                                       {'t', 1},   // assume full, if ballast set category from 1 to 0
                                       {'g', 5},   // assume full, if ballast set category from 5 to 4
                                       {'b', 3},   // assume full, if ballast set category from 3 to 2
                                       {'a', 8},   // information not available, assume others
                                       {'c', 6},
                                       {'r', 7} };

    if (const auto it = shipCharToIdx.find(shipTyChar); it == shipCharToIdx.end())
    {
        cout << "Error in Avx: unknown ship char" << endl;
        return 0;
    }
    else
    {
        const int idx{it->second};
        const double X{a[idx] * shipB + b[idx] * Loa + c[idx]}; // Equation 2 of [19]

        const int& f{form[idx]};

        if (f == 1)
            return X;
        else if (f == 2)
            return X * Loa;
        else if (f == 3)
            return X * shipB;
        else if (f == 4)
            return X * Loa * Loa;
        else if (f == 5)
            return X * Loa * shipB;
        else if (f == 6)
            return X * shipB * shipB;
        else
        {
            cout << "Error in Axv: Invalid form" << endl;
            return X;
        }
    }
}

/*
// Reference: [19] Kitamura, F., Ueno, M., Fujiwara, T., & Sogihara, N. (2017b). Estimation of above water structural parameters and wind loads on ships. Ships and Offshore Structures, 12(8), 1100–1108. https://doi.org/10.1080/17445302.2017.1316556
double wind::Axv(const double& shipB, const double& Loa, const char& shipTyChar)
{
    // Table 2 of [19]
    int idx{0};

    switch (shipTyChar)
    {
        case 't': idx = 1; break;          // assume full, if ballast set category from 1 to 0
        case 'g': idx = 5; break;          // assume full, if ballast set category from 5 to 4
        case 'b': idx = 3; break;          // assume full, if ballast set category from 3 to 2
        case 'a': idx = 8; break;          // information not available, assume others
        case 'c': idx = 6; break;
        case 'r': idx = 7; break;
        default:
            cout << "Error in Avx: unknown ship type" << endl;
    }

    // Equation 2 of [19]
    const double X{Af_a[idx] * shipB + Af_b[idx] * Loa + Af_c[idx]};
    double value{X};

    switch (Af_form[idx])
    {
        case 1:     value = X; break;
        case 2:     value = X * Loa; break;
        case 3:     value = X * shipB; break;
        case 4:     value = X * Loa * Loa; break;
        case 5:     value = X * Loa * shipB; break;
        case 6:     value = X * shipB * shipB; break;
        default:
            cout << "Error in Axv: Invalid form" << endl;
    }

    return value;
}
*/

// Reference: [13] ITTC 2021: Preparation, Conduct and Analysis of Speed/Power Trials, 7.5-04-01-01
// Section F.3: Datasets of wind resistance coefficients
// Convention: Follows from convention, ie. 0 is head wind
double wind::Cda(const double& angleOfAttack) const
{
    /*
    if (angleOfAttack <= wave::ep) return -CX.front();
    if (angleOfAttack >= 180.0) return -CX.back();

    const double binSize = 180.0 / (CX.size() - 1);
    const int binNum = static_cast<int>(angleOfAttack / binSize);


    if (binNum < 0 || binNum + 1 >= static_cast<int>(CX.size()))
    {
        cout << ("Error: angleOfAttack out of interpolation range") << endl;
        return wave::noValDouble;
    }

    return -wave::interpolate(binNum * binSize, (binNum + 1) * binSize, CX[binNum], CX[binNum + 1], angleOfAttack);
    */

    const auto idx{ship::idxOfUniformTable(angleOfAttack, Cx_lb, Cx_size, Cx_step)};

    if (idx.first == idx.second)
        return -Cx[idx.first];
    else
        return -wave::interpolate(Cx_lb + idx.first * Cx_step, Cx_lb + idx.second * Cx_step, Cx[idx.first], Cx[idx.second], angleOfAttack);
}

// Reference:  [15] Kim, Y., Steen, S., Kramel, D., Muri, H., & Strømman, A. H. (2023c). Modelling of ship resistance and power consumption for the global fleet: The MariTEAM model. Ocean Engineering, 281, 114758. https://doi.org/10.1016/j.oceaneng.2023.114758
// Eq 6 and 7
/*
double wind::getR_woAirDensity()
{
    // Cda:        air drag coefficient measured at wind tunnel
    // C_X:         wind resistance coefficient, note that Cda = -Cx
    // AXV:        area of maximum transverse section exposed to the wind estimated from Kitamura et al. (2017)
    // V_s:         ship’s speed over ground
    // V_WRef:      relative wind speed at reference height (normally 10m above the water surface)
    // V_WTref:     true wind speed at reference height
    // phi_Wref:    relative wind direction at reference height
    // phi_WTref:   true wind direction at reference height
    // phi_WT:      true wind direction at the vertical position of the anemometer(0 means headwind)
    // phi:         ship heading

    if (windMag < wave::smallMag)
        return wave::smallMag;

    const double V_s = V;               // tentative
    const double V_WTref = windMag;   // tentative
    const double windDegRelShip = wave::radtoDeg(windDegRelShipRad);         // phi_WT - phi
    const double phi_WRef = windDegRelShip;
    const double AXV = AVX(shipB, shipL_OA, shiptyChar);

    // Eq 5
    double V_WRef = sqrt(pow(V_WTref, 2) + pow(V_s, 2) + 2 * V_WTref * V_s * cos(windDegRelShipRad));
    // double V_WRef = windMag;

    // Eq 6
    return 0.5 * Cda(phi_WRef, shiptyChar) * AXV * pow(V_WRef, 2) - 0.5 * Cda(0, shiptyChar) * AXV * pow(V_s, 2);
    // cout << 0.5 * wind::airDensity * getCda(phi_WRef, curShip) * AXV * pow(V_WRef, 2) << endl;
    // cout << 0.5 * wind::airDensity * getCda(0, curShip) * AXV * pow(V_s, 2) << endl;
}
*/

double wind::getR_woAirDensity()
{
    // Cda:        air drag coefficient measured at wind tunnel
    // C_X:         wind resistance coefficient, note that Cda = -Cx
    // AXV:        area of maximum transverse section exposed to the wind estimated from Kitamura et al. (2017)
    // V_s:         ship’s speed over ground
    // V_WRef:      relative wind speed at reference height (normally 10 m above the water surface)
    // V_WTref:     true wind speed at reference height
    // phi_Wref:    relative wind direction at reference height
    // phi_WTref:   true wind direction at reference height
    // phi_WT:      true wind direction at the vertical position of the anemometer(0 means headwind)
    // phi:         ship heading

    if (windMagSq < wave::smallMagSq)
        return wave::smallMag;

    const double V_wRef_sq{windMagSq + Vsq + cos_windRelShip_2Mag * V}; // Eq 5

    return Cda_phi_half_Axv * V_wRef_sq - Cda_0_half_Axv * Vsq;          // Eq 6
}

void wind::makeTable()
{
    if (!table.empty())   //cannot call more than once
        return;

    table.resize(knotsSize);

    for (int i = 0; i < knotsSize; i++)
    {
        V = wave::knotsToMeterPerSec(knotsLb + i); //use member obj to loop new val is more concise
        Vsq = V * V;

        table[i] = getR_woAirDensity();
    }

    Cx.clear();
}
