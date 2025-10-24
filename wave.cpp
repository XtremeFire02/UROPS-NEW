//#include <assert.h>
#include <mutex>
#include <atomic>
#include <chrono>

#include <cmath>
#include <fstream>
#include <iostream>

#include "ship.h"

using namespace std;

constinit const int wave::angleDefaultSize{5}; //relative angle is [0, 180], 4 * 45 = 180, numBins = 4, table size = numBins + 1
constinit const int wave::angleDefaultStep{45}; //may be too coarse, maybe use 15

constinit const double wave::swhDefaultStep{0.5};
constinit const double wave::pp1dDefaultStep{1};

constinit const double wave::pp1d_limL{2};
constinit const double wave::pp1d_limU{20};
constinit const double wave::swh_limU{8};

//constinit const int wave::waveNumMaxIter{100};
//constinit const double wave::waveNumTol{1e-5}; //changed to 1e-5

//Estimation and Analysis of JONSWAP Spectrum Parameter Using Observed Data around Korean Coast, gamma = 1.4, 2.1
const d1 wave::jsGammaVec{1.0, 3.3};//must be increasing order, first value must be 1 (for PM)

const d2 wave::R_aw_a{
               d1{-10, 7.55, 0.17},            // tanker = 0
               d1{-9.99, 10.0, 8.16},          // gas = 1
               d1{-1.79, -0.23, 0.34},         // bulker = 2
               d1{ 1.95, 0.92, 0.31},          // cargo = 3
               d1{-9.21, 4.52, 2.00},          // container = 4
               d1{-9.35, -7.91, 6.54}};        // roro = 5

const d2 wave::R_aw_b{
               d1{0.52, 1.08, 1.08},          // tanker = 0
               d1{0.26, 2.0, 1.05},           // gas = 1
               d1{0.59, 1.84, 1.92},          // bulker = 2
               d1{0.11, 0.98, 0.74},          // cargo = 3
               d1{0.55, 0.88, 0.95},          // container = 4
               d1{1.01, 0.71, 0.84}};         // roro = 5

//R_aw_c[shipTy][R_aw_c=0 or 1]
const d2 wave::R_aw_c{
               d1{0.04, 1.41},                  // tanker = 0
               d1{-1.2, 1.6},                   // gas = 1
               d1{0.0, 0.01},                   // bulker = 2
               d1{0.21, 0.36},                  // cargo = 3
               d1{0.04, -0.08},                 // container = 4
               d1{1.49, 0.05}};                 // roro = 5

//should be constexpr but msvc does not support currently
const double wave::ceil_4_pitchGyration{ceil(wave::Rlp_pitchGyration * 4)};
const double wave::floor_4_pitchGyration{floor(wave::Rlp_pitchGyration * 4)};

ofstream s_waveLog;
mutex s_waveLogMtx;
atomic<bool> s_waveLogOn{false};
bool s_waveLogHeader{false};

//https://alexsm.com/cpp-closest-lower-bound/
int wave::nearestIdxOfSortedVec(const std::vector<double>& sortedVec, const double& x)
{
    if (sortedVec.empty())
        return -1;

    const auto iter_geq{lower_bound(sortedVec.begin(), sortedVec.end(), x)};

    if (iter_geq == sortedVec.begin())
        return 0;

    if (iter_geq == sortedVec.end())
        return static_cast<int>(sortedVec.size()) - 1;

    const auto& a{*(iter_geq - 1)};
    const auto& b{*(iter_geq)};

    if (fabs(x - a) < fabs(x - b))
        return static_cast<int>(iter_geq - sortedVec.begin()) - 1;

    return static_cast<int>(iter_geq - sortedVec.begin());
}

double wave::interpolate2d(const double& x, const double& y, const double& g00, const double& g01, const double& g10, const double& g11,
    const double& x0, const double& x1, const double& y0, const double& y1)
{
    //grid00 = f(x0, y0), grid01 = f(x1, y0), grid10 = f(x0, y1), grid11 = f(x1, y1)
    //middle 2 are diff from libInterpolate

    //https://www.adicot.com/bilinear-interpolation
    //double x0 = 75; double y0 = 63;
    //double x1 = 85; double y1 = 67;

    //double g00 = 55.04;
    //double g01 = 52.59;
    //double g10 = 59.00;
    //double g11 = 56.34;

    //cout << wave::interpolate2d(79, 65, g00, g01, g10, g11, x0, x1, y0, y1) << endl; //ans: 55.998

    if (x < x0 || x > x1 || y < y0 || y > y1)
    {
        return -1;
    }

    const double x_ratio{(x - x0) / (x1 - x0)};
    const double y_ratio{(y - y0) / (y1 - y0)};

    //const double v1 = grid[0][0] * (1.0 - x_ratio) * (1.0 - y_ratio);
    //const double v2 = grid[0][1] * x_ratio * (1.0 - y_ratio);
    //const double v3 = grid[1][0] * (1.0 - x_ratio) * y_ratio;
    //const double v4 = grid[1][1] * x_ratio * y_ratio;

    const double v1{g00 * (1.0 - x_ratio) * (1.0 - y_ratio)};
    const double v2{g01 * x_ratio * (1.0 - y_ratio)};
    const double v3{g10 * (1.0 - x_ratio) * y_ratio};
    const double v4{g11 * x_ratio * y_ratio};

    return v1 + v2 + v3 + v4;
}

//https://stackoverflow.com/questions/67894740/much-faster-interpolation-in-python
double wave::interpolate(const double& x1, const double& x2, const double& y1, const double& y2, const double& x)
{
    return ((y2 - y1) * x + x2 * y1 - x1 * y2) / (x2 - x1);
}

double wave::knotsToMeterPerSec(const double& V)
{
    return V * wave::meterPerSecOfOneKnot;
}

double wave::meterPerSecToKnots(const double& V)
{
    return V / wave::meterPerSecOfOneKnot;
}

double wave::celsiusToKelvin(const double& _tempCel)
{
    if (_tempCel < -273.15)
    {
        //cout << "Error: tempCel < 273.15";
        return wave::noValDouble;
    }
    else 
    {
        return max(0.0, _tempCel + 273.15);  //lowest value of Kelvin is 0 (absolute zero)
    }  
}

double wave::waterDensity(const int& tempCel)
{
    //static_assert(wave::waterTempCelTable.size() == wave::waterDensityTable.size(), "waterDensity table size wrong");

    return wave::waterDensityTable[wave::waterTempIdx(tempCel)];
}

double wave::waterAbsViscosity(const int& tempCel)
{
    //static_assert(wave::waterTempCelTable.size() == wave::waterAbsViscosityTable.size(), "waterAbsViscosity table size wrong");

    return wave::waterAbsViscosityTable[wave::waterTempIdx(tempCel)];
}

int wave::waterTempIdx(const int& tempCel)
{
    if (tempCel <= wave::waterTempCelTable.front())
        return 0;
    else
    {
        //the table is spaced 1 deg apart & arg is int, so there is a simple way
        const int idx{tempCel - wave::waterTempCelTable.front()};
        constexpr int lastIdx{static_cast<int>(wave::waterTempCelTable.size()) - 1};

        if (idx >= lastIdx)
            return lastIdx;
        else
            return idx;
    }
}

double wave::degToRad(const double& angle)
{
    constexpr double fac{numbers::pi / 180.0};
    return angle * fac;
}

double wave::radToDeg(const double& angle)
{
    constexpr double fac{180.0 / numbers::pi};
    return angle * fac;
}

double wave::swapAngleConvention(const double& angle)
{
    return wave::make360Deg(fmod(angle + 180, 360));
}

double wave::make360Deg(double angle)
{
    while (angle > 360)
        angle -= 360;

    while (angle < 0)
        angle += 360;

    return angle;
}

double wave::make180Deg(double angle)
{
    //double angle360 = wave::make360Deg(angle);
    if (angle > 180.0)
        return 360.0 - angle;
    else
        return angle;
}

// bearing and shipBearing must be same convention
// If bearing is of to_convention, relative_angle = 0 -> tail waves/tail winds
// If bearing is of from_convention, relative_angle = 0 -> head waves/head winds
double wave::relAngleDegMod180(const double& bearing, const double& shipBearing)
{
    double relAngle{wave::make360Deg(wave::make360Deg(bearing) - wave::make360Deg(shipBearing))};
    // relAngle must be 0 to 180

    // As an example, 181, 179 should be the same, which justifies the below
    // If (relAngle > 180)  relAngle = 360 - relAngle;

    if (relAngle < 0)
        relAngle += 360;

    if (relAngle > 180)
        relAngle = 360 - relAngle;

    return relAngle;
}

double wave::relAngleDegMod180_frConvention(const double& angleDeg_toConvention, const double& shipDeg_toConvention)
{
    const double angleFrom{wave::swapAngleConvention(angleDeg_toConvention)};
    const double shipFrom{wave::swapAngleConvention(shipDeg_toConvention)};

    return wave::relAngleDegMod180(angleFrom, shipFrom);
}

void wave::init(const ship* _shipP, const int& _knotsLb, const int& _knotsSize, double _swh, double _pp1d,
                const double& _waveDegRelShipLb, const int& _waveDegRelShipSize, const double& _windDegRelShipStep)
{
    //we just do warning
    if (_knotsLb <= 0)
    {
        cout << "knotsLb <= 0" << endl;
        return;
    }

    if (_knotsSize <= 0)
    {
        cout << "knotsSize <= 0" << endl;
        return;
    }

    if (_waveDegRelShipLb < -wave::ep)
    {
        cout << "_waveDegRelShipLb < 0" << endl;
        return;
    }

    if (_waveDegRelShipSize <= 0)
    {
        cout << "_waveDegRelShipSize <= 0" << endl;
        return;
    }

    if (_windDegRelShipStep <= 0)
    {
        cout << "_windDegRelShipStep <= 0" << endl;
        return;
    }

    if (_shipP == nullptr)
    {
        cout << "shipPtr is null" << endl;
        return;
    }

    //gammaIdx are clipped before passing in, we just do warning
    //if (gammaIdx < 0 || gammaIdx >= static_cast<int>(wave::jsGammaVec.size()))
        //cout << "gammaIdx out of range" << endl;
    //end warning

    switch (_shipP->shipTyChar) //must come first
    {
        case 't': curShip = tanker; break;
        case 'g': curShip = gas;    break;
        case 'b': curShip = bulker; break;
        case 'a': curShip = cargo;  break;
        case 'c': curShip = container; break;
        case 'r': curShip = roro;   break;
        default:
            cout << "Error: Unknown ship type" << endl;
    }

    shipB = _shipP->B;
    shipCb = _shipP->Cb;
    shipLpp = _shipP->Lpp;

    shipT = _shipP->T;
    shipTmax = _shipP->Tmax;

    //must come before pp1d
    //warning
    if (_swh > wave::swh_limU)
    {
        _swh = wave::swh_limU;
    }
    else if (_swh < wave::ep)
    {
        _swh = wave::ep;
        swhZero = true;
    }

    constexpr double pp1d_boundary{5};

    if (_pp1d < pp1d_boundary)
        _pp1d = 5 * sqrt(_swh); //1 param PM spectrum, Table 2 of Estimation of wave spectra from wave height and period, DJT Carter 1982

    //must come after trying 1 param PM
    //warning
    if (_pp1d > wave::pp1d_limU)
        _pp1d = wave::pp1d_limU;
    else if (_pp1d < wave::pp1d_limL)
        _pp1d = wave::pp1d_limL;

    //----
    knotsLb = _knotsLb;
    knotsSize = _knotsSize;

    waveDegRelShipLb = _waveDegRelShipLb;
    waveDegRelShipSize = _waveDegRelShipSize;
    waveDegRelShipStep = _windDegRelShipStep;

    knots = knotsLb;
    //----

    // Reference: [11] Kim, Y., Esmailian, E., & Steen, S. (2022). A meta-model for added resistance in waves. Ocean Engineering, 266, 112749. https://doi.org/10.1016/j.oceaneng.2022.112749
    // Table 1
    Le = shipCb * Le_regGradTable[curShip] + Le_regInterceptTable[curShip];

    // Reference: [11] Kim, Y., Esmailian, E., & Steen, S. (2022). A meta-model for added resistance in waves. Ocean Engineering, 266, 112749. https://doi.org/10.1016/j.oceaneng.2022.112749
    // Table 2
    Lr = shipCb * Lr_regGradTable[curShip] + Lr_regInterceptTable[curShip];

    pow_Lpp_div_B = pow(shipLpp / shipB, -2.66);
    sqrt_Lpp_div_g = sqrt(shipLpp / wave::gravityConst);

    // Reference: [12] Holthuijsen, L. H. (2007). Waves in oceanic and coastal waters. https://doi.org/10.1017/cbo9780511618536
    // Chapter 4: Statistics, Eq 4.2.20
    // swh = 4 sqrt(m_0) = 4 dst; As in sinuisoidal waves, std = amplitude / sqrt(2), swh = 2sqrt(2) A
    zetaSq = _swh * _swh / 8;

    rho_g_zeta2_shipB = wave::rho_g * zetaSq * shipB;
    rho_g_zeta2_shipB2_div_Lpp = rho_g_zeta2_shipB * shipB / shipLpp;

    Rlp_wBarFac = 2.142 * pow(wave::Rlp_pitchGyration, wave::oneThird) * (1 - 0.111 / shipCb * (log(shipB / shipTmax) - log(2.75))) * pow(shipCb / 0.65, 0.17);

    Rlp_a1Fac = 60.3 * pow(shipCb, 1.34) * pow(4 * wave::Rlp_pitchGyration, 2) / log(shipB / shipTmax);

    //must come after Rlp_a1Fac
    Rlp_a1_f0 = Rlp_a1_head(0, 0);
    Rlp_a1_f1 = Rlp_a1_head(numbers::pi, 0);
    //----

    spectrum_B = pow(wave::twoPi / _pp1d, 4) / 0.8;           //Eq 7 of [14]
    //same result as double B{691.18 / pow(0.7718 * _pp1d, 4)}; //Eq 5 of [14], mean = 0.7718 peak is from Table 1 of Estimation of wave spectra from wave height and period, DJT Carter 1982

    spectrum_A = 0.25 * spectrum_B * _swh * _swh;     // Eq 4 of [14]
    spectrum_b = pow(spectrum_B, -0.25);              // Scale factor in Eq 18 of [14], equivalent with Eq 20, 21 of [15]

    spectrum_m0 = 0.25 * spectrum_A / spectrum_B;              // Eq 3 of [14]

    //sct byTable cannot have
    //this->logInitParamsWithR(_knotsLb, _knotsSize, swh, _pp1d,
    //_waveDegRelShipLb, _waveDegRelShipSize,
    //_windDegRelShipStep, gammaIdx);
}

void wave::initSpectrum(const int& gammaIdx)
{
    calcFreqSpectrum_sct(gammaIdx);
    calcDirSpectrum(waveDegRelShipLb);
}

wave::wave(const ship* _shipP, const int& _knotsLb, const int& _knotsSize, const double& _swh, const double& _pp1d,
    const double& _waveDegRelShipLb, const int& _waveDegRelShipSize, const double& _windDegRelShipStep)
{
    init(_shipP, _knotsLb, _knotsSize, _swh, _pp1d, _waveDegRelShipLb, _waveDegRelShipSize, _windDegRelShipStep);
}

void wave::calcFreqSpectrum_sct(const int& gammaIdx)
{
    sVec.clear();
    wVec.clear();

    if (swhZero)
        return;

    constexpr int nBins{80};

    sVec.reserve(nBins + 1);                   // + 1, cos we are doing interval
    wVec.reserve(nBins + 1);

    const bool usePM = gammaIdx == 0 ? true : false;

    double gamma{1}, lambda{1};

    if (!usePM)
    {
        gamma = wave::jsGammaVec[gammaIdx];

        const double lgG{log10(gamma)};
        lambda = 1 + 0.261 * lgG - 0.0525 * lgG * lgG; //Eq 34 of [14]
    }

    //-----
    constexpr double xUb{4};               // Choice of 4 derived from Figure 4 of [14]
    const double xDelta{xUb / nBins};

    int iStart{2}, iEnd{nBins};

    constexpr double minS{3e-4}; //S is elevationSq, sVec.back() is around 3e-4 -> swh = 0.02

    for (int i = 2; i <= nBins; i++) //start from 2, keep away from 0, x = 0 -> getS_of_PM div by zero
    {
        const double x{i * xDelta};
        const double S = usePM ? wave::getS_of_PM(x, spectrum_m0) : wave::getS_of_JS(x, spectrum_m0, lambda);

        if (S > minS)
        {
            iStart = i;
            break;
        }
    }

    for (int i = nBins; i > iStart; i--)
    {
        const double x{i * xDelta};
        const double S = usePM ? wave::getS_of_PM(x, spectrum_m0) : wave::getS_of_JS(x, spectrum_m0, lambda);

        if (S > minS)
        {
            iEnd = i;
            break;
        }
    }

    if (iEnd - iStart >= 2) //integral needs at least 2
    {
        if (iEnd <= 0.6 * nBins)
        {
            for (int i = iStart; i <= iEnd; i++)
            {
                const double x{i * xDelta};
                const double S = usePM ? wave::getS_of_PM(x, spectrum_m0) : wave::getS_of_JS(x, spectrum_m0, lambda);

                sVec.emplace_back(S);
                wVec.emplace_back(x / spectrum_b);                           // Eq 16, 18, note divide by b
            }
        }
        else
        {
            const int tail{static_cast<int>(floor(0.7 * iEnd))};

            for (int i = iStart; i < tail; i++)
            {
                const double x{i * xDelta};
                const double S = usePM ? wave::getS_of_PM(x, spectrum_m0) : wave::getS_of_JS(x, spectrum_m0, lambda);

                sVec.emplace_back(S);
                wVec.emplace_back(x / spectrum_b);
            }

            for (int i = tail; i <= iEnd; i += 3) //tail changes slowly -> skip some points
            {
                const double x{i * xDelta};
                const double S = usePM ? wave::getS_of_PM(x, spectrum_m0) : wave::getS_of_JS(x, spectrum_m0, lambda);

                sVec.emplace_back(S);
                wVec.emplace_back(x / spectrum_b);
            }

            //cout << iStart << " " << iEnd << " sVec.size():" << sVec.size() << " " << sVec.back() << endl;
        }
    }

    //_sVec may be empty when minS is not reached. this is fine cos getR return smallMag when sVec.size() < 2
}

double wave::getS_of_PM(const double& x, const double& m0)
{
    const double xSq{x * x};
    const double x4{xSq * xSq};

    return (4.0 / (x4 * x) * exp(-1.0 / x4)) * m0;    // Eq 16, 18, note multiply by m0
}

double wave::getS_of_JS(const double& x, const double& m0_PM, const double& lambda)
{
    //const double x_m{ pow(4.0 / 5.0, 0.25) };    // Page 469 of [14] (Non-dimensional frequency where peak sectrum occurs)
    constexpr double x_m{0.9457416090031758};
    const double x_PM{x_m + lambda * (x - x_m)};

    return lambda * wave::getS_of_PM(x_PM, m0_PM); //Eq 35 of [14]
}

// Reference: [13] ITTC 2021: Preparation, Conduct and Analysis of Speed/Power Trials, 7.5-04-01-01
// Equation 13 at Directional Spectrum, Page 24, lower bound 0 to upper bound 360 degrees
/*
void wave::calcDirSpectrum(const double& waveDeg, d1& _gVec, d1& _alphaVec)
{
    const int nBins{40};

    // Reference: [13] ITTC 2021: Preparation, Conduct and Analysis of Speed/Power Trials, 7.5-04-01-01
    // Page 24: s = 1 for wind waves and s = 75 for swells

    constexpr double s{1.0};
    constexpr double two_s{s * 2};
    const double C = pow(2, two_s) / numbers::pi * pow(tgamma(s + 1), 2) / tgamma(2 * s + 1); // C(s) = 2^{2s}/pi * [gamma(s+1)]^2 / gamma(2s+1)

    _gVec.clear();
    _alphaVec.clear();

    _gVec.resize(nBins);
    _alphaVec.resize(nBins);

    const double dAlpha{360.0 / nBins};

    for (int i = 0; i < nBins; i++)
    {
        const double sigma{i * dAlpha};                        // input in degrees
        _alphaVec[i] = wave::degToRad(sigma);                 // save in radians, to_convention

        const double diff = relAngleRad(sigma, waveDeg);       // relative angle between main wave and spectrum wave

        //if (diff <= wave::halfPi)
            //_gVec[i] = C * pow(cos(diff), two_s);
        //else
            //_gVec[i] = 0.0;
    }

    waveDegRelShipFrConvention = waveDeg; //cleaner to put here
}
*/

void wave::calcDirSpectrum(const double& waveDeg)
{
    // Reference: [13] ITTC 2021: Preparation, Conduct and Analysis of Speed/Power Trials, 7.5-04-01-01
    // Page 24: s = 1 for wind waves and s = 75 for swells

    //constexpr double s{1.0};
    //constexpr double two_s{s * 2};
    //const double C = pow(2, two_s) / numbers::pi * pow(tgamma(s + 1), 2) / tgamma(2 * s + 1); // C(s) = 2^{2s}/pi * [gamma(s+1)]^2 / gamma(2s+1)

    constexpr int nBins{40};
    constexpr double C{0.6366197723675814}; //C(s) = 2^{2s}/pi * [gamma(s+1)]^2 / gamma(2s+1), s = 1

    //Eq 13 lower upper integral limits are -90 <= alpha - waveDeg <= 90
    constexpr double U{90}, L{-U};
    constexpr double delta{(U - L) / nBins};

    const int iStart{static_cast<int>(floor(L / delta))}; //floor should be constexpr, msvc yet to support it
    const int iEnd{static_cast<int>(ceil(U / delta))};

    alphaVec.resize(iEnd - iStart + 1);
    gVec.resize(alphaVec.size());

    for (int i = 0; i < static_cast<int>(alphaVec.size()); i++)
    {
        const double alphaLessWave{(iStart + i) * delta};
        const double cosine{cos(wave::degToRad(alphaLessWave))};

        alphaVec[i] = wave::degToRad(alphaLessWave + waveDeg);  //save in radians
        gVec[i] = C * cosine * cosine; //faster than emplace_back
    }

    waveDegRelShipFrConvention = waveDeg; //cleaner to put here

    //cout << waveDeg << " L = " << L << " U = " << U << " delta = " <<  delta << " iStart= " << iStart << " iEnd = " << iEnd << endl;
    //cout << wave::radToDeg(alphaVec.front()) << " " << wave::radToDeg(alphaVec.back()) << endl;
}

// Reference: [17] Jeng, D. (2001b). Wave dispersion equation in a porous seabed. Ocean Engineering, 28(12), 1585–1599. https://doi.org/10.1016/s0029-8018(00)00068-8
// Supplementary Reference: https://en.wikipedia.org/wiki/Dispersion_(water_waves)
// Given by Wave Dispersion Equation in Eq 22 of [17]: ww^2 = gktanh(kh)
double wave::waveNum(const double& ww)
{
    // Initial Value: Deep water approximation given in https://en.wikipedia.org/wiki/Dispersion_(water_waves)
    const double wwSq{ww * ww};
    double k{wwSq / wave::gravityConst};
    // cout << k << endl;

    /*
    for (int iter = 0; iter < wave::waveNumMaxIter; ++iter)
    {
        const double kh{k * wave::waveNumWaterDepth};
        const double tanh_kh{tanh(kh)};
        const double sech_kh{1.0 / cosh(kh)};

        const double f{wave::gravityConst * k * tanh_kh - wwSq};
        const double df{wave::gravityConst * (tanh_kh + kh * sech_kh * sech_kh)};

        // Newton's Step
        const double kNew{k - f / df};

        if (fabs(kNew - k) < wave::waveNumTol)
            return kNew;

        k = kNew;
    }
    */

    return k;
}

// Reference: [13] ITTC 2021: Preparation, Conduct and Analysis of Speed/Power Trials, 7.5-04-01-01
// Trapezoidal Method to calculate integral S(w) dw
double wave::integral1d(const d1& _sVec, const d1& _wVec)
{
    //warning
    if (_wVec.size() < 2 || _sVec.size() != _wVec.size())
    {
        cout << "Error in integral1d()";
        return wave::noValDouble;
    }

    double sum{0};

    for (int i = 0; i < static_cast<int>(_wVec.size()) - 1; i++)  // -1 cos we use i + 1 later
    {
        const double y{0.5 * (_sVec[i] + _sVec[i + 1])};
        const double dx{_wVec[i + 1] - _wVec[i]};                   // vec must be increasing

        sum += (dx * y);
    }

    return sum;
}

// Reference: [13] ITTC 2021: Preparation, Conduct and Analysis of Speed/Power Trials, 7.5-04-01-01
// Trapezoidal Method to calculate integral integral  S(w) G(alpha) R(w, alpha) dw dAlpha
/*
//sct alphaVec: angle between shiphead and component waves; 0 means head waves
double wave::integral2d()
{
    //warning
    if (wVec.size() < 2 || sVec.size() != wVec.size())
    {
        cout << "Error in integral2d()";
        return wave::noValDouble;
    }

    if (alphaVec.size() < 2 || gVec.size() != alphaVec.size())
    {
        cout << "Error in integral2d()";
        return wave::noValDouble;
    }

    d1 outerSum(gVec.size(), 0);

    for (int i = 0; i < static_cast<int>(alphaVec.size()); i++)
    {
        d1 sVal{sVec};

        for (auto& s : sVal)
            s *= gVec[i];

        double aa = alphaVec[i];

        if (aa > numbers::pi)
            aa = wave::twoPi - aa;
        else if (aa < 0)
            aa *= -1;

        const double angleRad = aa;

        //const double angleRad = wave::degToRad(wave::relAngleDegMod180(wave::radToDeg(alphaVec[i]), 0.0)); //i think not correct

        // const double beta_from = wave::relAngleRad(wave::radToDeg(alphaVec[i]), waveDegRelShipFrConvention); // 0 = head
        // const double angleRad  = numbers::pi - beta_from; // α (to-conv): π = head, 0 = following

        for (int j = 0; j < static_cast<int>(wVec.size()); j++)
        {
            sVal[j] *= R_aw(angleRad, wVec[j], LP_n_CTH);
        }

        outerSum[i] = wave::integral1d(sVal, wVec);
    }

    return wave::integral1d(outerSum, alphaVec) * 2 / zetaSq;
}
*/

double wave::froudeNum(const double& U, const double& Lpp)
{
    return (U / sqrt(Lpp * wave::gravityConst));
}

// Combined Wave Model
// Main Reference: [11] Kim, Y., Esmailian, E., & Steen, S. (2022). A meta-model for added resistance in waves. Ocean Engineering, 266, 112749. https://doi.org/10.1016/j.oceaneng.2022.112749
// Literature uses to_convention : 180 indicates head waves
double wave::g(const double& aa, const double& bb, const double& vv)
{
    return 0.5 * (1 + tanh(aa * (bb - vv)));
}

double wave::R_aw_combined(const double& aa, const double& ww) const
{
    constexpr double radOf45Deg{numbers::pi * 0.25}, radOf135Deg{numbers::pi * 0.75};

    const double lambda{wave::Rlp_lambda(ww)};
    const double lambdaDivL{lambda / shipLpp};

    const double Rcth{wave::R_aw(aa, ww, CTH)};
    const double Rlp{wave::R_aw(aa, ww, LP)};

    if (aa >= wave::halfPi)
    {
        const double gHead{wave::g(R_aw_a[static_cast<int>(curShip)][static_cast<int>(head)], R_aw_b[static_cast<int>(curShip)][static_cast<int>(head)], lambdaDivL)};
        const double gBeam{wave::g(R_aw_a[static_cast<int>(curShip)][static_cast<int>(beam)], R_aw_b[static_cast<int>(curShip)][static_cast<int>(beam)], lambdaDivL)};

        const double rHead{(1 - gHead) * Rcth + gHead * Rlp};
        const double rBeam{(1 - gBeam) * Rcth + gBeam * Rlp};

        //const double f = 1 - wave::g(R_aw_c[static_cast<int>(curShip)][0], radOf135Deg, aa);
        //return (1 - f) * rHead + f * rBeam;

        const double f{wave::g(R_aw_c[static_cast<int>(curShip)][0], radOf135Deg, aa)};
        return f * rHead + (1 - f) * rBeam;
    }
    else
    {
        const double gBeam{wave::g(wave::R_aw_a[static_cast<int>(curShip)][static_cast<int>(beam)], wave::R_aw_b[static_cast<int>(curShip)][static_cast<int>(beam)], lambdaDivL)};
        const double gFollow{wave::g(R_aw_a[static_cast<int>(curShip)][static_cast<int>(follow)], R_aw_b[static_cast<int>(curShip)][static_cast<int>(follow)], lambdaDivL)};

        const double rBeam{(1 - gBeam) * Rcth + gBeam * Rlp};
        const double rFollow{(1 - gFollow) * Rcth + gFollow * Rlp};

        //const double f = 1 - wave::g(R_aw_c[static_cast<int>(curShip)][1], radOf45Deg, aa);
        //return (1 - f) * rBeam + f * rFollow;

        const double f{wave::g(R_aw_c[static_cast<int>(curShip)][1], radOf45Deg, aa)};
        return f * rBeam + (1 - f) * rFollow;
    }
}

double wave::R_aw(const double& aa, const double& ww, const R_aw_choice& R_aw_method) const
{
    // If frequency is zero, no time-varied waves -> resistance = 0
    if (fabs(ww) < wave::ep)
        return 0;

    // Note that aa is of to convention, ie. 180 -> head waves
    switch (R_aw_method)
    {
        case LP_n_CTH:
            return wave::R_aw_combined(aa, ww);
            break;

        case CTH:
            return wave::Rcth_aw(aa, ww);
            break;

        case LP:
            return wave::Rlp_aw(aa, ww);
            break;

        default:
            cout << "Error: invalid argument for R_aw_method" << endl;
            return wave::noValDouble;
            break;
    }
}

double wave::getR()
{
    const bool WAVE_LOG_ALL = true;

    if (gVec.size() < 2 || sVec.size() < 2) {
        wave::logRow(this, "getR size < 2", 0.0, 0.0, wave::smallMag);
        return wave::smallMag;
    }

    if (WAVE_LOG_ALL) {
        const double R_TOTAL = integral2d_with(R_aw_choice::LP_n_CTH);
        wave::logRow(this, "getR_TOTAL_ONLY", 0.0, 0.0, R_TOTAL);
        return R_TOTAL;
    }
    else {
        const double R_LP = integral2d_with(R_aw_choice::LP);
        const double R_CTH = integral2d_with(R_aw_choice::CTH);
        const double R_TOTAL = integral2d_with(R_aw_choice::LP_n_CTH);
        wave::logRow(this, "getR", R_LP, R_CTH, R_TOTAL);
        return R_TOTAL;
    }
}

void wave::makeTable()
{
    if (!table.empty())   //cannot call more than once
        return;

    wave::makeVec3d(table, wave::jsGammaVec.size(), waveDegRelShipSize, knotsSize, 0.0);

    b1 skipWaveDeg(waveDegRelShipSize, false); //JS not used frequently, skip some waveDeg for speed

    //first and last point are never skipped
    //waveDegRelShipStep may be coarse -> skip every other point but not too many points, hence j += 2
    for (int j = 1; j < waveDegRelShipSize - 1; j += 2)
        skipWaveDeg[j] = true;

    for (int i = 0; i < static_cast<int>(wave::jsGammaVec.size()); i++) //better be first dim
    {
        calcFreqSpectrum_sct(i);

        for (int j = 0; j < waveDegRelShipSize; j++)
        {
            if (i > 0 && skipWaveDeg[j]) //i == 0 is PM spectrum, we do not skip
                continue;

            //cout << "i = " << i << " j= " << j << endl;

            calcDirSpectrum(waveDegRelShipLb + j * waveDegRelShipStep);

            for (int k = 0; k < knotsSize; k++)
            {
                knots = knotsLb + k;    //use member obj to loop new val is more concise
                table[i][j][k] = getR();
            }
        }
    }

    //fill skipped values by taking ave of left and right
    for (int i = 1; i < static_cast<int>(wave::jsGammaVec.size()); i++) //start from 1, cos we do not skip PM spectrum
        for (int j = 1; j < waveDegRelShipSize - 1; j += 2) //faster than using skipWaveDeg[j]
            for (int k = 0; k < knotsSize; k++)
                table[i][j][k] = 0.5 * (table[i][j - 1][k] + table[i][j + 1][k]);

    alphaVec.clear();
    gVec.clear();

    sVec.clear();
    wVec.clear();
}

double wave::integral2d_with(R_aw_choice choice) const
{
    //warning

    if (wVec.size() < 2 || sVec.size() != wVec.size())
    {
        cout << "Error in integral2d_with()";
        return wave::noValDouble;
    }

    //alphaVec & gVec are always >= 2, so not needed
    //if (alphaVec.size() < 2 || gVec.size() != alphaVec.size())
    //{
        //cout << "Error in integral2d_with()";
        //return wave::noValDouble;
    //}

    d1 outer(gVec.size(), 0.0);

    for (int i = 0; i < static_cast<int>(alphaVec.size()); ++i)
    {
        d1 col{sVec};                       // base spectrum across frequency

        for (auto& s : col)
            s *= gVec[i];

        // Change 1
        // compute the relative from convention angle for this bin
        //const double beta_from = wave::relAngleDegMod180_frConvention(wave::radToDeg(alphaVec[i]), waveDegRelShipFrConvention);
        // convert it to the to convention angle
        //const double angleRad = std::numbers::pi - wave::degToRad(beta_from);

        double aa = alphaVec[i];

        if (aa > numbers::pi)
            aa = wave::twoPi - aa;
        else if (aa < 0)    //R_aw needs [0, 180]. by symmetry, we use magnitude when angle is neg
            aa *= -1;

        //const double angleRad = aa;

        for (int j = 0; j < static_cast<int>(wVec.size()); j++)
            col[j] *= R_aw(aa, wVec[j], choice);

        outer[i] = wave::integral1d(col, wVec);
    }

    return wave::integral1d(outer, alphaVec) * 2.0 / zetaSq;  // same norm you use
}

void wave::enableLogger(const string& filepath, bool truncate)
{
    lock_guard<mutex> lk(s_waveLogMtx);

    if (s_waveLog.is_open())
        s_waveLog.close();

    s_waveLog.open(filepath, (truncate ? (ios::out | ios::trunc) : (ios::out | ios::app)));

    if (!s_waveLog.is_open())
    {
        cerr << "Error opening wave log file: " << filepath << "\n";
        return;
    }

    if (truncate || !s_waveLogHeader)
    {
        s_waveLog <<
            "timestamp_unix,tag,"
            "shipType,shipLpp,shipB,shipT,shipCb,"
            "knots,waveDegRelShipFrConvention"
            "R_LP,R_CTH,R_TOTAL\n";
        s_waveLogHeader = true;
    }

    s_waveLog.flush();
    s_waveLogOn.store(true);
}

void wave::disableLogger()
{
    lock_guard<mutex> lk(s_waveLogMtx);
    s_waveLogOn.store(false);

    if (s_waveLog.is_open())
        s_waveLog.close();
}

void wave::startInitLog(const string& filepath)
{
    // Just enable the new logger and truncate file once.
    wave::enableLogger(filepath, /*truncate=*/true);
}

void wave::logRow(const wave* self, const char* tag, double R_LP, double R_CTH, double R_TOTAL)
{
    if (!s_waveLogOn.load()) return;

    using namespace chrono;
    const auto ts = duration_cast<seconds>(system_clock::now().time_since_epoch()).count();

    lock_guard<mutex> lk(s_waveLogMtx);

    s_waveLog
        << ts << ',' << tag << ','
        << static_cast<int>(self->curShip) << ',' << self->shipLpp << ',' << self->shipB << ',' << self->shipT << ',' << self->shipCb << ','
        << self->knots << ',' << self->waveDegRelShipFrConvention << ',' << R_LP << ',' << R_CTH << ',' << R_TOTAL << '\n';

    s_waveLog.flush();
}

void wave::logInitParamsWithR(const int& knotsLb, const int& knotsSize,
    const double& swh, const double& pp1d,
    const double& waveDegRelShipLb, const int& waveDegRelShipSize,
    const double& windDegRelShipStep, const int& gammaIdx,
    const string& filepath) const
{
    ofstream ofs(filepath, ios::out | ios::app);

    if (!ofs.is_open())
    {
        cerr << "Error opening log file: " << filepath << '\n';
        return;
    }

    const double R_LP = integral2d_with(R_aw_choice::LP);
    const double R_CTH = integral2d_with(R_aw_choice::CTH);
    const double R_TOTAL = integral2d_with(R_aw_choice::LP_n_CTH);

    ofs << knotsLb << ','
        << knotsSize << ','
        << swh << ','
        << pp1d << ','
        << waveDegRelShipLb << ','
        << waveDegRelShipSize << ','
        << windDegRelShipStep << ','
        << gammaIdx << ','
        << R_LP << ','
        << R_CTH << ','
        << R_TOTAL << '\n';
}

/*
// Main Reference: [14] Pawlowski, M. (2011). Sea Spectra Revisited. In Fluid Mechanics and its Applications (pp. 573–587). https://doi.org/10.1007/978-94-007-1482-3_32
// Supplementary Reference: [15] Kim, Y., Steen, _sVec., Kramel, D., Muri, H., & Strømman, A. H. (2023c). Modelling of ship resistance and power consumption for the global fleet: The MariTEAM model. Ocean Engineering, 281, 114758. https://doi.org/10.1016/j.oceaneng.2023.114758
void wave::calcFreqSpectrum_PM(const int& nBins, d1& _sVec, d1& _wVec, const double& A, const double& B, const double& b)
{
    const double m0{0.25 * A / B};              // Eq 3 of [14]

    constexpr double xUb{4};               // Choice of 4 derived from Figure 4 of [14]
    const double xDelta{xUb / nBins};

    // Uses Nondimensional ITTC spectrum (Section 4.1 of [14]), lower bound 0 to upper bound 4
    for (int i = 0; i <= nBins; i++)            // Note <= so 1 more
    {
        const double x{i * xDelta};

        if (x < wave::ep)
        {
            _sVec.emplace_back(0);
            _wVec.emplace_back(0);
        }
        else
        {
            //this is a for loop so we optimise for speed
            const double xSq{ x * x };
            const double x4{ xSq * xSq };

            _sVec.emplace_back((4.0 / (x4 * x) * exp(-1.0 / x4)) * m0);   // Eq 16, 18, note multiply by m0
            _wVec.emplace_back(x / b);                           // Eq 16, 18, note divide by b
        }
    }
}

// Reference: [14] Pawlowski, M. (2011). Sea Spectra Revisited. In Fluid Mechanics and its Applications (pp. 573–587). https://doi.org/10.1007/978-94-007-1482-3_32
void wave::calcFreqSpectrum_JS(const int& nBins, d1& _sVec, d1& _wVec, const double& jsGamma, const double& b)
{
    // Used in calculation of JS spectrum
    const double sigmaA{0.07};              // Eq 10 of [14]
    const double sigmaB{0.09};              // Eq 10 of [14]

    const double x_m{pow(4.0 / 5.0, 0.25)};    // Page 469 of [14] (Non-dimensional frequency where peak sectrum occurs)

    const double m0{0.1366 + 0.8755 * jsGamma - 0.0075 * jsGamma * jsGamma};     // Eq 3 of [14]

    constexpr double xUb{4};               // Choice of 4 derived from Figure 4 of [14]
    const double xDelta{xUb / nBins};

    // Uses Nondimensional ITTC spectrum (Section 4.1 of [14]), lower bound 0 to upper bound 4
    for (int i = 0; i <= nBins; i++)             // Note <= so 1 more
    {
        const double x{i * xDelta};

        if (x < wave::ep)
        {
            _sVec.emplace_back(0);
            _wVec.emplace_back(0);
            continue;
        }
        else
        {
            //this is a for loop so we optimise for speed
            const double xSq{ x * x };
            const double x4{ xSq * xSq };

            // Generate PM Spectrum (Refer to implementation of calcFreqSpectrum_PM)
            const double s_pm = 4.0 / (x4 * x) * exp(-1.0 / x4);

            // Eq 19 of [14]: Determine JS peak enhancement (dimensionless)
            const double sigma = (x <= x_m) ? sigmaA : sigmaB;
            const double expo = - (x - x_m) * (x - x_m) / (2.0 * sigma * sigma * x_m * x_m);
            const double s_G = s_pm * pow(jsGamma, exp(expo));

            _sVec.emplace_back(s_G * m0);
            _wVec.emplace_back(x / b);
        }
    }
}

void wave::calcFreqSpectrum(const int& gammaIdx)
{
    sVec.clear();
    wVec.clear();

    constexpr int nBins{80};

    sVec.reserve(nBins + 1);                   // + 1, cos we are doing interval
    wVec.reserve(nBins + 1);

    if (gammaIdx == 0)
        wave::calcFreqSpectrum_PM(nBins, sVec, wVec, spectrum_A, spectrum_B, spectrum_b);
    else
        wave::calcFreqSpectrum_JS(nBins, sVec, wVec, wave::jsGammaVec[gammaIdx], spectrum_b);
}
*/
