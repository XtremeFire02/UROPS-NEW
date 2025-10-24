#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>

#include "prop.h"
#include "ship.h"
#include "wave.h"
#include "wind.h"
#include "rotorF.h"

using namespace std;

constinit const int ship::minKnotsAllShip{4};
constinit const int ship::maxKnotsAllShip{24};
constinit const int ship::maxRotorF{4};
constinit const int ship::maxSails{4};

constinit const double ship::Loa_limL{50};
constinit const double ship::Tmax_limL{4};
constinit const double ship::B_limL{10};
constinit const double ship::displaceVol_limL{5000};
constinit const double ship::sfoc_limL{120};
constinit const double ship::Abt_limL{4};

constinit const double ship::Loa_limU{410};              //suezmax is 400
constinit const double ship::Tmax_limU{22};              //suezmax is 21
constinit const double ship::B_limU{80};                 //suezmax is 78
constinit const double ship::displaceVol_limU{200000};   //suezmax is 160000
constinit const double ship::sfoc_limU{250};
constinit const double ship::Abt_limU{80};

constinit const double ship::ks_limU{50000};             //50000 -> 5cm
constinit const double ship::Cs_limU{4};

constinit const double ship::propEff_limL{0.2};
constinit const double ship::propEff_limU{0.8};

constinit const double ship::waveUpRatioCap{1};
constinit const double ship::waveDownRatioCap{1};
constinit const double ship::windUpRatioCap{1};
constinit const double ship::windDownRatioCap{0.33};

constinit const double ship::windAid_windMagStep{1.0};
constinit const int ship::windAid_windMagSize{30};

#ifdef SR_MODEL
#include "fmt/format.h"

void ship::printJonswapGamma() //sct
{
    for (const auto& it : wave::jsGammaVec)
        cout << it << endl;
}

void ship::printLimit()
{
    fmt::println("");
    fmt::println("---");
    fmt::println("");

    fmt::println("# {}", "Limits");

    fmt::println("{:<30}{:<15.1f}{:<25}{:.1f}   ", "Loa_lb: ", ship::Loa_limL, "Loa_ub: ", ship::Loa_limU);
    fmt::println("{:<30}{:<15.1f}{:<25}{:.1f}   ", "Tmax_lb: ", ship::Tmax_limL, "Tmax_ub: ", ship::Tmax_limU);
    fmt::println("{:<30}{:<15.1f}{:<25}{:.1f}   ", "B_lb: ", ship::B_limL, "B_ub: ", ship::B_limU);
    fmt::println("{:<30}{:<15.1f}{:<25}{:.1f}   ", "displaceVol_lb: ", ship::displaceVol_limL, "displaceVol_ub: ", ship::displaceVol_limU);
    fmt::println("{:<30}{:<15.1f}{:<25}{:.1f}   ", "sfoc_lb: ", ship::sfoc_limL, "sfoc_ub: ", ship::sfoc_limU);
    fmt::println("{:<30}{:<15.1f}{:<25}{:.1f}   ", "sfoc_lb: ", ship::Abt_limL, "sfoc_ub: ", ship::Abt_limU);

    fmt::println("{:<30}{:<15.1f}", "ks_ub: ", ship::ks_limU);
    fmt::println("{:<30}{:<15.1f}", "Cs_ub: ", ship::Cs_limU);
    fmt::println("{:<30}{:<15d}", "maxRotorF: ", ship::maxRotorF);
    fmt::println("{:<30}{:<15d}", "maxSails: ", ship::maxSails);

    fmt::println("{:<30}{:<15d}{:25}{:d}   ", "minKnotsAllShip: ", ship::minKnotsAllShip, "maxKnotsAllShip: ", ship::maxKnotsAllShip);

    fmt::println("");
}
#endif

/*
// Reference: https://www.afs.enea.it/project/neptunius/docs/fluent/html/ug/node294.htm
double ship::getAirAbsViscosity(const double& tempCel)
{
    constexpr double mu0{ 1.716e-5 };
    constexpr double T0{ 273.11 };
    constexpr double S{ 100.56 };

    const double T{ wave::celsiusToKelvin(tempCel) };

    return mu0 * pow(T / T0, 1.5) * (T0 + S) / (T + S);
}
*/

// Reference: https://en.wikipedia.org/wiki/Density_of_air
// At IUPAC standard temperature and pressure (0 °C and 100 kPa), dry air has a density of approximately 1.2754 kg/m3.
double ship::airDensityAtTemp(const double& tempCel)
{
    // Reference: https://en.wikipedia.org/wiki/Atmospheric_pressure
    // Sea level is approximately 1 atm, a unit of pressure defined as 101,325 Pa
    constexpr double p{101325};

    constexpr double kB{1.380649e-23};
    constexpr double m{4.81e-26};

    return p * m / (kB * wave::celsiusToKelvin(tempCel));
}

double ship::getWaveR(const int& knots, const double& waveDegRelShip, const double& swh, const double& pp1d, int gammaIdx) const
{
    if (waveDegRelShip > 180)
    {
        cout << "rel angle should be <= 180; the other half uses symmetry" << endl;
        return -1;
    }

    if (gammaIdx < 0)
        gammaIdx = 0;
    else if (gammaIdx >= static_cast<int>(wave::jsGammaVec.size()))
        gammaIdx = static_cast<int>(wave::jsGammaVec.size()) - 1;

    const double waveDegRelShipLb = waveDegRelShip; //for single wave resistance

    wave w(this, knots, 1, swh, pp1d, waveDegRelShipLb, wave::angleDefaultSize, wave::angleDefaultStep);
    w.initSpectrum(gammaIdx);

    w.enableLogger("wave_run.csv", /*truncate=*/false);
    return w.getR();
}

double ship::getWindR_woAirDensity(const int& knots, const double& windDegRelShip, const double& windMag) const
{
    if (windDegRelShip > 180)
    {
        cout << "rel angle should be <= 180; the other half uses symmetry" << endl;
        return -1;
    }

    wind wind(this, knots, 1, windMag, windDegRelShip);
    return wind.getR_woAirDensity();
}

double ship::cappedWith_Rcalm(const double& upRatio, const double& downRatio, const double& Rcalm, const double& value)
{
    const double upValue{upRatio * Rcalm};
    const double downValue{-downRatio * Rcalm};

    if (value > upValue)
        return upValue;
    else if (value < downValue)
        return downValue;
    else
        return value;
}

//knots/angle/temp -> int, windMag/swh/pp1d -> double
double ship::fcTonsPerNm(int knots, const int& tempCel, const double& shipDeg_toConvention, const double& windDeg_toConvention, const double& windMag,
    const double& waveDeg_toConvention, const double& swh, const double& pp1d, int gammaIdx, double effIn)
{
    if (knots < minSpeed)
        knots = minSpeed;
    else if (knots > maxSpeed)
        knots = maxSpeed;

    hmInit();
    R_calm = hm84_Fn_lt04(knots, tempCel);
    
    const double waveRelShip_frConvention{wave::relAngleDegMod180_frConvention(waveDeg_toConvention, shipDeg_toConvention)};
    R_wave = getWaveR(knots, waveRelShip_frConvention, swh, pp1d, gammaIdx);

    const double windDegRelShip_fromConvention{wave::relAngleDegMod180_frConvention(windDeg_toConvention, shipDeg_toConvention)};

    R_wind = ship::airDensityAtTemp(static_cast<double>(tempCel)) * getWindR_woAirDensity(knots, windDegRelShip_fromConvention, windMag);

    double R_total = R_calm +
        ship::cappedWith_Rcalm(ship::waveUpRatioCap, ship::waveDownRatioCap, R_calm, R_wave) +
        ship::cappedWith_Rcalm(ship::windUpRatioCap, ship::windDownRatioCap, R_calm, R_wind);

    // Reference: [8] Kim, Y., Steen, S., Kramel, D., Muri, H., & Strømman, A. H. (2023). Modelling of ship resistance and power consumption for the global fleet: The MariTEAM model. Ocean Engineering, 281, 114758. https://doi.org/10.1016/j.oceaneng.2023.114758

    cout << "nRotorF " << nRotorF << " nSail " << nSail << endl;

    //sct integration
    if (nRotorF > 0)
    {
        const double windFromDeg = wave::swapAngleConvention(windDeg_toConvention);

        rotorF rotorObj(nRotorF);

        const auto [Vapp_mps, windFromAppDeg] = rotorF::apparent_from_true(windMag, wave::knotsToMeterPerSec(knots), shipDeg_toConvention, windFromDeg);

        const double f = rotorScale * rotorObj.getForce(Vapp_mps, shipDeg_toConvention, windFromAppDeg);
        const double p = rotorScale * rotorObj.getPowerConsumption(Vapp_mps);

        cout << "rotor f = " << f << " p = " << p << endl;

        P_total = (R_total - f) * wave::meterPerSecOfOneKnot + p;
    }
    else if (nSail > 0)
    {
        sail sailObj;

        const double windFromDeg = wave::swapAngleConvention(windDeg_toConvention);

        const double f = sailScale * sailObj.getForce(wave::meterPerSecToKnots(windMag), knots, shipDeg_toConvention, windFromDeg, nSail);

        cout << "sail f = " << f << endl;

        P_total = (R_total - f) * wave::meterPerSecOfOneKnot;
    }
    else
    {
        P_total = R_total * wave::meterPerSecOfOneKnot;      //Eq 15 of [8], no knots give fcTonPerNm
    }

    prop p(this, knots, tempCel);
    propEff = p.getEff(max(R_calm, R_total));

    if (effIn > 1 || effIn < wave::ep)
    {
        effIn = propEff;
    }

    return P_total * 1e-9 * sfoc / effIn;      // Power / 1000 to get kWh, multiply by SFOC (g/kWh) to get fc_gPerHr * 1e-6
}

double ship::fcTonsPerHr(int knots, const int& tempCel,
                         const double& shipDeg_toConvention, const double& windDeg_toConvention, const double& windMag,
                         const double& waveDeg_toConvention, const double& swh, const double& pp1d,
                         const int& gammaIdx, double effIn)
{
    const double fc_per_nm = fcTonsPerNm(knots, tempCel,
                                         shipDeg_toConvention, windDeg_toConvention, windMag,
                                         waveDeg_toConvention, swh, pp1d,
                                         gammaIdx, effIn);
    return knots * fc_per_nm;
}

// Scaling fuel consumption by displacement
// Reference: [9] MAN Diesel & Turbo. (n.d.) Basic Principles of Ship Propulsion.
double ship::fcDwtScaling(const double& newLoadDivOld) const
{
    if (newLoadDivOld < wave::ep || newLoadDivOld >= 1)
        return -1;

    const double dlRatio = shipTyChar == 'c' ? 6 : 3;                       // Table 2, page 7 of [9]
    return pow((1 + newLoadDivOld * dlRatio) / (1 + dlRatio), 2 * wave::oneThird);   // page 13 of [9]
}

// Scaling fuel consumption by draught
// Reference: [9] MAN Diesel & Turbo. (n.d.) Basic Principles of Ship Propulsion.
double ship::fcDraftScaling(const double& newDraftDivOld) const
{
    if (newDraftDivOld < wave::ep || newDraftDivOld >= 1)
        return -1;

    return pow(newDraftDivOld, 2 * wave::oneThird);                      // page 13 of [9]
}

bool ship::paramOk()
{
    bool ok{true};

    auto check = [&](bool condition, const string& message)
    {
        if (!condition)
        {
            ok = false;
            cout << "Error: " << message << endl;
        }
    };

    // checks on unused value
    check(dwt >= 0,             "DWT < 0");
    check(fuelCap >= 0,             "FuelCap < 0");
    check(!label.empty(),   "ShipLabel is empty");

    // bound checks
    check(maxSpeed <= ship::maxKnotsAllShip,   "MaxSpeed > " + to_string(ship::maxKnotsAllShip));
    check(Tmax <= ship::Tmax_limU,                  "MaxDraft > " + to_string(ship::Tmax_limU));
    check(B <= ship::B_limU,                          "Beam > limit" + to_string(ship::B_limU));
    check(Loa <= ship::Loa_limU,                    "LenOverall > " + to_string(ship::Loa_limU));
    check(displaceVol <= ship::displaceVol_limU,      "DisplaceVol > " + to_string(ship::displaceVol_limU));
    check(Abt <= ship::Abt_limU,                    "BulbousBowTransverseArea > " + to_string(ship::Abt_limU));
    check(sfoc <= ship::sfoc_limU,                    "SFOC > " + to_string(ship::sfoc_limU));
    check(ks <= ship::ks_limU,                        "ks > " + to_string(ship::ks_limU));
    //check(Cs <= ship::Cs_limU,                        "Cs > " + to_string(ship::Cs_limU));

    //check(windBallastCurve <= 1,                        "WindBallastCurve > 1");
    check(nRotorF <= ship::maxRotorF,                        "nRotorF > " + to_string(ship::maxRotorF));
    check(rotorScale <= 1,                        "RotorScale > 1");
    check(nSail <= ship::maxSails,                        "nSails > " + to_string(ship::maxSails));
    check(sailScale <= 1,                        "SailsScale > 1");

    check(minSpeed >= ship::minKnotsAllShip,   "MinSpeed < " + to_string(ship::minKnotsAllShip));
    check(Tmax >= ship::Tmax_limL,                  "MaxDraft < " + to_string(ship::Tmax_limL));
    check(B >= ship::B_limL,                          "Beam < " + to_string(ship::B_limL));
    check(Loa >= ship::Loa_limL,                    "LenOverall < " + to_string(ship::Loa_limL));
    check(displaceVol >= ship::displaceVol_limL,      "DisplaceVol < " + to_string(ship::displaceVol_limL));
    check(Abt >= ship::Abt_limL,                    "BulbousBowTransverseArea < " + to_string(ship::Abt_limL));
    check(sfoc >= ship::sfoc_limL,                    "SFOC < " + to_string(ship::sfoc_limL));
    check(ks >= 0,                        "ks < 0");
    //check(Cs >= 0,                        "Cs < 0");
    //check(windBallastCurve >= 0,                        "windBallastCurve < 0");
    check(nRotorF >= 0,                        "nRotorF < 0");
    check(rotorScale >= 0,                        "RotorScale < 0");
    check(nSail >= 0,                        "nSail < 0");
    check(sailScale >= 0,                        "SailScale < 0");

    // logic checks
    check(Lpp >= 0.9 * Loa,           "LenPerpendiculars < 0.9 * LenOverall");
    check(T <= Tmax,                   "Draft > MaxDraft");
    check(maxSpeed >= minSpeed,         "MaxSpeed < MinSpeed");
    check(nRotorF <= 0 || nSail <= 0,         "Either nRotorF or nSails > 0, but not both");

    check(serviceSpeed >= minSpeed,         "ServiceSpeed < MinSpeed");
    check(maxSpeed >= serviceSpeed,         "MaxSpeed < ServiceSpeed");

    const double Fn{maxSpeed * wave::meterPerSecOfOneKnot / sqrt(max(wave::ep, Lpp) * wave::gravityConst)};
    check(Fn <= 0.55, "Error: FroudeNum > 0.55");

    //----inputs, using these ratio simplifies
    Tf = T;
    hb = 0.6 * Tf; //handymax 7 / 12 = 0.58, h_B must be < T_F
    Fni_fac = Tf - hb - 0.25 * sqrt(Abt);

    check(Fni_fac >= wave::ep, "Tf - hb - 0.25 * sqrt(Abt): <= 0");

    //if (Fni_fac < wave::ep)
        //Fni_fac = wave::ep; //prevents neg & div by 0

    //if (windBallastCurve == 1 && T == Tmax)
        //cout << "When windBallastCurve = 1, draft should be < maxDraft" << endl;

    //prop var, we do it rather than hmInit() cos prop may be called before hmInit()
    // Reference: [4] Kristensen, Lutzen. (2013). Prediction of Resistance and Propulsion Power of Ships
    // Page 3 / 8
    const double ldr{Lpp / pow(displaceVol, wave::oneThird)};
    prop_diameter = prop::diameter(shipTyChar, Tmax);

    Cb = displaceVol / (Lpp * B * T); //must come before wakeFrac, thrustDeductionFac

    prop_w = prop::wakeFrac(this, Cb, ldr, prop_diameter);
    prop_t = prop::thrustDeductionFac(this, Cb, ldr, prop_diameter);

    return ok;
}

string ship::getComments() const
{
    return comments;
}

string ship::getLabel() const
{
    return label;
}

double ship::getDraft() const
{
    return T;
}

double ship::getSFOC() const
{
    return sfoc;
}

int ship::getDwt() const
{
    return dwt;
}

int ship::getFuelCap() const
{
    return fuelCap;
}

int ship::getMinSpeed() const
{
    return minSpeed;
}

int ship::getMaxSpeed() const
{
    return maxSpeed;
}

int ship::getServiceSpeed() const
{
    return serviceSpeed;
}

int ship::getNRotorF() const { return nRotorF; }
int ship::getNSail() const { return nSail; }

// setters that invalidate rotor/sail tables
void ship::setNRotorF(const int& n)
{
    nRotorF = clamp(n, 0, ship::maxRotorF);

    if (nRotorF > 0)
        nSail = 0; // ensure only one assist type is active

    windAidVec.clear();
}

void ship::setNSail(const int& n)
{
    nSail = clamp(n, 0, ship::maxSails);

    if (nSail > 0)
        nRotorF = 0; // ensure only one assist type is active

    windAidVec.clear();
}

void ship::setRotorScale(const double& s)
{
    if (s <= 0.0)
        rotorScale = wave::ep;
    else if (s >= 1.0)
        rotorScale = 1.0;
    else
        rotorScale = s;
}

void ship::setSailScale(const double& s)
{
    if (s <= 0.0)
        sailScale = wave::ep;
    else if (s >= 1.0)
        sailScale = 1.0;
    else
        sailScale = s;
}

pair<int, int> ship::idxOfUniformTable(const double& v, const double& lb, const int& size, const double& step)
{
    //size = table size, size of 1 -> singleton table, value is always first idx
    //see wind::Cda
    if (v <= lb || size <= 1) //<= lb return the first idx (0)
        return {0, 0};

    const int lastIdx{size - 1};

    if (v >= lb + lastIdx * step) //>= ub return the last idx
        return {lastIdx, lastIdx};

    const double idxDouble = (v - lb) / step;           //double ver of idx
    const int idx = static_cast<int>(floor(idxDouble)); //idx should be int, so double -> int using floor cos interpolate with idx & idx + 1

    if (idx == lastIdx)
        return {lastIdx, lastIdx};
    else
        return {idx, idx + 1};
}

void ship::makeWindAidVec(const int& _knotsLb, const int& _knotsSize)
{
    if (!windAidVec.empty())
        return; // build tables once

    if (nRotorF <= 0 && nSail <= 0)
    {
        cout << "makeWindAide() prob" << endl;
        windAidVec.clear();
        return;
    }

    windKnotsLb = _knotsLb;
    windKnotsSize = _knotsSize;

    windAidVec.reserve(ship::windAid_windMagSize);

    const unique_ptr<rotorF> rotorPtr = nRotorF > 0 ? make_unique<rotorF>(nRotorF) : nullptr;
    const unique_ptr<sail> sailPtr = nSail > 0 ? make_unique<sail>() : nullptr;

    for (int wi = 0; wi < ship::windAid_windMagSize; wi++)
    {
        const double wMag_mps{wi * windAid_windMagStep};

        windAidVec.emplace_back(rotorF::shipDegSize, rotorF::windDegSize, _knotsSize);
        auto& tab = windAidVec.back();

        for (int si = 0; si < rotorF::shipDegSize; si++)
        {
            const double shipDeg_to{si * rotorF::shipDegStep};

            for (int ai = 0; ai < rotorF::windDegSize; ai++)
            {
                const double windDeg_from{ai * rotorF::windDegStep};

                for (int ki = 0; ki < _knotsSize; ki++)
                {
                    const int shipSpeed_knots{_knotsLb + ki};
                    const double windSpeed_knots{wave::meterPerSecToKnots(wMag_mps)};

                    if (rotorPtr != nullptr)
                    {
                        // compute apparent wind then call 3-arg rotor API
                        const double ship_mps{shipSpeed_knots * wave::meterPerSecOfOneKnot};
                        const auto [Vapp_mps, windFromAppDeg] = rotorF::apparent_from_true(wMag_mps, ship_mps, shipDeg_to, windDeg_from);

                        const double f = rotorScale * rotorPtr->getForce(Vapp_mps, shipDeg_to, windFromAppDeg);
                        const double p = rotorScale * rotorPtr->getPowerConsumption(Vapp_mps); //assumption, power may not scale linearly

                        tab.force[si][ai][ki] = max(0.0, f);
                        tab.power[si][ai][ki] = max(0.0, p);
                    }
                    else if (sailPtr != nullptr)
                    {
                        const double f = sailScale * sailPtr->getForce(windSpeed_knots, shipSpeed_knots, shipDeg_to, windDeg_from, nSail);

                        tab.force[si][ai][ki] = max(0.0, f);
                    }
                }
            }
        }
    }
}

//sct, not tested rigourous, if has problem revert to old ver
pair<double, double> ship::getWindAid_byTable(const int& knots, const double& shipDeg_to, const double& windDeg_from, const double& windMag)
{
    //sct integration
    if (windAidVec.empty())
    {
        cout << "need to call makeAllVec()" << endl;
        return {-1, -1};
    }

    //if (nRotorF == 0 && nSail == 0)
    //return {0.0, 0.0};

    const auto sIdx{ship::idxOfUniformTable(shipDeg_to, 0, rotorF::shipDegSize, rotorF::shipDegStep)};
    const auto aIdx{ship::idxOfUniformTable(windDeg_from, 0, rotorF::windDegSize, rotorF::windDegStep)};

    const int kIdx{ship::idxOfKnotsTable(knots, windKnotsLb)}; //knots is passed in as int, so no point interpolate

    //2d intepolate on force, power
    const auto interpolatePair = [&](const int& wIdx)->pair<double, double>
    {
        const auto& tab = windAidVec[wIdx];

        const auto interpolate = [&](const d3& T)->double
        {
            if (sIdx.first == sIdx.second && aIdx.first == aIdx.second) //4 cases, return extremeValue, 2d, 1d (first var, second var)
            {
                return T[sIdx.first][aIdx.first][kIdx]; //extreme value
            }
            else if (sIdx.first != sIdx.second && aIdx.first != aIdx.second) //2d interpolate
            {
                const double x0{sIdx.first * rotorF::shipDegStep};
                const double x1{sIdx.second * rotorF::shipDegStep};

                const double y0{aIdx.first * rotorF::windDegStep};
                const double y1{aIdx.second * rotorF::windDegStep};

                const double g00{T[sIdx.first][aIdx.first][kIdx]};
                const double g01{T[sIdx.second][aIdx.first][kIdx]};
                const double g10{T[sIdx.first][aIdx.second][kIdx]};
                const double g11{T[sIdx.second][aIdx.second][kIdx]};

                return wave::interpolate2d(shipDeg_to, windDeg_from, g00, g01, g10, g11, x0, x1, y0, y1);
            }
            else
            {
                const bool interpolWind = aIdx.first != aIdx.second ? true : false;  //1d interpolate

                if (interpolWind)
                {
                    const double x0{aIdx.first * rotorF::windDegStep};
                    const double x1{aIdx.second * rotorF::windDegStep};

                    const double y0{T[sIdx.first][aIdx.first][kIdx]};
                    const double y1{T[sIdx.first][aIdx.second][kIdx]};

                    return wave::interpolate(x0, x1, y0, y1, windDeg_from);
                }
                else
                {
                    const double x0{sIdx.first * rotorF::shipDegStep};
                    const double x1{sIdx.second * rotorF::shipDegStep};

                    const double y0{T[sIdx.first][aIdx.first][kIdx]};
                    const double y1{T[sIdx.second][aIdx.second][kIdx]};

                    return wave::interpolate(x0, x1, y0, y1, shipDeg_to);
                }
            }
        };

        if (nRotorF > 0)
            return {interpolate(tab.force), interpolate(tab.power)};
        else if (nSail > 0)
            return {interpolate(tab.force), 0};
        else
            return {0, 0};
    };

    const auto wIdx = ship::idxOfUniformTable(windMag, 0, ship::windAid_windMagSize, ship::windAid_windMagStep); //0 = lb

    const double wLowMag = wIdx.first * ship::windAid_windMagStep;
    const double wHighMag = wIdx.second * ship::windAid_windMagStep;

    if (wIdx.first == wIdx.second || std::fabs(wHighMag - wLowMag) < wave::ep)
    {
        const pair<double, double> v = interpolatePair(wIdx.first);

        return {v.first, v.second};
    }
    else
    {
        const pair<double, double> wLow = interpolatePair(wIdx.first);
        const pair<double, double> wHigh = interpolatePair(wIdx.second);

        const double f = wave::interpolate(wLowMag, wHighMag, wLow.first, wHigh.first, windMag);
        const double p = wave::interpolate(wLowMag, wHighMag, wLow.second, wHigh.second, windMag);

        return {f, p};
    }
}

/*
double ship::fcTonsRoute(const i1& knots, const i1& tempCel, const d1& shipDeg_toConvention, const d1& windDeg_toConvention, const d1& windMag,
    const d1& waveDeg_toConvention, const d1& swh, const d1& pp1d, const d1& time, const i1& gammaIdx, const double eff)
{
    double fcTons{0};

    for (size_t i = 0; i < time.size(); i++)
    {
        const double f = fcTonsPerNm(knots[i], tempCel[i], shipDeg_toConvention[i], windDeg_toConvention[i], windMag[i], waveDeg_toConvention[i], swh[i], pp1d[i], gammaIdx[i], eff);

        if (f > 0)
            fcTons += time[i] * f;
    }

    return fcTons;
}
*/
