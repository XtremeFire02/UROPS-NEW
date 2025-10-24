#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>

#include "ship.h"
#include "prop.h"
#include "wave.h"
#include "wind.h"

using namespace std;

void ship::makeTableGivenRange(const int& _knotsLb, const int& _knotsSize, const int& _tempLb, const int& _tempSize,
    const double& _windMagLb, const int& _windMagSize, const double& _windMagStep,
    const double& _windDegRelShipLb, const int& _windDegRelShipSize, const double& _windDegRelShipStep,
    const double& _swhLb, const int& _swhSize, const double& _swhStep,
    const double& _pp1dLb, const int& _pp1dSize, const double& _pp1dStep,
    const double& _waveDegRelShipLb, const int& _waveDegRelShipSize, const double& _waveDegRelShipStep)
{
    if (!tableEmpty())  //cannot call more than once
        return;

    if (_tempSize <= 0 || _windMagSize <= 0 || _swhSize <= 0|| _pp1dSize <= 0 || _knotsSize <= 0 || _windDegRelShipSize <= 0 || _waveDegRelShipSize <= 0)
    {
        cout << "size <= 0" << endl;
        return;
    }

    if (_windMagStep <= 0 || _swhStep <= 0|| _pp1dStep <= 0 || _windDegRelShipStep <= 0 || _waveDegRelShipStep <= 0)
    {
        cout << "step <= 0" << endl;
        return;
    }

    if (_tempLb < -wave::ep || _windMagLb < -wave::ep || _swhLb < -wave::ep || _pp1dLb < -wave::ep || _windDegRelShipLb < -wave::ep || _waveDegRelShipLb < -wave::ep)
    {
        cout << "lb < 0" << endl;
        return;
    }

    if (_knotsLb < minSpeed)
    {
        cout << "_knotsLb < minSpeed" << endl;
        return;
    }

    if (_knotsLb + _knotsSize - 1 > maxSpeed)
    {
        cout << "_knotsUb > minSpeed" << endl;
        return;
    }

    makeCalmPropTable(_knotsLb, _knotsSize, _tempLb, _tempSize);
    makeWindTable(_knotsLb, _knotsSize, _windMagLb, _windMagSize, _windMagStep, _windDegRelShipLb, _windDegRelShipSize, _windDegRelShipStep);
    makeWaveTable(_knotsLb, _knotsSize, _swhLb, _swhSize, _swhStep, _pp1dLb, _pp1dSize, _pp1dStep, _waveDegRelShipLb, _waveDegRelShipSize, _waveDegRelShipStep);

    //sct integration
    if (nRotorF > 0 || nSail > 0)
        makeWindAidVec(_knotsLb, _knotsSize);
}

void ship::clearTableIfKnotsExceedRange(const int& _knotsLb, const int& _knotsSize)
{
    if (tableEmpty())
        return;

    set<int> prevKnots, newKnots;

    for (int i = calmKnotsLb; i < calmKnotsSize; i++)
        prevKnots.insert(i);

    //for (int i = windKnotsLb; i < windKnotsSize; i++)
        //prevKnots.insert(i);

    //for (int i = waveKnotsLb; i < waveKnotsSize; i++)
        //prevKnots.insert(i);

    for (int i = _knotsLb; i < _knotsSize; i++)
        newKnots.insert(i);

    i1 diff;

    set_difference(newKnots.begin(), newKnots.end(), prevKnots.begin(), prevKnots.end(), inserter(diff, diff.begin()));

    if (!diff.empty()) //new knots not present previously are required -> clear previous table
    {
        propTable_R.clear();
        propTable_eff.clear();
        calmTable.clear();

        windTable.clear();
        waveTable.clear();
    }
}

//because of omp, using fullRange is still fast
void ship::makeTable(const int& _knotsLb, const int& _knotsSize)
{
    if (!tableEmpty())  //cannot call more than once
        return;

    if (_knotsLb < minSpeed)
    {
        cout << "_knotsLb < minSpeed" << endl;
        return;
    }

    if (_knotsLb + _knotsSize - 1 > maxSpeed)
    {
        cout << "_knotsUb > minSpeed" << endl;
        return;
    }

    const int tempSize{wave::waterTempCelTable.back() - wave::waterTempCelTable.front() + 1}; //about 30 deg but used only by calmProp, which is fast

    const int windMagSize{static_cast<int>(floor(wind::windMag_limU / wind::windMagDefaultStep)) + 1};
    const int wpSize{static_cast<int>(floor((wave::pp1d_limU - wave::pp1d_limL) / wave::pp1dDefaultStep)) + 1};
    const int whSize{static_cast<int>(floor(wave::swh_limU / wave::swhDefaultStep)) + 1};

    makeTableGivenRange(_knotsLb, _knotsSize, wave::waterTempCelTable.front(), tempSize,
        0.0, windMagSize, wind::windMagDefaultStep,
        0.0, wind::angleDefaultSize, wind::angleDefaultStep,
        0.0, whSize, wave::swhDefaultStep,
        wave::pp1d_limL, wpSize, wave::pp1dDefaultStep,
        0.0, wave::angleDefaultSize, wave::angleDefaultStep);
}

void ship::makeCalmPropTable(const int& _knotsLb, const int& _knotsSize, const int& _tempLb, const int& _tempSize)
{
    if (!calmTable.empty())  //cannot call more than once
        return;

    tempLb = _tempLb;
    tempSize = _tempSize;

    calmKnotsLb = _knotsLb;
    calmKnotsSize = _knotsSize;

    hmInit(); //must come before hm84_Fn_lt04()

    wave::makeVec2d(calmTable, calmKnotsSize, tempSize, 0.0);

    //#pragma omp parallel for collapse(2) num_threads(6)
    for (int i = 0; i < calmKnotsSize; i++)
        for (int j = 0; j < tempSize; j++)
            calmTable[i][j] = hm84_Fn_lt04(calmKnotsLb + i, tempLb + j);

    //prop
    wave::makeVec2d(propTable_eff, calmKnotsSize, prop::Rsize, 0.0);
    wave::makeVec2d(propTable_R, calmKnotsSize, prop::Rsize, 0.0);

    constexpr double frac{0.7};

    for (int i = 0; i < calmKnotsSize; i++)
    {
        prop p(this, calmKnotsLb + i, prop::assumedTemp);  //for simplicity, use fixed tempCel

        const double R_base{frac * hm84_Fn_lt04(calmKnotsLb + i, prop::assumedTemp)};

        //#pragma omp parallel for num_threads(6)
        for (int j = 0; j < prop::Rsize; j++)
        {
            const double R{R_base * (1 + j * prop::Rfrac)};
            const double eff{p.getEff(R)};

            propTable_eff[i][j] = eff;
            propTable_R[i][j] = R;
        }
    }
}

void ship::makeWindTable(const int& _knotsLb, const int& _knotsSize, const double& _windMagLb, const int& _windMagSize, const double& _windMagStep,
    const double& _windDegRelShipLb, const int& _windDegRelShipSize, const double& _windDegRelShipStep)
{
    if (!windTable.empty())  //cannot call more than once
        return;

    windMagLb = _windMagLb;
    windMagSize = _windMagSize;
    windMagStep = _windMagStep;

    windKnotsLb = _knotsLb;
    windKnotsSize = _knotsSize;
    windDegRelShipLb = _windDegRelShipLb;
    windDegRelShipSize = _windDegRelShipSize;
    windDegRelShipStep = _windDegRelShipStep;

    wave::makeVec2d(windTable, windMagSize, windDegRelShipSize, wind());

    for (int i = 0; i < windMagSize; i++)
    {
        for (int j = 0; j < windDegRelShipSize; j++)
            windTable[i][j].init(this, windKnotsLb, windKnotsSize, windMagLb + i * windMagStep, windDegRelShipLb + j * windDegRelShipStep);
    }

    //#pragma omp parallel for collapse(2) num_threads(6)
    for (int i = 0; i < windMagSize; i++)
    {
        for (int j = 0; j < windDegRelShipSize; j++)
            windTable[i][j].makeTable();
    }
}

#include <chrono>
void ship::makeWaveTable(const int& _knotsLb, const int& _knotsSize, const double& _swhLb, const int& _swhSize, const double& _swhStep, const double& _pp1dLb, const int& _pp1dSize, const double& _pp1dStep,
                         const double& _waveDegRelShipLb, const int& _waveDegRelShipSize, const double& _waveDegRelShipStep)
{
    if (!waveTable.empty())  //cannot call more than once
        return;

    waveR_atAngle.reserve(2);

    //must come before for loop
    swhLb = _swhLb;
    swhSize = _swhSize;
    swhStep = _swhStep;

    pp1dLb = _pp1dLb;
    pp1dSize = _pp1dSize;
    pp1dStep = _pp1dStep;

    waveKnotsLb = _knotsLb;
    waveKnotsSize = _knotsSize;
    waveDegRelShipLb = _waveDegRelShipLb;
    waveDegRelShipSize = _waveDegRelShipSize;
    waveDegRelShipStep = _waveDegRelShipStep;
    //

    wave::makeVec2d(waveTable, swhSize, pp1dSize, wave());

    for (int i = 0; i < swhSize; i++)
        for (int j = 0; j < pp1dSize; j++)
            waveTable[i][j].init(this, waveKnotsLb, waveKnotsSize, swhLb + i * swhStep, pp1dLb + j * pp1dStep, waveDegRelShipLb, waveDegRelShipSize, waveDegRelShipStep);
    
    using namespace std::chrono;
    const auto clockStart = steady_clock::now();

    //#pragma omp parallel for collapse(2) num_threads(6)
    //#pragma omp parallel for collapse(2) num_threads(6) schedule(dynamic,12) //freq nBins is dynamic
    for (int i = 0; i < swhSize; i++)
        for (int j = 0; j < pp1dSize; j++)
            waveTable[i][j].makeTable();

    const double time = duration_cast<duration<double> >(steady_clock::now() - clockStart).count();
    cout << "makeWaveVec time " << time << endl;
}

double ship::fcTonsPerNm_byTable(int knots, const int& tempCel, const double& shipDeg_toConvention, const double& windDeg_toConvention, const double& windMag,
    const double& waveDeg_toConvention, const double& swh, const double& pp1d, int gammaIdx, double effIn)
{
    //warning
    if (tableEmpty())
    {
        cout << "need to call makeAllVec()" << endl;
        return -1;
    }

    if (knots < minSpeed)
        knots = minSpeed;
    else if (knots > maxSpeed)
        knots = maxSpeed;
    //end warning

    //cout << "knots " << knots << " shipDeg " << shipDeg_toConvention << " tempCel " << tempCel << " windDeg " << windDeg_toConvention << " windMag " << windMag
    //<< " waveDeg" << waveDeg_toConvention << " swh " << swh << " pp1d " << pp1d << " gammaIdx " << gammaIdx << endl;

    const double R_calm{getCalmR_byTable(knots, tempCel)};

    const double waveRelShip_frConvention{wave::relAngleDegMod180_frConvention(waveDeg_toConvention, shipDeg_toConvention)};

    const double R_wave{getWaveR_byTable(knots, waveRelShip_frConvention, swh, pp1d, gammaIdx)};

    const double windDegRelShip_fromConvention{wave::relAngleDegMod180_frConvention(windDeg_toConvention, shipDeg_toConvention)};

    const double R_wind{ship::airDensityAtTemp(tempCel) * getWindR_woAirDensity_byTable(knots, windDegRelShip_fromConvention, windMag)};

    double R_total = R_calm +
        ship::cappedWith_Rcalm(ship::waveUpRatioCap, ship::waveDownRatioCap, R_calm, R_wave) +
        ship::cappedWith_Rcalm(ship::windUpRatioCap, ship::windDownRatioCap, R_calm, R_wind);

    //sct integration
    if (nRotorF > 0 || nSail > 0)
    {
        const double windFromDeg = wave::swapAngleConvention(windDeg_toConvention);

        const auto [f, p] = getWindAid_byTable(knots, shipDeg_toConvention, windFromDeg, windMag);

        cout << "byTable windAid f = " << f << " p = " << p << endl;

        if (nRotorF > 0)
            P_total = (R_total - f) * wave::meterPerSecOfOneKnot + p;
        else
            P_total = (R_total - f) * wave::meterPerSecOfOneKnot;
    }
    else
    {
        P_total = R_total * wave::meterPerSecOfOneKnot;      //Eq 15 of [8], no knots give fcTonPerNm
    }

    if (effIn > 1 || effIn < wave::ep)
        effIn = getProp_byTable(knots, R_total);

    cout << "\tR_Calm by table: " << R_calm << endl;
    cout << "\tR_Wave by table: " << R_wave << endl;
    cout << "\tR_Wind by table: " << R_wind << endl;

    cout << "\tR_Total by table: " << R_total << endl;
    cout << "\tProp efficiency: " << effIn << endl << endl;

    return P_total * 1e-9 * sfoc / effIn;      // Power / 1000 to get kWh, multiply by SFOC (g/kWh) to get fc_gPerHr * 1e-6 = fcTonsPerNm
}

bool ship::tableEmpty() const
{
    return propTable_R.empty() && propTable_eff.empty() && calmTable.empty() && windTable.empty() && waveTable.empty();
}

int ship::idxOfKnotsTable(const int& knots, const int& knotsLb)
{
    return knots - knotsLb;
}

double ship::getCalmR_byTable(const int& knots, const int& tempCel) const
{
    const int kIdx{ship::idxOfKnotsTable(knots, calmKnotsLb)};
    return calmTable[kIdx][tempCel - tempLb];
}

double ship::getProp_byTable(const int& knots, const double& R) const
{
    const int kIdx{ship::idxOfKnotsTable(knots, calmKnotsLb)};

    if (R <= propTable_R[kIdx].front())
        return propTable_eff[kIdx].front();
    else if (R >= propTable_R[kIdx].back())
        return propTable_eff[kIdx].back();
    {
        const int rIdx{static_cast<int>(floor((R - propTable_R[kIdx].front()) / (propTable_R[kIdx].front() * prop::Rfrac)))};

        const double x0{propTable_R[kIdx][rIdx]};
        const double x1{propTable_R[kIdx][rIdx + 1]};

        const double y0{propTable_eff[kIdx][rIdx]};
        const double y1{propTable_eff[kIdx][rIdx + 1]};

        return wave::interpolate(x0, x1, y0, y1, R);
    }
}

double ship::getWindR_woAirDensity_byTable(const int& knots, const double& windDegRelShip, const double& windMag)
{
    //warning
    if (windDegRelShip > 180)
    {
        cout << "rel angle should be <= 180; the other half uses symmetry" << endl;
        return -1;
    }

    const int kIdx{ship::idxOfKnotsTable(knots, windKnotsLb)};

    const auto mIdx{ship::idxOfUniformTable(windMag, windMagLb, windMagSize, windMagStep)};
    const auto aIdx{ship::idxOfUniformTable(windDegRelShip, windDegRelShipLb, windDegRelShipSize, windDegRelShipStep)};
    //----

    if (mIdx.first == mIdx.second && aIdx.first == aIdx.second)
    {
        return windTable[mIdx.first][aIdx.first].table[kIdx]; //extreme value
    }
    else if (mIdx.first != mIdx.second && aIdx.first != aIdx.second)
    {
        const double x0{windMagLb + mIdx.first * windMagStep};
        const double x1{windMagLb + mIdx.second * windMagStep};

        const double y0{windDegRelShipLb + aIdx.first * windDegRelShipStep};
        const double y1{windDegRelShipLb + aIdx.second * windDegRelShipStep}; //sct corrected aIdx.first -> aIdx.second

        const double g00{windTable[mIdx.first][aIdx.first].table[kIdx]};
        const double g01{windTable[mIdx.second][aIdx.first].table[kIdx]};
        const double g10{windTable[mIdx.first][aIdx.second].table[kIdx]};
        const double g11{windTable[mIdx.second][aIdx.second].table[kIdx]};

        return wave::interpolate2d(windMag, windDegRelShip, g00, g01, g10, g11, x0, x1, y0, y1);
    }
    else
    {
        const bool interpolAngle = aIdx.first != aIdx.second ? true : false;

        if (interpolAngle)
        {
            const double x0{windDegRelShipLb + aIdx.first * windDegRelShipStep};
            const double x1{windDegRelShipLb + aIdx.second * windDegRelShipStep};

            const double y0{windTable[mIdx.first][aIdx.first].table[kIdx]};
            const double y1{windTable[mIdx.first][aIdx.second].table[kIdx]};

            return wave::interpolate(x0, x1, y0, y1, windDegRelShip);
        }
        else
        {
            const double x0{windMagLb + mIdx.first * windMagStep};
            const double x1{windMagLb + mIdx.second * windMagStep};

            const double y0{windTable[mIdx.first][aIdx.first].table[kIdx]};
            const double y1{windTable[mIdx.second][aIdx.second].table[kIdx]};
            // Correction: aIdx.first -> aIdx.second
            //sct's reply: for interpolAngle == false, aIdx.first == aIdx.second, so no diff

            return wave::interpolate(x0, x1, y0, y1, windMag);
        }
    }
}

double ship::getWaveR_byTable(const int& knots, const double& waveDegRelShip, const double& swh, const double& pp1d, int gammaIdx)
{
    //warning
    if (waveTable.empty())
    {
        cout << "Wave Vec empty" << endl;
        return -1;
    }

    if (waveDegRelShip > 180)
    {
        cout << "rel angle should be <= 180; the other half uses symmetry" << endl;
        return -1;
    }

    if (gammaIdx < 0)
        gammaIdx = 0;
    else if (gammaIdx >= static_cast<int>(wave::jsGammaVec.size()))
        gammaIdx = static_cast<int>(wave::jsGammaVec.size()) - 1;
    //end warning

    //gammaIdx = 0;

    const int kIdx{ship::idxOfKnotsTable(knots, waveKnotsLb)};

    const auto hIdx{ship::idxOfUniformTable(swh, swhLb, swhSize, swhStep)};
    const auto pIdx{ship::idxOfUniformTable(pp1d, pp1dLb, pp1dSize, pp1dStep)};
    const auto aIdx{ship::idxOfUniformTable(waveDegRelShip, waveDegRelShipLb, waveDegRelShipSize, waveDegRelShipStep)};

    waveR_atAngle.clear();

    //cout << gammaIdx << endl;

    //we first interpolate on h & p, there are 4 cases
    //value is beyond range -> get the extreme value
    //2d interpolation on h, p
    //1d interpolation on h
    //1d interpolation on p

    //cout << "swh " << swh << " hIdx " << hIdx.first << endl;

    if (hIdx.first == hIdx.second && pIdx.first == pIdx.second)
    {
        //cout << "extreme" << endl;

        for (int a = aIdx.first; a <= aIdx.second; a++)
            waveR_atAngle.emplace_back(waveTable[hIdx.first][pIdx.first].table[gammaIdx][aIdx.first][kIdx]); //extreme value
    }
    else if (hIdx.first != hIdx.second && pIdx.first != pIdx.second)
    {
        for (int a = aIdx.first; a <= aIdx.second; a++)
        {
            const double x0{swhLb + hIdx.first * swhStep};
            const double x1{swhLb + hIdx.second * swhStep};

            const double y0{pp1dLb + pIdx.first * pp1dStep};
            const double y1{pp1dLb + pIdx.second * pp1dStep};

            const double g00{waveTable[hIdx.first][pIdx.first].table[gammaIdx][a][kIdx]};
            const double g01{waveTable[hIdx.second][pIdx.first].table[gammaIdx][a][kIdx]};
            const double g10{waveTable[hIdx.first][pIdx.second].table[gammaIdx][a][kIdx]};
            const double g11{waveTable[hIdx.second][pIdx.second].table[gammaIdx][a][kIdx]};

            waveR_atAngle.emplace_back(wave::interpolate2d(swh, pp1d, g00, g01, g10, g11, x0, x1, y0, y1));
        }
    }
    else
    {
        const bool interpolPeriod = pIdx.first != pIdx.second ? true : false;

        if (interpolPeriod)
        {
            for (int a = aIdx.first; a <= aIdx.second; a++)
            {
                const double x0{pp1dLb + pIdx.first * pp1dStep};
                const double x1{pp1dLb + pIdx.second * pp1dStep};

                const double y0{waveTable[hIdx.first][pIdx.first].table[gammaIdx][a][kIdx]};
                const double y1{waveTable[hIdx.first][pIdx.second].table[gammaIdx][a][kIdx]};

                waveR_atAngle.emplace_back(wave::interpolate(x0, x1, y0, y1, pp1d));
            }
        }
        else
        {
            for (int a = aIdx.first; a <= aIdx.second; a++)
            {
                const double x0{swhLb + hIdx.first * swhStep};
                const double x1{swhLb + hIdx.second * swhStep};

                const double y0{waveTable[hIdx.first][pIdx.first].table[gammaIdx][a][kIdx]};
                const double y1{waveTable[hIdx.second][pIdx.first].table[gammaIdx][a][kIdx]};

                waveR_atAngle.emplace_back(wave::interpolate(x0, x1, y0, y1, swh));
            }
        }
    }

    //interpolate on angle
    if (waveR_atAngle.size() == 1)
    {
        return waveR_atAngle.front();
    }
    else
    {
        double x0{waveDegRelShipLb + aIdx.first * waveDegRelShipStep};
        double x1{waveDegRelShipLb + aIdx.second * waveDegRelShipStep};

        double y0{waveR_atAngle.front()};
        double y1{waveR_atAngle.back()};

        return wave::interpolate(x0, x1, y0, y1, waveDegRelShip);
    }
}
