#include "csvFn.h"
#include "shipV.h"
#include "wave.h"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <algorithm>

using namespace std;

//22 Oct
//A) simulation.csv shipKnots are all 16. Previous ver has wide range -> byTable is slow

//10 Oct
//A) waveDegRelShipFrConvention is used by calcDirSpectrum. I think this is all that is needed, look for double aa = alphaVec[i];
//once the above is firm, comment out the below of wave.h (no longer needed)
//double waveDegRelShipFrConvention{0};
//B) ammended, look for // Change 2:
//C) wind::Axv uses map. see my query under the reference, just to make sure
//D) added wave::initSpectrum, this is to get rid of 'constructor calc for i == 0' below (which is cumbersome)
//for (int i = 0; i < static_cast<int>(wave::jsGammaVec.size()); i++)
//{
    //if (i > 0)  //constructor calc for i == 0
        //calcFreqSpectrum_sct(i);

    //for (int j = 0; j < waveDegRelShipSize; j++)
    //{
        //if (i > 0 && skipWaveDeg[j]) //i == 0 is PM spectrum, we do not skip
            //continue;

        //if (i > 0 && j > 0) //constructor calc for i == 0 && j == 0
        //{
            //const double waveDeg{ waveDegRelShipLb + j * waveDegRelShipStep };
            //wave::calcDirSpectrum(waveDeg);
        //}
//E) wave::Rcth_aw differs a lot for first few pts of simulation.csv. Last few points are the same. You decide which ver to use

//Brian's ver
    //R_Wave: 126212
    //R_Wave: 122804
    //R_Wave: 114544

//sct's ver
    //R_Wave: 77597.9
    //R_Wave: 75845.5
    //R_Wave: 73891.5

//F) Fix item F) under 3 Oct

//3 Oct
//D) to make coding simple, we use fcTonsPerNm(..., gammaIdx), arg corresponds to idx of wave::jsGammaVec{1.0, 3.3}. Usage example:
//double userGamma = 5.5
//gammaIdx = wave::nearestIdxOfSortedVec(wave::jsGammaVec, userGamma)
//fcTonsPerNm(..., gammaIdx)
//fcTonsPerNm(..., -1), fcTonsPerNm(..., 100) are fine, because out of range value are clipped
//PM spec is just JS with gamma = 1, so wave::jsGammaVec{1.0, ...} always starts with 1.0, user can just use gammaIdx = 0 for PM

//F) the below has div by zero problem, though R_aw avoids it
//double wave::Rlp_lambda(const double& ww)
//{
//    return wave::twoPi / wave::waveNum(ww);
//}

//double wave::R_aw(const double& aa, const double& ww, const R_aw_choice& R_aw_method) const
//{
// If frequency is zero, no time-varied waves -> resistance = 0
//if (fabs(ww) < wave::ep)
//return 0;

//solution: wave::waveNum to have a lower bound of say 0.001 wave per meter
//check whether there are more of such singularity (I dun think there is)

//H) prev calcDirSpectrum does not look right, new ver is clearer, note (int i = iStart; i <= iEnd; i++), includes both ends, like freq spectrum
//Rwave is higher than previous
//I retain storing alphaVec as radian, leave Brian to store in deg. When store in radian, we end up with the below radToDeg and then degToRad fn calls
//const double angleRad = wave::degToRad(wave::relAngleDegMod180(wave::radToDeg(alphaVec[i]), 0.0));
//wave::relAngleDegMod180 is named to indicate result will be [0, 180] deg. Clearer convention

//Brian to do the following (most important comes first)
//1)
//Check whether Rlp, Rcth use the same convention for alphaVec. If each model uses different convention, need to tidy up out the mess (& document it). Example as below.

// Combined Wave Model
// Literature uses to_convention : 180 indicates head waves
//double wave::g(const double& aa, const double& bb, const double& vv)

//3) I tuned freq spectrum nBins to be efficient, using freq 'distribution'. (It skip points at tail end).
//For dir spectrum, this is not possible, so we stick to uniform bins. Find good nBins will do, currently wave::angleDefaultStep{45}
//4) wave::getR() uses all 3 methods for logging. When not needed, use only the combined method.
//5) Once the above are all done, compare the below ships against their respective paper
//Panamax container: Bunker consumption of containerships considering sailing speed and wind conditions
//RoRo: On the estimation of ship’s fuel consumption and speed curve: A statistical approach

//6) Also compare the ships in L&P, Lang&Mao. These are easier than the above

//---------

//14 jun
//For info
//h) at prop.cpp, double fCth = min(0.65, 0.81 - 0.014 * C_th); min is wrong, should be max -> already fixed
//also, hullEff is around 1.27, need to limit it to 1 -> already fixed
//using Handymax Paper, A_E / A_O = 0.53, pitch ratio P/D = 0.74
//https://enghandbook.com/calculators/b-series-propeller-open-water-curves/ gives about 0.63 max
//https://github.com/mkergoat/bseries

//8)
//Ref [2] pg 415: Thus, in case of high required thrust levels, particularly for low-speed design sce-
//narios, i.e. high CTh coefficients, the installation of two propellers is preferable, be-
//cause this way the resulting CTh is reduced and the efficiency η0 increases.

//Perhaps for twin properller, can use CTh / 2 to calc eff?

//A) ballast condition
//Paper 15 says that ballast draft is 0.6, can test different to model ballast condition, via draft or admiralty coefficient
//IITC doc may have different CX figures for ballast for wind (use shipParam.csv new flag windBallastCurve)

//B) other ships, see Ref [2], pg 51, Fig 1.2 onwards
//a) panamax container
//b) roro
//c) tankers

//C) Quantifying voyage optimisation with wind propulsion for short-term CO2 mitigation in shipping
//Eq 3, SFOC as function of engine load

//D) Impact of Wind-Assisted Propulsion on Fuel Savings and Propeller Efficiency: A Case Study
//Figure 3, 4 - prop eff
//Figure 7 - counterpart of rotorF.h for sails, can add another class
//the sail size was approximated as three 1900 m2 sails
//Figure 12 - see can get back the savings figure
//https://blog.3ds.com/industries/marine-offshore/the-return-of-wind-assisted-propulsion-at-sea/
//The company plans to install two 76-meter masts on a cargo ship, achieving a total sail area of 3000 m2

//G) Useful info from ref [2]
    //i) pg 57, Eq 1.11
    //ii) pg 70, Table 2.1
    //iii) pg 178 Eq 2.105

//H) double ship::S_byMumford() const has limited ship, look for more

//I) Because flettner rotor consumes power, it is better to off it there is no wind. Just can use fc = min(withRotor, withoutRotor)

int main(int argc, char** argv)
{
    string shipFile = "shipHandy.csv"; // default ship list
    string wayPtFile = "simulation.csv"; // default route; pass a CSV to override

    if (argc >= 3)
    {
        shipFile = argv[1];
        wayPtFile = argv[2];
    }

    if (!filesystem::exists(shipFile))
    {
        std::cout << shipFile << " does not exist" << endl;
        return 1;
    }

    if (!filesystem::exists(wayPtFile))
    {
        std::cout << wayPtFile << " does not exist" << endl;
        return 1;
    }

    const csvFn sMat(wayPtFile);
    const auto& wayPts = sMat.getWaypoints();

    //----

    shipV shipList(shipFile, 0); // start with first row; we will switch rows below

    const int nRows = shipList.getRowCount();
    if (nRows <= 0)
    {
        std::cout << "No ships loaded from " << shipFile << endl;
        return 1;
    }


    auto avg_for_row = [&](int rowIdx) -> pair<string, double>
    {
        shipList.setNewRowIdx(rowIdx);
        ship* s = shipList.shipPtr();
        if (!s)
            return {string("Row ") + to_string(rowIdx), 0.0};

        // Optional: prebuild tables for rotor/sail cases to speed lookups
        // Build knots range from waypoints
        vector<int> shipKnots;
        shipKnots.reserve(wayPts.size());
        for (const auto& it : wayPts)
        {
            const int k = static_cast<int>(std::lround(it.shipKnots));
            if (k > 0) shipKnots.emplace_back(k);
        }
        if (!shipKnots.empty())
        {
            const int knotsLb = *min_element(shipKnots.begin(), shipKnots.end());
            const int knotsSize = *max_element(shipKnots.begin(), shipKnots.end()) - knotsLb + 1;
            s->makeTable(knotsLb, knotsSize);
        }

        double sum{0};
        int count{0};

        for (size_t wpIdx = 0; wpIdx < wayPts.size(); wpIdx++)
        {
            const auto& wp = wayPts[wpIdx];
            const int knots = static_cast<int>(std::round(wp.shipKnots));
            if (knots <= 0) continue;

            const double shipDeg_toConvention = wp.shipBearing;
            const double windDeg_toConvention = wp.windBearing;
            const double waveDeg_toConvention = wp.waveBearing;
            const int tempCel = static_cast<int>(std::round(wp.temp2m));
            const double windMag = wp.windMag;
            const double swh = wp.swh;
            const double pp1d = wp.pp1d;
            const int gammaIdx = 0; // PM (idx of 1.0 in jsGammaVec)

            const double fcTonsPerNm = s->fcTonsPerNm(knots, tempCel, shipDeg_toConvention, windDeg_toConvention, windMag, waveDeg_toConvention, swh, pp1d, gammaIdx);
            const double fcTonsPerHr = fcTonsPerNm * knots;

            if (fcTonsPerHr > wave::ep)
            {
                sum += fcTonsPerHr;
                count++;
            }
        }

        const string label = s->getLabel();
        const double avgHr = count > 0 ? sum / count : 0.0;
        std::cout << "Ship: " << label << "   wayPtCount: " << count
                  << "   Avg tonsPerHr: " << avgHr
                  << "   Avg tonsPerDay: " << (avgHr * 24.0) << endl;

        return {label, avgHr};
    };

    // Run for each row and report. If there are >=3 rows, assume 0=baseline, 1=rotor, 2=sail (shipHandy.csv layout)
    vector<pair<string, double>> results;
    for (int i = 0; i < nRows; i++)
        results.emplace_back(avg_for_row(i));

    if (!results.empty())
    {
        const double base = results.front().second;
        std::cout << "\nSummary (relative to baseline: " << results.front().first << ")" << endl;
        for (size_t i = 0; i < results.size(); i++)
        {
            const auto& [lbl, v] = results[i];
            if (i == 0)
            {
                std::cout << "  " << lbl << ": " << v << " tons/hr" << endl;
            }
            else
            {
                const double savingPct = base > wave::ep ? 100.0 * (base - v) / base : 0.0;
                std::cout << "  " << lbl << ": " << v << " tons/hr   Saving: " << savingPct << " %" << endl;
            }
        }
    }

    std::cout << endl << argv[0] << " " << shipFile << " " << wayPtFile << endl;

    return 0;
}
