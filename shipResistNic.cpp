#include "csvFn.h"
#include "shipV.h"
#include "wave.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

using namespace std;

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
    //shipV::printShipCsv();
    //return 1;

    string shipFile = "shipHandy.csv";
    string wayPtFile = "simulation.csv"; // "simulation.csv";

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

	//Irrelevant code commented out
    //shipV shipList(shipFile, 0);
    //const int nRows = shipList.getRowCount();

    //const int idx{0}; //use row0 of file

    //NEW CODE START (To select ship row from shipHandy file)
    int rowIdx = 0;
    if (argc >= 4) {
        rowIdx = std::stoi(argv[3]);
    }
    shipV shipList(shipFile, 0);
    const int nRows = shipList.getRowCount();
    if (rowIdx < 0 || rowIdx >= nRows) {
        std::cout << "Invalid row index\n";
        return 1;
    }
    const int idx{ rowIdx };
    //NEW CODE END

    shipList.setNewRowIdx(idx);
    std::cout << "Selected ship row: " << idx << "\n";

    ship* shipP = shipList.shipPtr();

    if (shipP == nullptr)
    {
        std::cout << "invalid ship" << endl;
        return 1;
    }

    std::vector<int> shipKnots;
    shipKnots.reserve(wayPts.size());

    for (const auto& wp : wayPts)
    {
        const int k = static_cast<int>(std::lround(wp.shipKnots));

        if (k > 0)
            shipKnots.emplace_back(k);
    }

    if (shipKnots.empty())
    {
        std::cout << "No waypoints with positive speed" << endl;
        return 1;
    }

    const auto [minIt, maxIt] = std::minmax_element(shipKnots.begin(), shipKnots.end());
    const int knotsLb = *minIt;
    const int knotsSize = *maxIt - knotsLb + 1;

    const string shipLabel = shipP->getLabel();

    const std::string var = "knots";
    //shipParams.writeFuelCurveCsv(var, waypoints, shipParams, "output.csv");

    //return 1; //sct make it return here

    //NEW CODE START
    // --- replace your single-pass loop with this two-pass helper + calls ---

    auto runCase = [&](const char* title, int nRotor, int nSail) -> double
        {
            // Fresh ship instance so internal precomputed tables build per device type
            shipV shipList2(shipFile, 0);
            shipList2.setNewRowIdx(idx);
            ship* s = shipList2.shipPtr();
            if (!s) { std::cout << "invalid ship\n"; return 0.0; }

            // Select device (IMPORTANT: exactly one of these must be > 0)
            s->setNRotorF(nRotor);
            s->setNSail(nSail);
            s->setRotorScale(0.9);
            s->setSailScale(0.9);

            s->makeTable(knotsLb, knotsSize);

            double sum = 0.0;
            int count = 0;

            std::cout << "\n==== " << title << " run ====\n";

            // Now the full average fuel loop (same as your original, but using 's')
            for (size_t wpIdx = 0; wpIdx < wayPts.size(); ++wpIdx)
            {
                const auto& wp = wayPts[wpIdx];

                const int knots = static_cast<int>(std::round(wp.shipKnots));
                if (knots <= 0) continue;

                const double shipDeg_toConvention = wp.shipBearing;
                const double windDeg_toConvention = wp.windBearing;
                const double waveDeg_toConvention = wp.waveBearing;
                const int    tempCel = static_cast<int>(std::round(wp.temp2m));
                const double windMag = wp.windMag;
                const double swh = wp.swh;
                const double pp1d = wp.pp1d;
                const int gammaIdx = 0; // PM (idx of 1.0 in jsGammaVec)

                const double fcTonsPerNm = s->fcTonsPerNm(
                    knots, tempCel,
                    shipDeg_toConvention, windDeg_toConvention, windMag,
                    waveDeg_toConvention, swh, pp1d, gammaIdx);

                const double fcTonsPerHr = fcTonsPerNm * knots;

                if (fcTonsPerHr > wave::ep) {
                    sum += fcTonsPerHr;
                    count++;

                    const double filterLb{ 1.9 };
                    if (fcTonsPerHr < filterLb)
                        std::cout << "wpIdx " << wpIdx << " filter val < " << filterLb
                        << ", fcTonsPerHr: " << fcTonsPerHr << "\n";
                }
            }

            if (count > 0) {
                const double avgHr = sum / count;
                std::cout << "wayPtCount: " << count
                    << "   Average tonsPerHr: " << avgHr
                    << "   Average tonsPerDay: " << (avgHr * 24.0) << "\n";
                return avgHr;
            }
            else {
                std::cout << "No valid waypoints.\n";
                return 0.0;
            }
        };

    // Run twice: once for rotor, once for hard-sail.
    // Set counts to what you want (e.g., 2 rotors or 2 sails).
    const double base_hr = runCase("Baseline", 0, 0);
    const double rotor_hr = runCase("RotorF", 2, 0);
    const double sail_hr = runCase("HardSail", 0, 2);

    if (base_hr > 0) {
        const double rotor_saving_pct = 100.0 * (base_hr - rotor_hr) / base_hr;
        const double sail_saving_pct = 100.0 * (base_hr - sail_hr) / base_hr;

        std::cout << "\nSummary\n";
        std::cout << "Baseline tonsPerHr: " << base_hr << "\n";
        std::cout << "RotorF  tonsPerHr: " << rotor_hr << "    Saving: " << rotor_saving_pct << " %\n";
        std::cout << "HardSail tonsPerHr: " << sail_hr << "    Saving: " << sail_saving_pct << " %\n";
    }

    std::cout << "\n" << argv[0] << " " << shipFile << " " << wayPtFile << "\n";
    return 0;
}
