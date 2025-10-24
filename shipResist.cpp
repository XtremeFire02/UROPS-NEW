#include "csvFn.h"
#include "shipV.h"
#include "wave.h"

#include <cmath>
#include <filesystem>
#include <iostream>

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
    string shipFile = "shipHandy.csv"; //sct integration
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

    shipV shipList(shipFile, 1); //sct integration
    //const int nRows = shipList.getRowCount();

    //const int idx{0}; //use row0 of file
    //shipList.setNewRowIdx(idx);

    ship* shipP = shipList.shipPtr();

    if (shipP == nullptr)
    {
        std::cout << "invalid ship" << endl;
        return 1;
    }

    const string shipLabel = shipP->getLabel(); //sct integration

    const std::string var = "knots";
    //shipParams.writeFuelCurveCsv(var, waypoints, shipParams, "output.csv");
    

    {
        int knots{16};
        int temp{26};

        double shipB{0};
        double windB{90};
        double windM{5.56612};
        double waveB{90};
        double swh{4};
        double pp1d{1};

        vector<int> shipKnots;

        for (const auto& it : wayPts)
        {
            const int k = static_cast<int>(std::lround(it.shipKnots));

            if (k > 0)
            {
                //cout << "k " << k << endl;
                shipKnots.emplace_back(k);
            }
        }

        const int knotsLb = *min_element(shipKnots.begin(), shipKnots.end());
        const int knotsSize = *max_element(shipKnots.begin(), shipKnots.end()) - knotsLb + 1;

        std::cout << "byTable knotsLb: " << knotsLb << " knotsSize: " << knotsSize << endl;
        shipP->makeTable(knotsLb, knotsSize);

        const int gammaIdx = 0;

        cout << "shipLabel " << shipP->getLabel() << endl; //sct integration

        const auto fcTonsPerNm = shipP->fcTonsPerNm(knots, temp, shipB, windB, windM, waveB, swh, pp1d, gammaIdx);
        std::cout << "fcTonsPerHr " << fcTonsPerNm * knots << endl;

        std::cout << "----" << endl;

        const auto fcTonsPerNmVec = shipP->fcTonsPerNm_byTable(knots, temp, shipB, windB, windM, waveB, swh, pp1d, gammaIdx);
        std::cout << "fcTonsPerHrByVec " << fcTonsPerNmVec * knots << endl;

        std::cout << "----" << endl;
    }
    
    return 1; //sct make it return here


    double sum{0};
    int count{0};
    size_t wpMax = 1; //use just 1 wayPt

    for (int wpIdx = 0; wpIdx < wayPts.size(); wpIdx++)
    {
        const auto& wp = wayPts[wpIdx];

        const int knots = static_cast<int>(round(wp.shipKnots));

        const double shipDeg_toConvention = wp.shipBearing;
        const double windDeg_toConvention = wp.windBearing;
        const double waveDeg_toConvention = wp.waveBearing;

        const int tempCel = static_cast<int>(round(wp.temp2m));

        const double windMag = wp.windMag;

        const double swh = wp.swh;
        const double pp1d = wp.pp1d;
        const double speed = wp.shipKnots;

        if (knots <= 0)
            continue;

        const int gammaIdx = 0; //wave::jsGammaVec{1.0, 3.3}, so 0 = PM, 1 = gamma3.3

        //cout << "pp1d: " << pp1d << endl;

        //all bearing in waypoint uses to_convention, fcTonsPerNm will take in to_convention
        //std::cout << "shipBearing: " << shipDeg_toConvention << " windBearing: " << windDeg_toConvention << " waveBearing: " << waveDeg_toConvention << endl;

        const double fcTonsPerNm = shipP->fcTonsPerNm(knots, tempCel, shipDeg_toConvention, windDeg_toConvention, windMag, waveDeg_toConvention, swh, pp1d, gammaIdx);

        const double fcTonsPerHr = fcTonsPerNm * knots;

        if (fcTonsPerHr > wave::ep)
        {
            std::cout << speed << endl;
            std::cout << fcTonsPerHr * 24 << endl;
            sum += fcTonsPerHr;
            count++;

            const double filterLb{1.9};

            if (fcTonsPerHr < filterLb)
                std::cout << "wpIdx " << wpIdx << " filter val < " << filterLb << ", fcTonsPerHr: " << fcTonsPerHr << endl;
        }
    }

    std::cout << endl << "wayPtCount: " << count << " Average tonsPerHr: " << sum / count << " Average tonsPerDay: " << sum / count * 24 << endl;
    std::cout << endl << argv[0] << " " << shipFile << " " << wayPtFile << endl;

    return 0;
}
