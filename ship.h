#ifndef shipH
#define shipH

#include "wave.h"
#include "wind.h"
#include <string>

class fp //fp = force, power
{
public:

    d3 force, power; //dims are [shipDeg][windDeg][knots]

    //sct shorten this class
    fp(const int& _shipDegSize, const int& _windDegSize, const int& _knotsSize)
    {
        wave::makeVec3d(force, _shipDegSize, _windDegSize, _knotsSize, 0.0);
        wave::makeVec3d(power, _shipDegSize, _windDegSize, _knotsSize, 0.0);
    }
};

class ship
{
public:

    ship() {};

    static constinit const int minKnotsAllShip, maxKnotsAllShip, maxRotorF, maxSails, windAid_windMagSize;
    static constinit const double Loa_limU, Tmax_limU, B_limU, displaceVol_limU, sfoc_limU, Abt_limU, ks_limU, Cs_limU, propEff_limU;
    static constinit const double Loa_limL, Tmax_limL, B_limL, displaceVol_limL, sfoc_limL, Abt_limL, propEff_limL;
    static constinit const double waveUpRatioCap, waveDownRatioCap, windUpRatioCap, windDownRatioCap, windAid_windMagStep;

    #ifdef SR_MODEL
        static void printLimit(), printJonswapGamma();
    #endif

    //because we allow user to access ship directly, we make all fields private
    //we friend classes which does peripheral calc so that they can access private fields
    friend class shipV;
    friend class prop;
    friend class wind;
    friend class wave;

    //----
    bool paramOk(), tableEmpty() const;
    int getMinSpeed() const, getMaxSpeed() const, getServiceSpeed() const, getDwt() const, getFuelCap() const;

    double getSFOC() const, getDraft() const, fcDwtScaling(const double& loadFac) const, fcDraftScaling(const double& newDraftDivOld) const;

    //double fcTonsRoute(const std::vector<int>& knots, const i1& tempCel, const d1& shipDeg_toConvention, const d1& windDeg_toConvention,
        //const d1& windMag, const d1& waveDeg_toConvention, const d1& swh, const d1& pp1d, const d1& time, const i1& gammaIdx, const double eff);

    //NEW CODE START
    //int  getNRotorF() const, getNSail() const;
    //void setNRotorF(const int& n), setNSail(const int& n), setRotorScale(const double& s), setSailScale(const double& s);
    //NEW CODE END

    // Compatibility wrapper for callers that expect per hour output
    double fcTonsPerHr(int knots, const int& tempCel, const double& shipDeg_toConvention, const double& windDeg_toConvention, const double& windMag,
                       const double& waveDeg_toConvention, const double& swh, const double& pp1d, const int& gammaIdx, double effIn);
    //NEW CODE END

    double fcTonsPerNm(int knots, const int& tempCel, const double& shipDeg_toConvention, const double& windBearing, const double& windMag,
                       const double& waveDeg_toConvention, const double& swh, const double& pp1d, int gammaIdx = 0, double effIn = -1);

    double fcTonsPerNm_byTable(int knots, const int& tempCel, const double& shipDeg_toConvention, const double& windBearing, const double& windMag,
                               const double& waveDeg_toConvention, const double& swh, const double& pp1d, int gammaIdx = 0, double effIn = -1);

    std::string getLabel() const, getComments() const;

    void makeTable(const int& _knotsLb, const int& _knotsSize), clearTableIfKnotsExceedRange(const int& _knotsLb, const int& _knotsSize);

private:

    //static constexpr double airDensity{1.225};      //air density kg/m3 https://en.wikipedia.org/wiki/Density_of_air
    static double airDensityAtTemp(const double& tempCel); //, getAirAbsViscosity(const double& tempCel);

    static std::pair<int, int> idxOfUniformTable(const double& v, const double& lb, const int& size, const double& step);
    static int idxOfKnotsTable(const int& knots, const int& knotsLb);

    static constexpr double Sapp{0}, SappCoeff{0}, At{0}, lcb{0.5};     // lcb as a percentage of Lpp

    static double cappedWith_Rcalm(const double& upRatio, const double& downRatio, const double& Rcalm, const double& value);

    double bioFoul_Cf(const double& V, const double& Cf, const double& Rn, const int& tempCel) const;
    double bioFoul_townsin(const double& Rn) const;
    double bioFoul_granville(const double& V, const double& Cf, const double& Rn, const int& tempCel) const;

    double getWaveR(const int& knots, const double& waveDegRelShip, const double& swh, const double& pp1d, int gammaIdx) const;
    double getWindR_woAirDensity(const int& knots, const double& windDegRelShip, const double& windMag) const;

    double S_byMumford() const, hm84_Fn_lt04(const int& knots, const int& tempCel) const;
    //----

    double Loa{190};                   // Length overall (m)
    double Lpp{185};                   // Length between perpendiculars (m)
    double B{32};                      // Beam [Width] (m)
    double T{12}, Tmax{12};            // Draft, Maximum Draft (m)
    double displaceVol{53000};         // Displacement (m^3)
    double Abt{12};                    // Transverse bulbous bow area (m^2)
    double sfoc{170};                  // Specific fuel oil consumption (g/kWh)
    double ks{150};                    // Roughness Height of hull (um), val from: What Did We Learn from the Ship Scale Blind CfD Validation Exercise, HullPIC’24
    double Cs{0.26};                   // Cs accounts for different roughness fn
    double rotorScale{1};              // scaling factor for rotoR
    double sailScale{1};               // scaling factor for sails

    double Cwp{0}, Cb{0}, Ca{0}, S{0}, Tf{0}, hb{0}, Fni_fac{0}, Rb_fac{0}, Rw_fac{0}, lmbd{0}, c_15_fac{0}, m_1{0}, k1{0}; //used by hmInit
    double prop_w{0}, prop_t{0}, prop_diameter{0}; //used by prop class

    double R_calm{-1}, R_wave{-1}, R_wind{-1}, P_total{-1}, propEff{-1};

    int minSpeed{10}, maxSpeed{16}, serviceSpeed{minSpeed}, dwt{44000}, fuelCap{150}, nRotorF{0}, nSail{0}, windBallastCurve{0};

    bool hmInitDone{false};
    std::string label{"HandymaxBulker"}, comments{"NoComments"};
    char shipTyChar{'b'};               //'t': tanker, 'g': gas, 'b': bulker, 'a': cargo, 'c': container, 'r': roro

    enum sternTy {s_P = -25, s_V = -10, s_U = 10, s_N = 0};   //numeric values are important, we cast to int later
    sternTy C_stern{s_U};               // s_P:	Very Pronounced, s_V: V-shaped, s_N: Normal, s_U: U-shaped

    enum hullFormTy {h_U = -2, h_V = 2, h_N = 0};   //numeric values are important, we cast to int later
    hullFormTy hullForm{h_U};           // h_V: V-shaped, h_N: Normal, h_U: U-shaped

    void hmInit();

    //fc_byTable
    double getProp_byTable(const int& knots, const double& R) const;
    double getCalmR_byTable(const int& knots, const int& tempCel) const;
    double getWindR_woAirDensity_byTable(const int& knots, const double& windDegRelShip, const double& windMag);
    double getWaveR_byTable(const int& knots, const double& waveDegRelShip, const double& swh, const double& pp1d, int gammaIdx);

    void makeCalmPropTable(const int& _knotsLb, const int& _knotsSize, const int& _tempLb, const int& _tempSize);
    void makeWaveTable(const int& _knotsLb, const int& _knotsSize, const double& _swhLb, const int& _swhSize, const double& _swhStep, const double& _pp1dLb,
        const int& _pp1dSize, const double& _pp1dStep, const double& _waveDegRelShipLb, const int& _waveDegRelShipSize, const double& _waveDegRelShipStep);

    void makeWindTable(const int& _knotsLb, const int& _knotsSize, const double& _windMagLb, const int& _windMagSize, const double& _windMagStep,
        const double& _windDegRelShipLb, const int& _windDegRelShipSize, const double& _windDegRelShipStep);

    void makeTableGivenRange(const int& _knotsLb, const int& _knotsSize, const int& _tempLb, const int& _tempSize, const double& _windMagLb,
        const int& _windMagSize, const double& _windMagStep, const double& _windDegRelShipLb, const int& _windDegRelShipSize, const double& _windDegRelShipStep,
        const double& _swhLb, const int& _swhSize, const double& _swhStep, const double& _pp1dLb, const int& _pp1dSize,
        const double& _pp1dStep, const double& _waveDegRelShipLb, const int& _waveDegRelShipSize, const double& _waveDegRelShipStep);

    //waveTable
    d1 waveR_atAngle;

    int swhSize{0}, pp1dSize{0}, waveKnotsLb{0}, waveKnotsSize{0}, waveDegRelShipSize{0};
    double swhStep{0}, pp1dStep{0}, waveDegRelShipStep{0}, waveDegRelShipLb{0}, swhLb{0}, pp1dLb{0};
    //std::vector<std::vector<std::vector<wave> > > waveTable;
    std::vector<std::vector<wave> > waveTable;

    //windTable
    int windMagSize{0}, windKnotsLb{0}, windKnotsSize{0}, windDegRelShipSize{0};
    double windMagStep{0}, windDegRelShipStep{0}, windDegRelShipLb{0}, windMagLb{0};
    std::vector<std::vector<wind> > windTable;

    //prop, calmTable
    int tempLb{0}, tempSize{0}, calmKnotsLb{0}, calmKnotsSize{0};
    d2 calmTable, propTable_R, propTable_eff;

    //NEW CODE START
    std::vector<fp> windAidVec; //windAidVec[windMag]
    void makeWindAidVec(const int& _knotsLb, const int& _knotsSize);
    std::pair<double, double> getWindAid_byTable(const int& knots, const double& shipDeg_to, const double& windDeg_from, const double& windMag);
    //NEW CODE END
};
#endif
