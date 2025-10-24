#ifndef WAVE_H
#define WAVE_H

#include <array>
#include <numbers>
#include <vector>
#include <string> 

class ship;

using b1 = std::vector<bool>;
using i1 = std::vector<int>;

using d1 = std::vector<double>;
using d2 = std::vector<d1>;
using d3 = std::vector<d2>;

class wave
{
public:

    friend class ship;

    wave() {};

    wave(const ship* _shipP, const int& _knotsLb, const int& _knotsSize, const double& _swh, const double& _pp1d,
        const double& _waveDegRelShipLb, const int& _waveDegRelShipSize, const double& _windDegRelShipStep);

    void init(const ship* _shipP, const int& _knotsLb, const int& _knotsSize, double _swh, double _pp1d,
        const double& _waveDegRelShipLb, const int& _waveDegRelShipSize, const double& _windDegRelShipStep);

    void makeTable();

    double getR();

    static double knotsToMeterPerSec(const double& V), meterPerSecToKnots(const double& V), celsiusToKelvin(const double& tempCel);
    static double waterDensity(const int& tempCel), waterAbsViscosity(const int& tempCel);

    static double relAngleDegMod180(const double& bearing, const double& shipBearing);
    static double relAngleDegMod180_frConvention(const double& waveDeg_toConvention, const double& shipDeg_toConvention);
    static double make360Deg(double angle), make180Deg(double angle), swapAngleConvention(const double& angle);
    static double degToRad(const double& angle), radToDeg(const double& angle);

    static int nearestIdxOfSortedVec(const std::vector<double>& sortedVec, const double& target);

    static double interpolate2d(const double& x, const double& y, const double& g00, const double& g01, const double& g10, const double& g11,
        const double& x0, const double& x1, const double& y0, const double& y1);

    static double interpolate(const double& x1, const double& x2, const double& y1, const double& y2, const double& x);

    template <typename T>
    static void makeVec1d(std::vector<T>& v, const size_t& sizeS, const T& val)
    {
        if (sizeS == 0)
            return;

        v = std::vector<T>(sizeS, val);
    }

    template <typename T>
    static void makeVec2d(std::vector<std::vector<T> >& v, const size_t& sizeS, const size_t& sizeD, const T& val)
    {
        if (sizeS == 0 || sizeD == 0)
            return;

        std::vector<T> vTmp;
        makeVec1d(vTmp, sizeD, val);

        v = std::vector<std::vector<T> >(sizeS, vTmp);
    }

    template <typename T>
    static void makeVec3d(std::vector<std::vector<std::vector<T> > >& v, const size_t& sizeS, const size_t& sizeD, const size_t& sizeT, const T& val)
    {
        if (sizeS == 0 || sizeD == 0 || sizeT == 0)
            return;

        std::vector<std::vector<T> > vTmp;
        makeVec2d(vTmp, sizeD, sizeT, val);

        v = std::vector<std::vector<std::vector<T> > >(sizeS, vTmp);
    }

    template <typename T>
    static bool isEqual(const T& a, const T& b)
    {
        if (fabs(a - b) < 1e-9)
            return true;
        else
            return false;
    }
    //----

    static constinit const int angleDefaultStep, angleDefaultSize;
    static constinit const double swhDefaultStep, pp1dDefaultStep, swh_limU, pp1d_limL, pp1d_limU;

    static constexpr double smallMag{0.001};
    static constexpr double smallMagSq{smallMag * smallMag};

    static constexpr double noValDouble{1e-6};
    static constexpr double ep{1e-5}; //larger than noValDouble so that we can do if (val < wave::ep)

    static constexpr double gravityConst{9.81}, meterPerSecOfOneKnot{0.5144444444};
    static constexpr double oneThird{1.0 / 3}, twoPi{std::numbers::pi * 2}, halfPi{std::numbers::pi * 0.5};

    static constexpr std::array<int, 30> waterTempCelTable{
        1, 2, 3, 4, 5,
        6, 7, 8, 9, 10,
        11, 12, 13, 14, 15,
        16, 17, 18, 19, 20,
        21, 22, 23, 24, 25,
        26, 27, 28, 29, 30};

    static const d1 jsGammaVec;

    void initSpectrum(const int& gammaIdx);
    void logInitParamsWithR(const int& knotsLb, const int& knotsSize,
        const double& swh, const double& pp1d,
        const double& waveDegRelShipLb, const int& waveDegRelShipSize,
        const double& waveDegRelShipStep, const int& gammaIdx,
        const std::string& filepath = "wave_init_log.csv") const;

    static void enableLogger(const std::string& filepath = "wave_eval_log.csv", bool truncate = true);
    static void disableLogger();
    void startInitLog(const std::string& filepath = "wave_init_log.csv");

private:

    static constexpr double waveNumWaterDepth{600}, Rlp_pitchGyration{0.25}, Ta_minus_Tf{0}; //assume 0 trim, ie Tf - Tf = 0
    static constexpr double Rlp_waterDensity{1025.8}; //at tempCel = 16
    static constexpr double rho_g{wave::Rlp_waterDensity * wave::gravityConst};

    static constexpr std::array<double, 6> Le_regGradTable{-0.7833, -0.4258, -0.4904, -1.061, -0.7414, -0.655};
    static constexpr std::array<double, 6> Le_regInterceptTable{0.8158, 0.5828, 0.5814, 1.049, 0.787, 0.7583};
    static constexpr std::array<double, 6> Lr_regGradTable{-0.6875, -0.8447, -1.04, -0.6722, 1.247, 2.731};
    static constexpr std::array<double, 6> Lr_regInterceptTable{0.7821, 0.8244, 1.081, 0.6952, -0.6726, -1.28};

    static constexpr std::array<double, 9> Rcth_betaTable{0, 30, 45, 60, 90, 120, 135, 150, 180};
    static constexpr std::array<double, 9> Rcth_betaCorrectionTable{1, 0.925, 0.9, 0.8, 0.75, 0.7, 0.7, 0.7, 0.6};

    //ITTC2011, Fresh Water and Seawater Properties, Table 3
    static constexpr std::array<double, 30> waterDensityTable{
        1028.0941, 1028.0197, 1027.9327, 1027.8336, 1027.7225,
        1027.6000, 1027.4662, 1027.3214, 1027.1659, 1027.0000,
        1026.8238, 1026.6376, 1026.4416, 1026.2360, 1026.0210,
        1025.7967, 1025.5633, 1025.3210, 1025.0700, 1024.8103,
        1024.5421, 1024.2656, 1024.9808, 1024.6881, 1023.3873,
        1023.0788, 1022.7626, 1022.4389, 1022.1078, 1021.7694};

    static constexpr std::array<double, 30> waterAbsViscosityTable{
        1.843e-06, 1.783e-06, 1.726e-06, 1.671e-06, 1.620e-06,
        1.571e-06, 1.524e-06, 1.480e-06, 1.438e-06, 1.397e-06,
        1.359e-06, 1.322e-06, 1.286e-06, 1.252e-06, 1.220e-06,
        1.189e-06, 1.159e-06, 1.131e-06, 1.103e-06, 1.077e-06,
        1.051e-06, 1.027e-06, 1.004e-06, 9.81e-07, 9.59e-07,
        9.38e-07, 9.18e-07, 8.98e-07, 8.79e-07, 8.61e-07};

    static const double ceil_4_pitchGyration, floor_4_pitchGyration;
    static const d2 R_aw_a, R_aw_b, R_aw_c;

    //static constinit const double waveNumTol;
    //static constinit const int waveNumMaxIter;

    enum shipTy {tanker = 0, gas = 1, bulker = 2, cargo = 3, container = 4, roro = 5};

    enum dir {head = 0, beam = 1, follow = 2}; //must be 0, 1, 2 cos we cast to int to access the vector
    enum R_aw_choice {LP_n_CTH = 0, LP = 1, CTH = 2};
    // combined:A meta-model for added resistance in waves by Kim et al. (2022)
    // LP:      Regression analysis of experimental data for added resistance in waves of arbitrary heading and development of a semi-empirical formula by Liu & Papanikolaou (2020)
    // CTH:     A Practical Speed Loss Prediction Model at Arbitrary Wave Heading for Ship Voyage Optimization by Lang & Mao (2021)

    //static void calcFreqSpectrum_PM(const int& nBins, d1& _sVec, d1& _wVec, const double& A, const double& B, const double& b);
    //static void calcFreqSpectrum_JS(const int& nBins, d1& _sVec, d1& _wVec, const double& jsGamma, const double& b);
    //void calcFreqSpectrum(const int& gammaIdx);

    void calcFreqSpectrum_sct(const int& gammaIdx);
    void calcDirSpectrum(const double& waveDeg);

    double R_aw(const double& ab, const double& ww, const R_aw_choice& R_aw_method) const;
    double R_aw_combined(const double& ab, const double& ww) const;

    double Rlp_aw(const double& ab, const double& ww) const;
    double Rlp_a1_following(const double& U, const double& Vg) const;
    double Rlp_a1(const double& aa, const double& _U, const double& lambda) const;
    double Rlp_a1_head(const double& aa, const double& Fn) const;

    double Rcth_aw(const double& ab, const double& ww) const;
    double Rcth_aw_original(const double& ab, const double& ww) const;
    double Rcth_corrected_wBar(const double& beta, const double& wBar) const;

    //double integral2d();
    //----
    shipTy curShip{bulker};

    bool swhZero{false};
    int knotsLb{0}, knotsSize{0}, waveDegRelShipSize{0};

    double shipB{0}, shipCb, shipLpp{0}, shipT{0}, shipTmax{0}; //copied from ship cos we may parallel this class
    double waveDegRelShipLb{0}, waveDegRelShipStep{0}, knots{0}; //use double for knots
    double Le{0}, Lr{0};
    double waveDegRelShipFrConvention{0};

    double zetaSq{0}, rho_g_zeta2_shipB2_div_Lpp{0}, rho_g_zeta2_shipB{0}, pow_Lpp_div_B{0}, sqrt_Lpp_div_g{0};
    double Rlp_a1Fac, Rlp_a1_f0{0}, Rlp_a1_f1{0}, Rlp_wBarFac{0}, spectrum_A{0}, spectrum_B{0}, spectrum_b{0}, spectrum_m0{0};

    d1 sVec, gVec, wVec, alphaVec;
    d3 table;

    static double getS_of_JS(const double& x, const double& m0_PM, const double& lambda);
    static double getS_of_PM(const double& x, const double& m0);

    static int waterTempIdx(const int& tempCel);

    static double froudeNum(const double& U, const double& Lpp);
    static double waveNum(const double& ww), Rlp_lambda(const double& ww);
    static double integral1d(const d1& _sVec, const d1& _wVec), g(const double& aa, const double& bb, const double& vv);

    static void logRow(const wave* self, const char* tag, double R_LP, double R_CTH, double R_TOTAL);
    double integral2d_with(R_aw_choice choice) const;
};
#endif
