#include <cmath>
#include <iostream>
#include <optional>

#include "wave.h"

using namespace std;

// Wave Model by Lang & Mao
// Main Reference: [18] Lang, X., & Mao, W. (2021). A practical speed loss prediction model at arbitrary wave heading for ship voyage optimization. Journal of Marine Science and Application, 20(3), 410–425. https://doi.org/10.1007/s11804-021-00224-z
// Literature uses from_convention : 0 indicates head waves

/*
double wave::Rcth_aw(const double& aa, const double& ww) const
{
    const double PI = std::numbers::pi;
    const double HPI = std::numbers::pi * 0.5;

    // 90–180 stays unchanged
    if (aa >= HPI - wave::ep) {
        return Rcth_aw_original(aa, ww);
    }
    const double cth_beam = Rcth_aw_original(HPI, ww);
    const double cth_head = Rcth_aw_original(PI, ww);

    const double rlp_tail = Rlp_aw(0.0, ww); //Rlp must be computed before this fn
    const double rlp_beam = Rlp_aw(HPI, ww);
    const double rlp_head = Rlp_aw(PI, ww);

    //std::cout << std::fixed
    //    << "[RCTH dbg] aa=" << aa << " rad, ww=" << ww << " rad/s\n"
    //    << "  cth_beam=" << cth_beam << ", cth_head=" << cth_head << '\n'
    //    << "  rlp_tail=" << rlp_tail << ", rlp_beam=" << rlp_beam << ", rlp_head=" << rlp_head
    //    << std::endl;

    if (!(std::isfinite(cth_beam) && std::isfinite(cth_head) &&
        std::isfinite(rlp_tail) && std::isfinite(rlp_beam) && std::isfinite(rlp_head))) {
        return Rcth_aw_original(aa, ww);
    }

    double cth_tail_1;
    if (std::abs(rlp_head) > wave::ep) cth_tail_1 = cth_head * (rlp_tail / rlp_head);
    else cth_tail_1 = std::nan("");

    double cth_tail_2;
    if (std::abs(rlp_beam) > wave::ep) cth_tail_2 = cth_beam * (rlp_tail / rlp_beam);
    else cth_tail_2 = std::nan("");

    const bool t1ok = std::isfinite(cth_tail_1);
    const bool t2ok = std::isfinite(cth_tail_2);
    if (!t1ok && !t2ok) {
        return Rcth_aw_original(aa, ww);
    }
    const double cth_tail = (t1ok && t2ok) ? 2 * (cth_tail_1 * cth_tail_2) / (cth_tail_1 + cth_tail_2)
        : (t1ok ? cth_tail_1 : cth_tail_2);
    if (!std::isfinite(cth_tail)) {
        return Rcth_aw_original(aa, ww);
    }

    double t = aa / HPI;                       
    if (t < 0.0) t = 0.0; else if (t > 1.0) t = 1.0;
    return cth_tail + t * (cth_beam - cth_tail);
}
*/

//Brain's code: std::isfinite == false for NaN, and mostly true otherwise, so it is only useful for NaN (basically it does not do much)
//using cth_orig for some cases & not for others will introduce discontinuity
//we use meanFn(cth_orig, cth_byRlpHead, cth_byRlpBeam) and to drop cth_orig ONLY if it is out of range
//this reduces discontinuity, cos Rlp is 'always' there
//for meanFn, use AM, shifted GM (make all pos), HM -> div by 0 or small value problem (we do not use it)

//given that rlp is always fine, we use rlpMax = max(rlp_tail, rlp_beam, rlp_head) to determine whether cth_tail is out of range
//syntax-wise, we use optional<double> to drop cth_orig when it is out of range
double wave::Rcth_aw(const double& aa, const double& ww) const
{
    // 90–180 stays unchanged
    if (aa >= wave::halfPi - wave::ep)
        return Rcth_aw_original(aa, ww);

    const double cth_beam = Rcth_aw_original(wave::halfPi, ww);
    const double cth_head = Rcth_aw_original(numbers::pi, ww);

    const double rlp_tail = Rlp_aw(0.0, ww);     //Rlp must be computed before this fn
    const double rlp_beam = Rlp_aw(wave::halfPi, ww);
    const double rlp_head = Rlp_aw(numbers::pi, ww);

    optional<double> cth_tail = Rcth_aw_original(aa, ww);

    constexpr double fac{3}; //3 may be too loose, hv other way to tighten it?

    const double rlpMax = max(max(fabs(rlp_tail), fabs(rlp_beam)), fabs(rlp_head));

    if (fabs(cth_tail.value()) > rlpMax * fac || fabs(cth_tail.value()) < rlpMax / fac)
        cth_tail.reset(); //make cth_tail no value

    const double rlpTailDivHead = rlp_tail / (wave::ep + rlp_head); //ep avoids NaN
    const double rlpTailDivBeam = rlp_tail / (wave::ep + rlp_beam);

    optional<double> cth_tail_byRlpHead, cth_tail_byRlpBeam;

    if (fabs(rlpTailDivHead) < 10)
        cth_tail_byRlpHead = cth_head * rlpTailDivHead; //when rlp_head = 0, cth_tail_byRlpHead -> no value

    if (fabs(rlpTailDivBeam) < 10)
        cth_tail_byRlpBeam = cth_beam * rlpTailDivBeam;

    const vector<optional<double> > tailValues{cth_tail, cth_tail_byRlpHead, cth_tail_byRlpBeam};

    double mostNegVal{0}; //shift to pos then use GM, must init with 0

    for (const auto& it : tailValues)
    {
        if (it.has_value() && it.value() < 0) //note < 0 is needed
        {
            if (mostNegVal > it.value())
                mostNegVal = it.value();
        }
    }

    int nValueInMean{0};

    double sum{0};
    double product{1};

    for (const auto& it : tailValues)
    {
        if (it.has_value())
        {
            sum += it.value();
            product *= (it.value() - mostNegVal); //minus cos mostNegVal is always a neg val
            nValueInMean++;
        }
    }

    double cth_tail_byMean{rlp_tail}; //rlp_tail is the backup value, when nValueInMean == 0

    if (nValueInMean > 0)
    {
        cth_tail_byMean = sum / nValueInMean;
        //cth_tail_byMean = pow(product, 1.0 / nValidVal) + mostNegVal; //sign is reversed from prev
    }

    double t = aa / wave::halfPi;

    if (t < 0.0)
        t = 0.0;
    else if (t > 1.0)
        t = 1.0;

    return cth_tail_byMean + t * (cth_beam - cth_tail_byMean);
}

double wave::Rcth_corrected_wBar(const double& beta_deg, const double& wBar) const
{
    // Correction in accordance with Eq 14 (Lang & Mao, 2021)
    //static_assert(wave::RCTHBetaTable.size() == wave::Rcth_betaCorrFacTable.size(), "Rcth_betaCorrFacTable size wrong");

    //warning
    if (beta_deg < 0.0 || beta_deg > 180.0)
    {
        cout << "Error: invalid beta range" << endl;
        return wave::noValDouble;
    }

    for (size_t i = 0; i < Rcth_betaTable.size() - 1; i++)
    {
        const double betaLow{Rcth_betaTable[i]};
        const double betaHigh{Rcth_betaTable[i + 1]};

        if (wave::isEqual(beta_deg, betaLow)) return wBar * Rcth_betaCorrectionTable[i];
        else if (wave::isEqual(beta_deg, betaHigh)) return wBar * Rcth_betaCorrectionTable[i + 1];
        else if (beta_deg > betaLow && beta_deg < betaHigh)
        {
            const double fLow{ Rcth_betaCorrectionTable[i] };
            const double fHigh{ Rcth_betaCorrectionTable[i + 1] };

            return wBar * interpolate(betaLow, betaHigh, fLow, fHigh, beta_deg);
        }
    }

    cout << "Error in Rcth_corrected_wBar: beta out of range";
    return wave::noValDouble;
}

// All equations from [18]
double wave::Rcth_aw_original(const double& aa, const double& ww) const
{
    const double k{wave::waveNum(ww)};
    const double lambda{wave::Rlp_lambda(ww)};
    const double U{wave::knotsToMeterPerSec(knots)};
    const double Fn{wave::froudeNum(U, shipLpp)};

    const double aa_deg{radToDeg(aa)};
    const double beta_deg{wave::make180Deg(wave::swapAngleConvention(aa_deg))};
    const double beta{degToRad(beta_deg)};
    const double cosBeta{cos(beta)};

    // R_awr
    const double E{atan(shipB / (2 * Le))};
    const double Bf{2.25 * pow(sin(E), 2)};                // Eq 3
    const double omega{ww * U / wave::gravityConst};
    const double ke{k * pow(1 + omega * cosBeta, 2)};

    const double Cu{max(-310 * Bf + 68, 10.0)};           // Eq 5

    //const double alpha_T{1 - exp(-2 * k_e * shipT)};       // Eq 4
    const double alphaU_alphaT{(1.0 + Cu * Fn) * (1 - exp(-2 * ke * shipT))}; //(1 + alphaU) * alphaT

    const double R_awr_0{0.5
                         * rho_g_zeta2_shipB
                         * Bf
                         * alphaU_alphaT
                         * (0.19 / shipCb)
                         * pow(lambda / shipLpp, Fn - 1.11)}; // Eq 2

    // R_awm
    const double a_1{60.3 * pow(shipCb, 1.34) * pow(shipCb, -(1 + Fn))}; // Eq 7

    const double a_2 = Fn < 0.12 ?
                       0.0072 + 0.24 * Fn :
                       pow(Fn, -1.05 * shipCb + 2.3) * exp(Fn * (-2 - wave::ceil_4_pitchGyration - wave::floor_4_pitchGyration)); // Eq 8

    const double inv_c_1{1 / (0.4567 * shipCb / wave::Rlp_pitchGyration + 1.689)}; // Eq 11

    constexpr double pow05_143{0.6515574435275924}; //pow(0.05, 0.143)

    const double wBar = Fn < 0.05 ?
                        ww * (sqrt_Lpp_div_g * pow(wave::Rlp_pitchGyration, inv_c_1) * pow05_143) / (1.09 + wave::ceil_4_pitchGyration * 0.08) :
                        ww * (sqrt_Lpp_div_g * pow(wave::Rlp_pitchGyration, inv_c_1) * pow(Fn, 0.143)) / (1.09 + wave::ceil_4_pitchGyration * 0.08); // Eq 11

    const double wBarCorrected{wave::Rcth_corrected_wBar(beta_deg, wBar)}; // Eq 14

    const double b_1 =
                    (wBarCorrected < 1 && shipCb < 0.75) ? (19.77 * shipCb / wave::Rlp_pitchGyration - 36.39) / wave::ceil_4_pitchGyration :
                    (wBarCorrected < 1 && shipCb >= 0.75) ? 11.0 / wave::ceil_4_pitchGyration :
                    (wBarCorrected >= 1 && shipCb < 0.75) ? -12.5 / wave::ceil_4_pitchGyration :
                    -5.5 / wave::ceil_4_pitchGyration;                      // Eq 9

    const double d_1 =
                    (wBarCorrected < 1 && shipCb < 0.75) ? 14 :
                    (wBarCorrected < 1 && shipCb >= 0.75) ? 566 * pow_Lpp_div_B * 2 :
                    -566 * pow_Lpp_div_B * 6;                               // Eq 10

    const double R_awm_0{4
                         * rho_g_zeta2_shipB2_div_Lpp
                         * pow(wBarCorrected, b_1)
                         * exp(b_1 / d_1 * (1 - pow(wBarCorrected, d_1)))
                         * a_1 * a_2};                                      // Eq 6

    // Added angle
    const double R_awr = beta <= wave::halfPi ?
                        R_awr_0 * pow(Fn, Fn * (floor(cosBeta) - ceil(cosBeta))) * cosBeta :
                        R_awr_0 * pow(Fn, -1.5 * Fn * (floor(cosBeta) + ceil(cosBeta))) * cosBeta; // Eq 13

    const double R_awm{R_awm_0
                        * exp(-pow(shipB / numbers::pi, 4 * sqrt(Fn)))
                        + rho_g_zeta2_shipB2_div_Lpp
                        * pow(lambda / shipB * max(cosBeta, 0.45), -6 * Fn)
                        * sin(beta)};                                       // Eq 15, note + rho..., whereas Eq6 is * rho...
    
    /*
    std::cout << "CTH dbg: beta_deg=" << beta_deg
        << " cosBeta=" << cosBeta
        << " R_awr0=" << R_awr_0
        << " R_awr=" << R_awr
        << " R_awm=" << R_awm 
        << " R_CTH" << R_awm + R_awr
        << "\n";
    */ 
    //constexpr double kNegFrac = 0.75;           // allow R_awr to be at most 75% of R_awm in magnitude
    //const double R_awr_clamped = std::max(R_awr, -kNegFrac * R_awm);
    //return R_awm + R_awr_clamped;

    return R_awm + R_awr;
}
