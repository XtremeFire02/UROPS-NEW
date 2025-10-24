#include <cmath>
#include <iostream>

#include "ship.h"
#include "wave.h"

using namespace std;

// Calm Water Resistance Implementation: Holtrop and Mennon Method
// Reference: [1] Holtrop, J., & Mennen, G. (1982). An approximate power prediction method. International Shipbuilding Progress, 29(335), 166–170. https://doi.org/10.3233/isp-1982-2933501
double ship::hm84_Fn_lt04(const int& knots, const int& tempCel) const
{
    const double V{wave::knotsToMeterPerSec(knots)};
    const double V2{V * V};

    const double rho{wave::waterDensity(tempCel)};
    const double halfRhoV2{0.5 * rho * V2};

    const double Rn{V * Lpp / wave::waterAbsViscosity(tempCel)};     // Reynolds Number
    const double Fn{wave::froudeNum(V, Lpp)};

    // R_F 
    // Reference: [5] ITTC Resistance Test 7.5-02-02-01, 2.1: Data Reduction Equations
    const double Cf_tmp{0.075 / pow(log10(Rn) - 2, 2) * k1};           // [5]: ITTC 57 Model-Ship Correlation Line
    const double Cf{Cf_tmp + bioFoul_Cf(V, Cf_tmp, Rn, tempCel)};

    const double Rf{halfRhoV2 * S * Cf};                          // [5]: Total Resistance Coefficient

    //most water is deep, we skip Shallow Water Correction of Reference: [8] ITTC Preparation, Conduct and Analysis of Speed/Power Trials

    const double Rapp{halfRhoV2 * ship::Sapp * ship::SappCoeff * Cf};

    // R_W
    constexpr double d{-0.9};

    const double m_4{c_15_fac * exp(-0.034 * pow(Fn, -3.29))};

    // Handle Fn < 0.4 and 0.4 < Fn < 0.55 separately
    /*
    double R_W;
    if (Fn <= 0.4) R_W = c_1 * c_2 * c_5 * displaceVol * rho * wave::gravityConst * exp(m_1 * pow(Fn, d) + m_4 * cos(lmbd * pow(Fn, -2)));
    else if (Fn > 0.4 && Fn <= 0.55)
    {
        const double m_4_Fn4 = c_15 * 0.4 * exp(-0.034 * pow(0.4, -3.29));
        const double R_W_Fn4 = c_1 * c_2 * c_5 * displaceVol * rho * wave::gravityConst * exp(m_1 * pow(Fn, d) + m_4_Fn4 * cos(lmbd * pow(Fn, -2)));

        const double m_4_Fn55 = c_15 * 0.4 * exp(-0.034 * pow(0.55, -3.29));
        const double R_W_Fn55 = c_1 * c_2 * c_5 * displaceVol * rho * wave::gravityConst * exp(m_1 * pow(Fn, d) + m_4_Fn55 * cos(lmbd * pow(Fn, -2)));

        R_W = R_W_Fn4 + ((10 * Fn - 4) * (R_W_Fn55 - R_W_Fn4) / 1.5);
    }
    else
    {
        cout << "Error calculating R_calm: Fn > 0.55";
        return wave::noValDouble;
    }*/

    const double cosLmbd_div_Fn2{cos(lmbd / (Fn * Fn))};

    double Rw;

    if (Fn <= 0.4)
    {
        Rw = Rw_fac * rho * exp(m_1 * pow(Fn, d) + m_4 * cosLmbd_div_Fn2);
    }
    else if (Fn > 0.4 && Fn <= 0.55)
    {
        const double m_4_Fn4{c_15_fac * exp(-0.034 * pow(0.4, -3.29))};
        const double Rw_Fn4{Rw_fac * rho * exp(m_1 * pow(Fn, d) + m_4_Fn4 * cosLmbd_div_Fn2)};

        const double m_4_Fn55{c_15_fac * exp(-0.034 * pow(0.55, -3.29))};
        const double Rw_Fn55{Rw_fac * rho * exp(m_1 * pow(Fn, d) + m_4_Fn55 * cosLmbd_div_Fn2)};

        Rw = Rw_Fn4 + ((10 * Fn - 4) * (Rw_Fn55 - Rw_Fn4) / 1.5);
    }
    else
    {
        cout << "Error calculating R_calm: Fn > 0.55";
        return wave::noValDouble;
    }

    // R_B
    const double Fni{V / sqrt(wave::gravityConst * Fni_fac + 0.15 * V2)};
    const double Fni_sq{Fni * Fni};

    const double Rb{Rb_fac * Fni_sq * Fni * rho / (1 + Fni_sq)}; // #bulbous bow resistance [N]

    // R_TR
    double Rtr;
    if (ship::At > 0)
    {
        const double FnT{V / sqrt(2 * wave::gravityConst * ship::At / (B + B * Cwp))};
        const double c_6 = FnT < 5 ? 0.2 * (1 - 0.2 * FnT) : 0;

        Rtr = halfRhoV2 * ship::At * c_6; // #transom resistance [N]
    }
    else
        Rtr = 0;

    // R_A
    const double Ra{halfRhoV2 * S * Ca};

    //cout << "calm Rf " << Rf << endl;             // Frictional Resistance
    //cout << "calm Rw " << Rw << endl;             // Wave-making Resistance
    //cout << "calm Rb " << Rb << endl;             // Bulbous Bow Resistance
    //cout << "calm Ra " << Ra << endl;             // Air Resistance

    //these are 0
    //cout << "calm Rapp " << Rapp << endl;         // Appended Resistance
    //cout << "calm Rtr " << Rtr << endl;           // Transom Resistance
    
    // R_total
    return Rf + Rapp + Rw + Rb + Rtr + Ra; // #[N]
}

//part of hm84_Fn_lt04 shifted here, to speed up calmTable
void ship::hmInit()
{
    if (hmInitDone) //do only once
        return;

    S = S_byMumford();

    // Reference: [2] Papanikolaou A. 2014. Ship design: methodologies of preliminary design. Dordrecht: Springer.
    // Section 2.9: Selection of Hull Form Coefficients, Page 140-142
    if (hullForm == h_U)
        Cwp = 0.778 * Cb + 0.248;     // Eq 2.81
    else if (hullForm == h_N)
        Cwp = (1 + 2 * Cb) / 3.0;     // Eq 2.83
    else if (hullForm == h_V)
        Cwp = 0.743 * Cb + 0.297;     // Eq 2.85
    else
    {
        cout << "C_WP: unknown hull form type" << endl;
    }

    // Reference: [3] Schneekluth H, Bertram V. 1998. Ship design for efficiency and economy (Vol. 218). Oxford: Butterworth- Heinemann.
    // Section 1.5: Midship section area coefficient and midship section design, Table 2
    const double Cm{1.006 - 0.0056 * pow(Cb, -3.56)};
    const double Cp{min(0.99, displaceVol / (Lpp * Cm * B * T))}; //must not exceed 1 cos ship is a block
    //----

    const double n_pow_OneThird{pow(displaceVol, wave::oneThird)};
    const double L3_div_n{Lpp * Lpp * Lpp / displaceVol};

    const double BT{B * T};
    const double B_div_L{B / Lpp};
    const double L_div_B{Lpp / B};

    const double sqrt_Abt{sqrt(Abt)};
    const double Pb{0.56 * sqrt_Abt / (Tf - 1.5 * hb)};
    const double Abt_pow_3div2{pow(Abt, 1.5)};

    Rb_fac = 0.11 * Abt_pow_3div2 * exp(-3 / (Pb * Pb)) * wave::gravityConst;

    const double Lr{Lpp * (1 - Cp + 0.06 * Cp * ship::lcb / (4 * Cp - 1.0))};
    const double c_14{1 + 0.011 * static_cast<int>(C_stern)};
    k1 = 0.93 + 0.487118 * c_14 * pow(B_div_L, 1.06806) * pow(T / Lpp, 0.46106) * pow(Lpp / Lr, 0.121563) * pow(L3_div_n, 0.36486) * pow(1 - Cp, -0.604247);

    lmbd = L_div_B <= 12 ? 1.446 * Cp - 0.03 * L_div_B : 1.446 * Cp - 0.36;

    double c_15;

    if (L3_div_n <= 512)
        c_15 = -1.69385;
    else if (512 < L3_div_n && L3_div_n < 1726.91)
        c_15 = -1.69385 + (Lpp / n_pow_OneThird - 8) / 2.36;
    else
        c_15 = 0;

    c_15_fac = c_15 * 0.4;

    const double c_16 = Cp < 0.8 ?
                        8.07981 * Cp - 13.8673 * pow(Cp, 2) + 6.984388 * pow(Cp, 3) :
                        1.73014 - 0.7067 * Cp;

    double c_7;

    if (B_div_L <= 0.11)
        c_7 = 0.229577 * pow(B_div_L, wave::oneThird);
    else if (0.11 < B_div_L && B_div_L < 0.25)
        c_7 = B_div_L;
    else
        c_7 = 0.5 - 0.0625 * L_div_B;

    const double i_E{1 + 89 * exp(-pow(L_div_B, 0.80856) * pow(1 - Cwp, 0.30484) * pow(1 - Cp - 0.0225 * ship::lcb, 0.6367)
                                    * pow(Lr / B, 0.34574) * pow(100 / L3_div_n, 0.16302))};

    const double c_1{2223105 * pow(c_7, 3.78613) * pow(T / B, 1.07961) * pow(90 - i_E, -1.37565)};
    const double c_3{0.56 * Abt_pow_3div2 / (BT * (0.31 * sqrt_Abt + Tf - hb))};
    const double c_2{exp(-1.89 * sqrt(c_3))};
    const double c_5{1 - 0.8 * ship::At / (BT * Cm)};

    Rw_fac = c_1 * c_2 * c_5 * displaceVol * wave::gravityConst;

    m_1 = 0.0140407 * Lpp / T - 1.75254 * n_pow_OneThird / Lpp - 4.79323 * B_div_L - c_16;

    const double Tf_div_L{Tf / Lpp};
    const double c_4 = Tf_div_L <= 0.04 ? Tf_div_L : 0.04;

    Ca = 0.006 * pow(Lpp + 100, -0.16) - 0.00205 + 0.003 * sqrt(Lpp / 7.5) * pow(Cb, 4) * c_2 * (0.04 - c_4);

    hmInitDone = true;
}

// Reference: [4] Kristensen, H.O, Lützen, M. "Prediction of Resistance and Propulsion Power of Ships", May 2013.
// Page 4
double ship::S_byMumford() const
{
    enum mumfordShipTy {ave = 0, bulkerTanker = 1, container = 2, twinScrew = 3, twinSkeg = 4, doubleFerry = 5};
    mumfordShipTy mumfordShip{ave};

    switch (shipTyChar)
    {
        case 't': mumfordShip = bulkerTanker; break;
        case 'b': mumfordShip = bulkerTanker; break;
        case 'c': mumfordShip = container; break;
        default:
            mumfordShip = ave;      // no info, assumed average
    }

    constexpr array<double, 6> outer{1.025, 0.99, 0.995, 1.53, 1.2, 1.11};
    constexpr array<double, 6> inner{1.700, 1.90, 1.900, 0.55, 1.5, 1.70};

    const int idx{static_cast<int>(mumfordShip)};

    return outer[idx] * (displaceVol / T + inner[idx] * Lpp * T);
}
