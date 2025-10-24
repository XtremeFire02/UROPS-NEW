#include <cmath>
#include <iostream>

#include "wave.h"

using namespace std;

// Wave Model by Liu & Papa
// Main Reference: [16] Liu, S., & Papanikolaou, A. (2020). Regression analysis of experimental data for added resistance in waves of arbitrary heading and development of a semi-empirical formula. Ocean Engineering, 206, 107357. https://doi.org/10.1016/j.oceaneng.2020.107357
// Literature uses to_convention: 180 indicates head waves, same as convention in program

// Eq 16 of [16], case when alpha = 0
// Methodology outlined in 5.3. Added resistance due to radiation effect in following waves (Page 8 - 10) of [16]
double wave::Rlp_a1_following(const double& U, const double& Vg) const
{
    const double s{U / Vg};       // dimensionless speed

    constexpr double fHalf{0.0};   // Speed = Vg/2 -> Resistance = 0

    // Speed = 0 -> Froude Number = Angle = 0, substitute Eq 14 of [16]
    //const double f0 = 60.3
    //                  * pow(shipC_B, 1.34)
    //                  * pow(4 * wave::Rlp_pitchGyration, 2)
    //                  * pow(0.87 / shipC_B, -1)
    //                  / log(shipB / shipTmax)
    //                  * -1 * oneThird;

    // Speed = Vg -> Froude Number = 0, Angle = 180, substitute Eq 14 of [16]
    //const double f1 = 60.3
    //                  * pow(shipC_B, 1.34)
    //                  * pow(4 * wave::Rlp_pitchGyration, 2)
    //                  * 0.87 / shipC_B
    //                  / log(shipB / shipTmax);

    if (s < 0)
    {
        //warning
        cout << "Error: Invalid U in Rlp_calc_a1_following" << endl;
        return wave::noValDouble;
    }
    else if (s <= 0.5)
    {
        return interpolate(0.0, 0.5, Rlp_a1_f0, fHalf, s);
    }
    else if (s <= 1.0)
    {
        return interpolate(0.5, 1.0, fHalf, Rlp_a1_f1, s);
    }
    else
    {
        const double Fn{wave::froudeNum(U - Vg, shipLpp)};
        return Rlp_a1_head(numbers::pi, Fn); // Speed = Vg -> Speed = U - Vg, angle = 180, substitute Eq 14 of [16]
    }
}

// Eq 14 of [16]
double wave::Rlp_a1_head(const double& aa, const double& Fn) const
{
    const double cos_aa{cos(aa)};

    return Rlp_a1Fac * pow(0.87 / shipCb, -1 * cos_aa * (1 + Fn)) * (1 - (2 * cos_aa)) * wave::oneThird;
}

// Eq 16 of [16]
double wave::Rlp_a1(const double& aa, const double& U, const double& lambda) const
{
    const double Vc{sqrt(wave::gravityConst * lambda / wave::twoPi)}; // Vc = sqrt(g * lambda / (2 * pi)) , Pg 9 of [16]

    const double Vg{0.5 * Vc};    // Eq 8 of [16]

    const double Fn{wave::froudeNum(U, shipLpp)};
    
    const double a1_follow{Rlp_a1_following(U, Vg)};

    if (aa >= wave::halfPi)         // when pi/2 <= alpha <= pi
    {
        return Rlp_a1_head(aa, Fn);
    }
    else if (fabs(aa) <= wave::ep) // when alpha == 0
    {
        return a1_follow;
    }    
    else if (aa > wave::ep)        // when 0 <= alpha <= pi/2
    {
        const double a1_beam = Rlp_a1_head(wave::halfPi, Fn);
        return interpolate(0.0, wave::halfPi, a1_follow, a1_beam, aa); //interpolate aa = 0 and aa = pi/2
    }
    else
    {
        //warning
        cout << "Error: aa input value invalid" << endl;
        return wave::noValDouble;
    }
}

// Eq 16 & 17 of [16]
double wave::Rlp_aw(const double& aa, const double& ww) const    //aa is in radians
{
    const double U{wave::knotsToMeterPerSec(knots)};
    const double Fn{wave::froudeNum(U, shipLpp)};
    const double alpha{aa};

    const double lambda{wave::Rlp_lambda(ww)};

    if (aa < 0.0)
    {
        cout << "aa < 0" << endl;
        return wave::noValDouble;
    }

    if (aa > numbers::pi)
    {
        cout << "aa > pi" << endl;
        return wave::noValDouble;
    }

    const double cosAlpha{cos(alpha)};

    // R_awm (Eq 16 of [16])
    const double a_1{wave::Rlp_a1(alpha, U, lambda)};

    constexpr double oneFourteenth{1.0 / 14};

    //const double w_bar{2.142
    //                   * pow(wave::Rlp_pitchGyration, wave::oneThird)
    //                   * sqrt(shipLpp / lambda)
    //                   * (1 - 0.111 / shipC_B * (log(shipB / shipTmax) - log(2.75)))
    //                   * pow(shipC_B / 0.65, 0.17)
    //                   * ((-1.377 * pow(Fn, 2) + 1.157 * Fn) * abs(cos(alpha)) + (0.618 * (13 + cos(2 * alpha)) * oneFourteenth))};

    const double wBar{Rlp_wBarFac * sqrt(shipLpp / lambda)
                    * ((-1.377 * pow(Fn, 2) + 1.157 * Fn) * fabs(cosAlpha) + (0.618 * (13 + cos(2 * alpha)) * oneFourteenth))};

    const double a_2 = Fn < 0.12 ?
                       0.0072 + 0.1676 * Fn :
                       pow(Fn, 1.5) * exp(-3.5 * Fn);

    const double a_3{1.0 + 28.7 * atan(wave::Ta_minus_Tf / shipLpp)};

    const double b_1 = wBar < 1 ? 11.0 : -8.5;

    const double d_1 = wBar < 1 ?
                        566 * pow(shipLpp * shipCb / shipB, -2.66) :
                        -566 * pow(shipLpp / shipB, -2.66) * (4 - 125 * atan(wave::Ta_minus_Tf) / shipLpp);

    const double R_awm{4
                        * rho_g_zeta2_shipB2_div_Lpp
                        * a_1 * a_2 * a_3
                        * pow(wBar, b_1)
                        * exp(b_1 / d_1 * (1 - pow(wBar, d_1)))};

    // R_awr (Eq 17 of [16])
    const double Ecommon{0.495 * shipB};

    const double E_1{atan(Ecommon / Le)};
    const double E_2{atan(Ecommon / Lr)};

    const double fAlpha = alpha >= numbers::pi - E_1 ? -cosAlpha : 0;

    const double Tstar_12{shipTmax};

    const double Tstar_34 = shipCb <= 0.75 ?
                            shipTmax * (4 + sqrt(fabs(cosAlpha))) * 0.2 :
                            shipTmax * (2 + sqrt(fabs(cosAlpha))) * wave::oneThird;

    const double a_Tstar_12 = lambda / shipLpp <= 2.5 ?
        1 - exp(-4 * numbers::pi * (Tstar_12 / lambda - Tstar_12 / (2.5 * shipLpp))) : 0;

    const double a_Tstar_34 = lambda / shipLpp <= 2.5 ?
        1 - exp(-4 * numbers::pi * (Tstar_34 / lambda - Tstar_34 / (2.5 * shipLpp))) : 0;

    /*
    const double w_0{ww}; // Incident wave frequency

    const double R_awr_1{ alpha > E_1
        ? 2.25 * 0.25 * wave::Rlp_waterDensity * wave::gravityConst * shipB * pow(zetaA, 2) * a_Tstar_12 * (pow(sin(E_1 - alpha),2) + 2 * w_0 * U / wave::gravityConst * (cos(E_1) * cos(E_1 - alpha) - cos(alpha)))
        * pow(0.87 / shipC_B, (1 + 4 * pow(Fn, 0.5)) * f_alpha)
        : 0 };
    const double R_awr_2{alpha > numbers::pi - E_1
        ? 2.25 * 0.25 * wave::Rlp_waterDensity * wave::gravityConst * shipB * pow(zetaA, 2) * a_Tstar_12 * (pow(sin(E_1 + alpha),2) + 2 * w_0 * U / wave::gravityConst * (cos(E_1) * cos(E_1 + alpha) - cos(alpha)))
        * pow(0.87 / shipC_B, (1 + 4 * pow(Fn, 0.5)) * f_alpha)
        : 0 };
    const double R_awr_3{alpha < numbers::pi - E_2
        ? -2.25 * 0.25 * wave::Rlp_waterDensity * wave::gravityConst * shipB * pow(zetaA, 2) * a_Tstar_34 * (pow(sin(E_2 + alpha),2) + 2 * w_0 * U / wave::gravityConst * (cos(E_2) * cos(E_2 + alpha) - cos(alpha)))
        : 0 };
    const double R_awr_4{alpha < E_2
        ? -2.25 * 0.25 * wave::Rlp_waterDensity * wave::gravityConst * shipB * pow(zetaA, 2) * a_Tstar_34 * (pow(sin(E_2 - alpha),2) + 2 * w_0 * U / wave::gravityConst * (cos(E_2) * cos(E_2 - alpha) - cos(alpha)))
        : 0 };
    */

    const double two_w0_U_div_g{2 * ww * U / wave::gravityConst}; // Incident wave frequency

    const double cosE1{cos(E_1)};
    const double cosE2{cos(E_2)};

    const double fFac{pow(0.87 / shipCb, (1 + 4 * sqrt(Fn)) * fAlpha)};

    const double R_awr_1 = alpha > E_1 ?
         2.25 * 0.25 * rho_g_zeta2_shipB * a_Tstar_12 * (pow(sin(E_1 - alpha), 2) + two_w0_U_div_g * (cosE1 * cos(E_1 - alpha) - cosAlpha)) * fFac : 0;

    const double R_awr_2 = alpha > numbers::pi - E_1 ?
         2.25 * 0.25 * rho_g_zeta2_shipB * a_Tstar_12 * (pow(sin(E_1 + alpha), 2) + two_w0_U_div_g * (cosE1 * cos(E_1 + alpha) - cosAlpha)) * fFac : 0;

    const double R_awr_3 = alpha < numbers::pi - E_2 ?
        -2.25 * 0.25 * rho_g_zeta2_shipB * a_Tstar_34 * (pow(sin(E_2 + alpha), 2) + two_w0_U_div_g * (cosE2 * cos(E_2 + alpha) - cosAlpha)) : 0;

    const double R_awr_4 = alpha < E_2 ?
        -2.25 * 0.25 * rho_g_zeta2_shipB * a_Tstar_34 * (pow(sin(E_2 - alpha), 2) + two_w0_U_div_g * (cosE2 * cos(E_2 - alpha) - cosAlpha)) : 0;

    const double R_awr = R_awr_1 + R_awr_2 + R_awr_3 + R_awr_4;

    return R_awr + R_awm;
}

// Note that k = 2 * pi / lambda, Eq 22 of [17]
double wave::Rlp_lambda(const double& ww)
{
    return wave::twoPi / wave::waveNum(ww); //potential div by 0. since we have bound w away from 0, not a issue
}
