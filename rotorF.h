#ifndef ROTOR_F_H
#define ROTOR_F_H

#include "wave.h"
#include "wind.h"

#include <iostream>

class sail
{
public:

    static constexpr double rhoAir{1.225};
    static constexpr double halfRhoAir{0.5 * sail::rhoAir}; //sct

    sail()
    {
        CL_s = CL * areaPerSail;  //sct
        CD_s = CD * areaPerSail;
    }

    // Minimal rigid sail model for table precomputation
    double getForce(const double& windSpeed_knots, const int& shipSpeed_knots,
        const double& shipDeg_to, const double& windDeg_from, const int& numSail) const
    {
        const double Vship{wave::knotsToMeterPerSec(shipSpeed_knots)};
        const double Vtrue{wave::knotsToMeterPerSec(windSpeed_knots)};

        const double rel_from{sail::relAngleDeg(windDeg_from, shipDeg_to)};
        const double rel_to{wave::swapAngleConvention(rel_from)};
        const double beta{wave::make180Deg(rel_to)};

        //const double alpha{(180.0 - beta) * std::numbers::pi / 180.0;

        const double alpha{wave::degToRad(180.0 - beta)};

        //const double Vapp = std::sqrt(Vtrue * Vtrue + Vship * Vship
        //- 2.0 * Vtrue * Vship * std::cos(std::numbers::pi - alpha));

        //sct
        const double Vapp = std::sqrt(Vtrue * Vtrue + Vship * Vship
                                      + 2.0 * Vtrue * Vship * std::cos(alpha)); //sct, cos (180 - a) = -a

        const double T = sail::halfRhoAir * Vapp * Vapp * (CL_s * std::sin(alpha) - CD_s * std::cos(alpha));

        return std::max(0.0, T * numSail);
    }

    static double relAngleDeg(const double& from_deg, const double& to_deg)
    {
        // "from" and "to" are bearings in degrees (0..360). Return relative angle [0,360).
        double x{std::fmod(from_deg - to_deg, 360.0)};

        if (x < 0)
            x += 360.0;

        return x;
    }

    private:

    double CL{1.20}, CD{0.12}, areaPerSail{1500}; // m^2 per sail
    double CL_s{0}, CD_s{0};
};

//Flettner rotorF has the following effects:
//increase air resistance
//consumes power
//produces force

class rotorF
{
public:

    rotorF(const int& r)
    {
        if (r >= 1 && r <= 4)
        {
            numRotor = r;
        }
        else
        {

            std::cout << r << " numRotor must be in [1, 4]; using 1" << std::endl;
            numRotor = 1;
        }

        const double A_side{std::numbers::pi * diameter * height}; //sct

        //const double Re{0.5 * diameterEff};
        //const double A_plates{wave::twoPi * Re * Re};

        const double A_plates{wave::halfPi * diameterEff * diameterEff};

        A_wetted = A_side + A_plates;

        // use the same constants as power
        const double S = height * diameter;

        // simple CL and CD model
        const double CL = 1.1 * SR;
        constexpr double CD{0.10};

        CL_S = CL * S;
        CD_S = CD * S;
    }

    static constexpr int shipDegSize{74};   // 0, 5, 10, ..., 180 sct make it 0..360
    static constexpr double shipDegStep{5.0};

    static constexpr int windDegSize{74};   // 0, 5, 10, ..., 180
    static constexpr double windDegStep{5.0};

    // Convert true wind + ship motion into apparent wind magnitude (m/s)
    // and apparent wind FROM bearing (deg, 0..360, nautical convention).
    static std::pair<double, double> apparent_from_true(const double& Vtrue_mps, const double& Vship_mps,
        const double& shipToDeg, const double& windFromDeg)
    {
        const double th_ship = wave::degToRad(shipToDeg);
        const double th_from = wave::degToRad(windFromDeg);

        // True wind vector points toward where it is going, opposite of “from”.
        const double wx = -Vtrue_mps * std::sin(th_from);
        const double wy = -Vtrue_mps * std::cos(th_from);

        // Ship velocity vector
        const double sx = Vship_mps * std::sin(th_ship);
        const double sy = Vship_mps * std::cos(th_ship);

        // Apparent wind (felt on ship) = true wind - ship velocity
        const double ax = wx - sx;
        const double ay = wy - sy;

        const double Vapp = std::sqrt(ax * ax + ay * ay);

        // “from” bearing is opposite the vector direction
        //const double bx{-ax}, by{-ay};
        //double bear_deg = wave::radToDeg(std::atan2(bx, by)); // atan2(x,y) for bearing

        double bear_deg = wave::radToDeg(std::atan2(-ax, -ay)); // atan2(x,y) for bearing

        if (bear_deg < 0)
            bear_deg += 360.0;

        return {Vapp, bear_deg};
    }

    //Eq9 of Wind-Assisted Ship Propulsion: Matching Flettner Rotors with Diesel Engines and Controllable Pitch Propellers
    // First argument is apparent wind speed in m s^?1
    // The two bearing arguments are intentionally unused
    double getPowerConsumption(const double& windMag_mps)
    {
        if (numRotor <= 0)
            return 0.0;

        const double Vapp{std::max(0.0, windMag_mps)};
        const double U_tan{SR * Vapp};                   // rotorF surface speed

        // frictional/drive power per rotorF, corrected by efficiency
        const double P_one = Cf * sail::halfRhoAir * U_tan * U_tan * U_tan * A_wetted / etaMR;
        return P_one * numRotor;
    }

    //Given the number of rotorF, we use a look up table approach, see Table 6.2.3.2 of the below
    //China Classification Society Guidelines for Survey of Marine Wind-rotorF Assisted Propulsion System 2023
    //Figures can be extracted from
    //Design, operation and analysis of wind-assisted cargo ships Table4, 5
    //Note that TWA, true wind angle, head wind = 0 -> from_Convention, but counter-clockwise, so to convert need 360 - angle
    //Flettner rotorF has diff dimension, this paper is 5m by 30m. Below has same diameter, so we can scale by height

    //https://www.offshore-energy.biz/dealfeng-rotorF-sails-for-hung-zes-14000-dwt-tankers/
    //As informed, each vessel will feature a 5m x 24m Dealfeng rotorF sail installed on its forecastle deck.

    //Also see Figure 4 of Wind-Assisted Ship Propulsion: Matching Flettner Rotors with Diesel Engines and Controllable Pitch Propellers
    double getForce(const double& windMag, const double& shipBearing, const double& windBearing)
    {
        if (numRotor <= 0)
            return 0.0;

        double V{windMag};

        if (V <= 0)
            V = wave::ep;
        else if (V > wind::windMag_limU)
            V = wind::windMag_limU;

        //relative angle FROM wind to ship TO, 0..360
        //double rel = std::fmod(windBearing - shipBearing, 360.0);

        //if (rel < 0)
            //rel += 360.0;

        // reduce to 0..180 and convert to inflow angle about ship axis
        //double beta = (rel <= 180.0) ? rel : (360.0 - rel);

        double beta{wave::relAngleDegMod180(windBearing, shipBearing)}; //sct
        double alpha{wave::degToRad(180.0 - beta)};

        const double q{sail::halfRhoAir * V * V};
        const double T_one = q * (CL_S * std::sin(alpha) - CD_S * std::cos(alpha));

        // do not return negative assist
        return std::max(0.0, T_one) * numRotor;
    }

private:

    int numRotor{0};
    // --- constants for the power model ---

    const double SR{3.0};       // spin ratio (matches your thrust model)
    const double Cf{0.008};     // friction/drive coefficient (~few x 1e-3)
    const double etaMR{0.60};      //mech/electrical efficiency, must not be 0 else div by 0

    // Rotor geometry (close to your comments: 5m x 24~30m + end plates)
    const double diameter = 5.0;       // diameter (m)
    const double height = 24.0; // height scales with 'scale'
    const double diameterEff = 5.0;       // end-plate effective diameter (m)

    double A_wetted{0}, CL_S{0}, CD_S{0};
};
#endif
