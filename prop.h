#ifndef PROP_H
#define PROP_H

#include "ship.h"

class prop
{
public:

    prop(ship* shipPtr, const int& knots, const int& tempCel);

    static double makeFrac(const double& v);
    static bool isFrac(const double& v);

    double getEff(const double& R) const;

    static constexpr double shaftEff{0.98};    //assumed value as found in pg 16 of [21]
    static constexpr double rotationEff{1};    //assumption in pg16 of [21]

    static constinit const double Rfrac;
    static constinit const int Rsize, assumedTemp;

    static double diameter(const char& shipTyChar, const double& shipTmax);
    static double wakeFrac(ship* shipP, const double& Cb, const double& ldr, const double& Dprop), thrustDeductionFac(ship* shipP, const double& Cb, const double& ldr, const double& Dprop);

private:

    double Cth_fac{0}, hullEff{0};
    double Bseries(const double& R) const, Bseries_gFn(const double& R) const;
};
#endif
