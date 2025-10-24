#ifndef WIND_H
#define WIND_H

#include <array>
#include <vector>

class ship;

class wind
{
public:

    friend class ship;

    wind() {};
    wind(const ship* _shipP, const int& _knotsLb, const int& _knotsSize, const double& _windMag, const double& _windDegRelShip);

    void init(const ship* _shipP, const int& _knotsLb, const int& _knotsSize, const double& _windMag, const double& _windDegRelShip);
    void makeTable();

    double getR_woAirDensity();

    static constinit const int angleDefaultSize, angleDefaultStep;
    static constinit const double windMag_limU, windMagDefaultStep;

    static double mergeWindWaveBf(const double& windMag, const double& swh);

private:

    int knotsLb{0}, knotsSize{0}, Cx_size{0};
    double V{0}, Vsq{0}, windMagSq{0}, cos_windRelShip_2Mag{0}, Cda_phi_half_Axv{0}, Cda_0_half_Axv{0}, Cx_lb{0}, Cx_step{0};

    std::vector<double> table, Cx;

    double Cda(const double& angleOfAttack) const;

    static double Axv(const double& shipB, const double& Loa, const char& curShip);
    static double interpolBf(const std::vector<double>& windTable, const std::vector<double>& bfTable, const double& wind);

    //----
    static const std::vector<double> kdwt280TankerConvBowLadenCx, lngPrismaticIntegratedCx;
    static const std::vector<double> generalCargoCx, teu6800ContainerLadenCx;
    static const std::vector<double> bfTable, windBfTable, swhBfTable;
    
    //static constexpr std::array<double, 9> Af_a{0.000E+00, 0.000E+00, 1.792E+01, 2.606E+01, 1.132E+00, 1.018E+00, 0.000E+00, -2.765E-02, 5.127E-01};
    //static constexpr std::array<double, 9> Af_b{-3.211E-05, -3.303E-05, 1.140E+00, -2.447E+00, -6.409E-02, -5.437E-02, 8.992E-02, 5.182E-03, 8.065E-02};
    //static constexpr std::array<double, 9> Af_c{2.571E-02, 2.094E-02, -7.515E+01, 2.183E+02, 4.221E+00, 4.203E+00, 8.500E+00, 8.259E-01, -6.399E-01};
    //static constexpr std::array<int, 9> Af_form{4, 4, 1, 1, 3, 3, 3, 6, 3};
};
#endif
