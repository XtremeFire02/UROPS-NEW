#include "csvFn.h"

#include <fstream>
#include <filesystem>
#include <iostream>

using namespace std;

void csvFn::clear()
{
    rowIdx.reset();
    waypoints_.clear();
}

csvFn::csvFn(const std::string& fname)
{
    clear();

    if (!filesystem::exists(fname))
    {
        std::cout << "No Waypt file found." << endl;
    }

    else
    {
        ifstream fIN(fname);
        bool fileError{ false };
        string line;

        getline(fIN, line);

        int lineNum{ 1 };
        while (getline(fIN, line) && !line.empty() && !fileError)
        {
            lineNum++;
            replace(line.begin(), line.end(), ',', ' ');

            // << "Reading line " << lineNum << endl;
            // std::std::cout << endl;
            const auto c = parseLine(line);
            if (c.has_value())
                waypoints_.emplace_back(c.value());
            else
                std::cout << "Invalid parameters for line " << lineNum << endl;
        }

        if (waypoints_.empty()) fileError = true;

        if (fileError)
        {
            std::cout << "Read file failed" << endl;
            clear();
        }
    }
}

const std::vector<Waypoint>& csvFn::getWaypoints() const {
    return waypoints_;
}

int csvFn::getRowCount() const
{
    return static_cast<int>(waypoints_.size());
}

optional<Waypoint> csvFn::parseLine(const string& line)
{

    stringstream stream(line);

    // std::std::cout << "[DEBUG] Raw line = " << line << "\n";

    double _lat{-1}, _lon{-1}, _z{-1}, _forceOcean{-1}, _shipKnots{-1}, _nmTravelled{-1}, _hoursUsed{-1},
        _fuelTonsConsumed{-1}, _fuelTonsPerNM{-1}, _shipBearing{-1},
        _windBearing{-1}, _windMag{-1}, _swh{-1},
        _waveBearing{-1}, _pp1d{-1}, _temp2m{-1};

    stream >> _lat >> _lon >> _z >> _forceOcean >> _shipKnots
        >> _nmTravelled >> _hoursUsed >> _fuelTonsConsumed >> _fuelTonsPerNM
        >> _shipBearing >> _windBearing >> _waveBearing >> _windMag >> _swh >> _pp1d >> _temp2m;

    Waypoint wp;
       
    wp.lat = _lat;
    wp.lon = _lon;
    wp.z = _z;
    wp.forceOcean = _forceOcean;
    wp.shipKnots = _shipKnots;
    wp.nmTravelled = _nmTravelled;
    wp.hoursUsed = _hoursUsed;
    wp.fuelTonsConsumed = _fuelTonsConsumed;
    wp.fuelTonsPerNM = _fuelTonsPerNM;
    wp.shipBearing = _shipBearing;
    wp.windBearing = _windBearing;
    wp.waveBearing = _waveBearing;
    wp.windMag = _windMag;
    wp.swh = _swh;
    wp.pp1d = _pp1d;
    wp.temp2m = _temp2m;

    return make_optional(wp);
}
