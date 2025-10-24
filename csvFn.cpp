#include "csvFn.h"

#include <fstream>
#include <filesystem>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <unordered_map>

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

        // Read header and build column index mapping
        if (!getline(fIN, line))
        {
            std::cout << "Empty Waypt file." << endl;
            return;
        }
        hasHeader = parseHeader(line);

        int lineNum{ 1 };
        while (getline(fIN, line) && !fileError)
        {
            lineNum++;
            if (line.empty()) continue;

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

// Split CSV on commas (no quote handling required for our inputs)
vector<string> csvFn::splitCSV(const string& l)
{
    vector<string> out;
    string cur;
    cur.reserve(l.size());
    for (char ch : l)
    {
        if (ch == ',')
        {
            out.emplace_back(cur);
            cur.clear();
        }
        else
            cur.push_back(ch);
    }
    out.emplace_back(cur);
    return out;
}

static inline string trim(const string& s)
{
    size_t b = 0, e = s.size();
    while (b < e && isspace(static_cast<unsigned char>(s[b]))) b++;
    while (e > b && isspace(static_cast<unsigned char>(s[e - 1]))) e--;
    return s.substr(b, e - b);
}

bool csvFn::parseHeader(const string& headerLine)
{
    const auto cols = splitCSV(headerLine);
    unordered_map<string, int> idx;
    for (int i = 0; i < static_cast<int>(cols.size()); i++)
    {
        idx[trim(cols[i])] = i;
    }

    // Required columns
    const vector<string> names{
        "lat","lon","z","forceOcean","shipKnots",
        // optional: arrivalUTC (ignored if present)
        "nmTravelled","hoursUsed","fuelTonsConsumed","fuelTonsPerNM",
        "shipBearing_toConvention","windBearing_toConvention","waveBearing_toConvention",
        "windMag","swh","pp1d","2m_Temp"};

    // Build ordered index list for fast parse
    colIdx.clear();
    colIdx.reserve(16);

    // First five are always present
    for (int i = 0; i < 5; i++)
    {
        if (!idx.count(names[i])) return false;
        colIdx.push_back(idx[names[i]]);
    }

    // Handle optional arrivalUTC between shipKnots and nmTravelled
    // We do not add it to colIdx; simply skip by mapping the rest by name
    for (int i = 5; i < static_cast<int>(names.size()); i++)
    {
        if (!idx.count(names[i])) return false;
        colIdx.push_back(idx[names[i]]);
    }

    return true;
}

optional<Waypoint> csvFn::parseLine(const string& line)
{
    if (hasHeader && colIdx.size() == 16)
    {
        const auto t = splitCSV(line);
        auto getd = [&](int ci) -> optional<double>
        {
            if (ci < 0 || ci >= static_cast<int>(t.size())) return nullopt;
            const string s = trim(t[ci]);
            if (s.empty()) return nullopt;
            try { return stod(s); } catch (...) { return nullopt; }
        };

        Waypoint wp{};
        const vector<optional<double>> v{
            getd(colIdx[0]),  getd(colIdx[1]),  getd(colIdx[2]),  getd(colIdx[3]),
            getd(colIdx[4]),  getd(colIdx[5]),  getd(colIdx[6]),  getd(colIdx[7]),
            getd(colIdx[8]),  getd(colIdx[9]),  getd(colIdx[10]), getd(colIdx[11]),
            getd(colIdx[12]), getd(colIdx[13]), getd(colIdx[14]), getd(colIdx[15]) };

        // Require the first 5 and bearings to be present; others can default to 0
        for (int req : {0,1,2,3,4,9,10,11})
            if (!v[req].has_value()) return nullopt;

        wp.lat = v[0].value();
        wp.lon = v[1].value();
        wp.z = v[2].value();
        wp.forceOcean = v[3].value();
        wp.shipKnots = v[4].value();
        wp.nmTravelled = v[5].value_or(0.0);
        wp.hoursUsed = v[6].value_or(0.0);
        wp.fuelTonsConsumed = v[7].value_or(0.0);
        wp.fuelTonsPerNM = v[8].value_or(0.0);
        wp.shipBearing = v[9].value();
        wp.windBearing = v[10].value();
        wp.waveBearing = v[11].value();
        wp.windMag = v[12].value_or(0.0);
        wp.swh = v[13].value_or(0.0);
        wp.pp1d = v[14].value_or(0.0);
        wp.temp2m = v[15].value_or(0.0);

        return make_optional(wp);
    }

    // Fallback: old behavior (for legacy files without header)
    string lineSp = line;
    replace(lineSp.begin(), lineSp.end(), ',', ' ');
    stringstream stream(lineSp);

    double _lat{-1}, _lon{-1}, _z{-1}, _forceOcean{-1}, _shipKnots{-1}, _nmTravelled{-1}, _hoursUsed{-1},
        _fuelTonsConsumed{-1}, _fuelTonsPerNM{-1}, _shipBearing{-1},
        _windBearing{-1}, _windMag{-1}, _swh{-1},
        _waveBearing{-1}, _pp1d{-1}, _temp2m{-1};

    // Expect comma-less values (caller previously replaced commas)
    stream >> _lat >> _lon >> _z >> _forceOcean >> _shipKnots
        >> _nmTravelled >> _hoursUsed >> _fuelTonsConsumed >> _fuelTonsPerNM
        >> _shipBearing >> _windBearing >> _waveBearing >> _windMag >> _swh >> _pp1d >> _temp2m;

    if (!stream)
        return nullopt;

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
