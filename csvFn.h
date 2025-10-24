#ifndef CSV_FN_H
#define CSV_FN_H

#include <string>
#include <vector>
#include <optional>

using s1 = std::vector<std::string>;
using s2 = std::vector<s1>;

struct Waypoint
{
    double lat;
    double lon;
    double z;
    bool forceOcean;
    double shipKnots;
    double nmTravelled;
    double hoursUsed;
    double fuelTonsConsumed;
    double fuelTonsPerNM;
    double shipBearing;
    double windBearing;
    double waveBearing;
    double windMag;
    double swh;
    double pp1d;
    double temp2m;
};

class csvFn
{
public:
    
    csvFn(const std::string& filename);

    const std::vector<Waypoint>& getWaypoints() const;
    int getRowCount() const;

    std::vector<Waypoint> waypoints_;
    std::optional<int> rowIdx;

private:

    void clear();
    std::optional<Waypoint> parseLine(const std::string& line);

    // Parse header to build column index mapping
    bool parseHeader(const std::string& headerLine);
    static std::vector<std::string> splitCSV(const std::string& line);

    // Column indices for required fields; -1 if not present
    std::vector<int> colIdx; // ordered as: lat,lon,z,forceOcean,shipKnots,nmTravelled,hoursUsed,fuelTonsConsumed,fuelTonsPerNM,shipBearing,windBearing,waveBearing,windMag,swh,pp1d,temp2m
    bool hasHeader{false};
};
#endif
