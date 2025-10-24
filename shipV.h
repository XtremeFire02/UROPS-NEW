#ifndef SHIP_V_H
#define SHIP_V_H

#include "csvFn.h"
#include "ship.h"

#include <optional>
#include <string>
#include <vector>

class shipV
{
public:

    shipV(const std::string& fname, const int& rIdx);

    static const std::vector<std::string> defaultShips;

    static void printShipCsv();
    static std::vector<std::vector<std::string> > readCsv(const std::string& fname);

    static double aveAngle(const std::vector<double>& angles, const std::vector<double>& weights), ave(const std::vector<double>& values, const std::vector<double>& weights);

    void writeFuelCurveCsv(const std::string& paramName, const std::vector<Waypoint>& waypoints, shipV& shipList, const std::string& filename);
    void setNewRowIdx(const int& rIdx);

    bool isLabelUnique() const;
    int getRowCount() const, getRowIdx() const;

    double getSFOC(), getFcDwtScaling(const double& newLoadDivOld), getFcDraftScaling(const double& newDraftDivOld);
    std::tuple<double, double, double, double, double, double> fcTonsPerNm(int knots, const int& tempCel, const double& shipDeg_toConvention, const double& windDeg_toConvention, const double& windMag,
        const double& waveDeg_toConvention, const double& swh, const double& pp1d, int gammaIdx = 0, double effIn = -1);

    ship* shipPtr();

    std::string getShipLabel();

private:

    std::optional<int> rowIdx;
    std::vector<ship> vec;

    static const std::vector<std::string> header;
    static const std::string csvExt, delim, llmToolGuideStr;

    static double aveAngleNormalise(const double& angle);
    static std::optional<ship> readShipLine(const std::string& line);

    void clear();
    void readShipCsv(const std::string& fname);    
};
#endif
