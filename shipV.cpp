#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>

#include <filesystem>

#include "wave.h"
#include "shipV.h"
#include "prop.h"
#include "rotorF.h"

#ifdef SR_MODEL
    #include "crypt.h"
#endif

using namespace std;

//sct, keep this short, smaller than 4096 model context
const string shipV::llmToolGuideStr{ "{'title': 'ship-fuel-rate tool guide', 'swh limit in meters': {'Panamax bulker': 3, 'Panamax container': 4}}" };

const vector<string> shipV::header{"ShipLabel", "Comments", "ShipTy(t/g/b/a/c/r)", "LenOverall(m)", "LenPerpendiculars(m)", "Beam(m)", "Draft(m)", "MaxDraft(m)", "DisplacementVol(m^3)", "BulbousBowTransverseArea(m^2)",
    "SFOC(g/kwh)", "ks(micrometer;hullRoughness)", "SternTy(p/v/u/n)", "HullForm(v/u/n)", "DWT(tons;int)", "FuelCap(tons;int;neg=unused)", "MinSpeed(kts;int)", "MaxSpeed(kts;int)",
    "ServiceSpeed(kts;int)", "NumRotorF(0..4;int)", "RotorScale(0..1)", "NumSail(0..2;int)", "SailScale(0..1)"};

const vector<string> shipV::defaultShips{
    "HandymaxBulker1,NoComments,b,190,185,32,12,12,53000,12,170,150,u,u,44000,150,10,16,14,0,1,0,1",
    "HandymaxBulker2,NoComments,b,190,185,32,12,12,53000,12,170,150,u,u,44000,150,10,16,14,0,1,0,1"};

const string shipV::csvExt{"csv"};
const string shipV::delim{","};

shipV::shipV(const string& fname, const int& rIdx)
{
    if (wave::jsGammaVec.empty())
    {
        cout << "wave::jsGammaVec.empty()" << endl;
        return;
    }

    //warning byTable
    if (prop::Rfrac > 0.3 - wave::ep)
        cout << "prop::Rfrac too large" << endl;

    if (prop::Rsize < 15)
        cout << "prop::Rsize too small" << endl;

    //----

    if (wind::windMag_limU < 19.9)
        cout << "wind::windMag_limU too small" << endl;

    if (wind::windMagDefaultStep > 3.1)
        cout << "wind::windMagDefaultStep too large" << endl;

    //----
    if (wave::swh_limU < 5.9)
        cout << "wave::swh_limU  too small" << endl;

    if (wave::pp1d_limU < 19.9)
        cout << "wave::pp1d_limU too small" << endl;

    if (wave::pp1d_limL > 3.1)
        cout << "wave::pp1d_limL too large" << endl;

    //----
    if (wave::angleDefaultStep > 45)
        cout << "wave::angleDefaultStep too large" << endl;

    if ((wave::angleDefaultSize - 1) * wave::angleDefaultStep != 180)
        cout << "last wave angle != 180" << endl;

    if (wind::angleDefaultStep > 45)
        cout << "wind::angleDefaultStep too large" << endl;

    if ((wind::angleDefaultSize - 1) * wind::angleDefaultStep != 180)
        cout << "last wind angle != 180" << endl;


    if (rotorF::shipDegSize * rotorF::shipDegStep < 359)
        cout << "rotorF::shipDegSize * rotorF::shipDegStep < 359" << endl;

    if (rotorF::windDegSize * rotorF::windDegStep < 359)
        cout << "rotorF::windDegSize * rotorF::windDegStep < 359" << endl;

    if (ship::windAid_windMagSize < 15)
        cout << "ship::windAid_windMagSize < 15" << endl;


    //----

    readShipCsv(fname);

    if (!vec.empty() && rIdx >= 0 && rIdx < static_cast<int>(vec.size()))
    {
        rowIdx = rIdx;
    }
    else
    {
        clear();
        cout << "Ship parameter parsing failed." << endl;
    }
}

ship* shipV::shipPtr()
{
    if (!vec.empty() && rowIdx.has_value())
        return &vec[rowIdx.value()];
    else
        return nullptr;
}

tuple<double, double, double, double, double, double> shipV::fcTonsPerNm(int knots, const int& tempCel, const double& shipDeg_toConvention, const double& windDeg_toConvention, const double& windMag,
    const double& waveDeg_toConvention, const double& swh, const double& pp1d, int gammaIdx, const double effIn)
{
    #ifdef SR_MODEL
        crypt::exitIfNoLic(true);
    #endif

    if (shipPtr() != nullptr)
    {
        const double fc = shipPtr()->fcTonsPerNm(knots, tempCel, shipDeg_toConvention, windDeg_toConvention, windMag, waveDeg_toConvention, swh, pp1d, gammaIdx, effIn);
        return {fc, shipPtr()->P_total, shipPtr()->propEff, shipPtr()->R_calm, shipPtr()->R_wind, shipPtr()->R_wave};
    }
    else
        return {-1, -1, -1, -1, -1, -1};
}

double shipV::getSFOC()
{
    if (shipPtr() != nullptr)
        return shipPtr()->getSFOC();
    else
        return -1;
}

double shipV::getFcDwtScaling(const double& newLoadDivOld)
{
    if (shipPtr() != nullptr)
        return shipPtr()->fcDwtScaling(newLoadDivOld);
    else
        return -1;
}

double shipV::getFcDraftScaling(const double& newDraftDivOld)
{
    if (shipPtr() != nullptr)
        return shipPtr()->fcDwtScaling(newDraftDivOld);
    else
        return -1;
}

string shipV::getShipLabel()
{
    if (shipPtr() != nullptr)
        return shipPtr()->getLabel();
    else
        return "";
}

void shipV::setNewRowIdx(const int& rIdx)
{
    if (rIdx >= 0 && rIdx < static_cast<int>(vec.size()))
        rowIdx = rIdx;
}

int shipV::getRowIdx() const
{
    if (!vec.empty() && rowIdx.has_value())
        return rowIdx.value();
    else
        return -1;
}

void shipV::clear()
{
    rowIdx.reset();
    vec.clear();
}

void shipV::writeFuelCurveCsv(const string& paramName, const vector<Waypoint>& waypoints, shipV& shipList, const string& filename)
{
    ofstream fout(filename);
    fout << paramName << ",FuelTonsPerHr" << endl;

    for (const auto& wp : waypoints)
    {
        const int knots = static_cast<int>(wp.shipKnots);
        const double shipBearing = wp.shipBearing;
        const int tempCel = static_cast<int>(wp.temp2m);
        const double windBearing = wp.windBearing;
        const double windMag = wp.windMag;
        const double waveBearing = wp.waveBearing;
        const double swh = wp.swh;
        const double pp1d = wp.pp1d;

        ship* s = shipList.shipPtr();

        if (s != nullptr)
        {
            int gammaIdx = 0;
            const double fcTonsPerNm = s->fcTonsPerNm(knots, tempCel, shipBearing, windBearing, windMag, waveBearing, swh, pp1d, gammaIdx);

            if (fcTonsPerNm > 0 && knots > 0)
            {
                const double fcTonsPerHr = fcTonsPerNm * knots;
                double paramVal = 0.0;

                if (paramName == "swh") paramVal = swh;
                else if (paramName == "windMag") paramVal = windMag;
                else if (paramName == "knots") paramVal = knots;
                else if (paramName == "pp1d") paramVal = pp1d;
                else if (paramName == "temp") paramVal = tempCel;
                else if (paramName == "waveBearing") paramVal = waveBearing;
                else
                {
                    cout << "Unsupported parameter: " << paramName << endl;
                    return;
                }

                fout << fixed << setprecision(3) << paramVal << "," << fcTonsPerHr << endl;
            }
        }
    }

    fout.close();
}

void shipV::printShipCsv()
{
    const string outFile = string("ship.") + shipV::csvExt;
    ofstream fout(outFile);

    if (!fout)
    {
        cout << "Failed to open " << outFile << " for writing.\n";
        return;
    }

    // Print header
    for (size_t i = 0; i < shipV::header.size(); ++i)
    {
        fout << shipV::header[i] << (i + 1 < shipV::header.size() ? shipV::delim : "\n");
    }

    ship c;
    constexpr int nRow{2};

    for (int i = 1; i <= nRow; i++)
    {
        fout << c.label + to_string(i) << shipV::delim
            << c.comments << shipV::delim
            << c.shipTyChar << shipV::delim
            << c.Loa << shipV::delim
            << c.Lpp << shipV::delim
            << c.B << shipV::delim
            << c.T << shipV::delim
            << c.Tmax << shipV::delim
            << c.displaceVol << shipV::delim
            << c.Abt << shipV::delim
            << c.sfoc << shipV::delim
            << c.ks << shipV::delim
            //<< c.Cs << shipV::delim
            << 'u' << shipV::delim     //must be char not enum
            << 'u' << shipV::delim
            //<< c.windBallastCurve << shipV::delim
            << c.dwt << shipV::delim
            << c.fuelCap << shipV::delim
            << c.minSpeed << shipV::delim
            << c.maxSpeed << shipV::delim
            << c.serviceSpeed << shipV::delim
            << c.nRotorF << shipV::delim
            << c.rotorScale << shipV::delim
            << c.nSail << shipV::delim
            << c.sailScale << endl;
    }
}

optional<ship> shipV::readShipLine(const string& line)
{
    stringstream stream(line);

    double _L_OA{-1}, _Lpp{-1}, _B{-1}, _T{-1}, _Tmax{-1},
        _n{-1}, _A_BT{-1}, _sfoc{-1}, _ks{-1}, _rotorScale{-1}, _sailsScale{-1};

    int _dwt{-1}, _fuelCap{-1}, _minSpeed{-1}, _maxSpeed{-1}, _serviceSpeed{-1}, _nRotorF{-1}, _nSails{-1};

    string shipLabel, comments, shipTyStr, sternStr, hullStr;

    stream >> shipLabel >> comments >> shipTyStr
        >> _L_OA >> _Lpp >> _B >> _T >> _Tmax
        >> _n >> _A_BT >> _sfoc >> _ks
        >> sternStr >> hullStr
        >> _dwt >> _fuelCap
        >> _minSpeed >> _maxSpeed >> _serviceSpeed >> _nRotorF >> _rotorScale >> _nSails >> _sailsScale;

    //cout << line << endl;
    //cout << _nRotorF << endl;

    ship s;

    // Mandatory values
    s.label = shipLabel;
    s.comments = comments;
    s.Loa = _L_OA;
    s.Lpp = _Lpp;
    s.B = _B;
    s.T = _T;
    s.Tmax = _Tmax;
    s.displaceVol = _n;
    s.Abt = _A_BT;
    s.sfoc = _sfoc;
    s.ks = _ks;
    //s.Cs = _Cs;
    //s.windBallastCurve = _windBallastCurve;
    s.dwt = _dwt;
    s.nRotorF = _nRotorF;
    s.rotorScale = _rotorScale;

    s.nSail = _nSails;
    s.sailScale = _sailsScale;

    // Optional values
    s.minSpeed = _minSpeed;
    s.maxSpeed = _maxSpeed;
    s.serviceSpeed = _serviceSpeed;
    s.fuelCap = _fuelCap;

    // Enumerations (stern, hull) and ship type
    if (!sternStr.empty())
    {
        const char sternChar = tolower(sternStr.front());

        if (sternChar == 'p' || sternChar == 'v' || sternChar == 'u')
            s.C_stern = (sternChar == 'p' ? ship::s_P : sternChar == 'v' ? ship::s_V : ship::s_U);
        else
            cout << "Invalid sternStr (" << sternStr << "), keeping default (" << s.C_stern << ")" << endl;
    }

    if (!hullStr.empty())
    {
        const char hullChar = tolower(hullStr.front());

        if (hullChar == 'v' || hullChar == 'u')
            s.hullForm = (hullChar == 'v' ? ship::h_V : ship::h_U);
        else
            cout << "Invalid hullStr (" << hullStr << "), keeping default (" << s.hullForm << ")" << endl;
    }

    if (!shipTyStr.empty())
    {
        const char shipChar = shipTyStr.front();

        if (string("tgbacrp").find(shipChar) != string::npos)
            s.shipTyChar = shipChar;
        else
            cout << "Invalid shipTyStr (" << shipTyStr << "), keeping default (" << s.shipTyChar << ")" << endl;
    }

    if (!s.paramOk())
    {
        cout << "Final range check failed (paramOk), discarding line" << endl;
        return nullopt;
    }

    return make_optional<ship>(s);
}

bool shipV::isLabelUnique() const
{
    vector<string> s;
    s.reserve(vec.size());

    for (const auto& it : vec)
        s.emplace_back(it.label);

    sort(s.begin(), s.end());

    if (s.empty())          return false;
    if (s.front().empty())  return false;
    else                    return adjacent_find(s.begin(), s.end()) == s.end();   
}

void shipV::readShipCsv(const string& fname)
{
    clear();

    if (!filesystem::exists(fname))
    {
        cout << "Using default ships" << endl;

        for (size_t i = 0; i < shipV::defaultShips.size(); i++)
        {
            string line = shipV::defaultShips[i];
            replace(line.begin(), line.end(), ',', ' ');

            if (const auto l = readShipLine(line); l.has_value())
                vec.emplace_back(l.value());
            else
                cout << "Invalid parameter for line " << i + 1 << " in default ship list" << endl;
        }
    }
    else
    {
        ifstream fIN(fname);
        bool fileError{false};
        string line;

        getline(fIN, line); 

        int lineNum{1};

        while (getline(fIN, line) && !line.empty() && !fileError)
        {
            lineNum++;
            replace(line.begin(), line.end(), ',', ' ');

            cout << "Reading line " << lineNum - 1 << endl;
            cout << endl;

            if (const auto l = readShipLine(line); l.has_value())
                vec.emplace_back(l.value());
            else
                cout << "Invalid parameters for line " << lineNum << endl;
        }

        if (vec.empty())
            fileError = true;
       
        if (fileError)
        {
            cout << "Read file failed" << endl;
            clear();
        }
    }

    if (!vec.empty())
        if (!isLabelUnique())
            cout << "shipLabel not unique" << endl;
}


int shipV::getRowCount() const
{
    return static_cast<int>(vec.size());
}

// Not used for now
double shipV::ave(const vector<double>& values, const vector<double>& weights)
{
    if (values.empty())
    {
        cout << "values is empty" << endl;
        return wave::noValDouble;
    }

    if (values.size() != weights.size())
    {
        cout << "values.size() != weights.size()" << endl;
        return wave::noValDouble;
    }

    for (const auto& it : weights)
    {
        if (it < wave::ep)
        {
            cout << "weights < " << wave::ep << endl;
            return wave::noValDouble;
        }
    }

    double result{0};
    for (size_t i = 0; i < values.size(); i++)
        result += (values[i] * weights[i]);

    return result / accumulate(weights.begin(), weights.end(), 0.0);
}

// Not used for now
double shipV::aveAngleNormalise(const double& angle)
{
    double result{angle};

    while (result >= 360)
        result = fmod(result, 360);
    while (result < 0)
        result = 360 + fmod(result, 360);

    return result;
}

// Not used for now
double shipV::aveAngle(const vector<double>& angles, const vector<double>& weights)
{
    // r = right (0 to 180 degrees)
    // l = left (180 to 360 degrees)

    if (angles.empty())
    {
        cout << "angles is empty" << endl;
        return wave::noValDouble;
    }

    if (angles.size() != weights.size())
    {
        cout << "angles.size() != weights.size()" << endl;
        return wave::noValDouble;
    }

    for (const auto& it : weights)
    {
        if (it < wave::ep)
        {
            cout << "weights < " << wave::ep << endl;
            return wave::noValDouble;
        }
    }

    double rAve{0}, lAve{0}, rw_sum{0}, lw_sum{0};

    for (size_t i = 0; i < angles.size() ; i++)
    {
        const double norm{shipV::aveAngleNormalise(angles[i])};

        if (norm >= 180)
        {
            lAve += (norm * weights[i]);
            lw_sum += weights[i];
        }
        else
        {
            rAve += (norm * weights[i]);
            rw_sum += weights[i];
        }
    }

    if (lw_sum > wave::ep)  //important or else div by 0 will happen
        lAve /= lw_sum;

    if (rw_sum > wave::ep)  //important or else div by 0 will happen
        rAve /= rw_sum;

    if (rAve > lAve + 180)
        lAve += 360;

    if (lAve > rAve + 180)
        rAve += 360;

    const double rPart{rAve * rw_sum / (rw_sum + lw_sum)};
    const double lPart{lAve * lw_sum / (rw_sum + lw_sum)};

    return aveAngleNormalise(rPart + lPart);
}

/*
https://stackoverflow.com/questions/491738/how-do-you-calculate-the-average-of-a-set-of-circular-data
double sr::aveAngle(const vector<double>& angles)
{
    // r = right (0 to 180 degrees)
    // l = left (180 to 360 degrees)

    if (angles.empty())
    {
        cout << "angles is empty" << endl;
        return wave::noValDouble;
    }

    double rTotal{0}, lTotal{0};
    int rCount{0}, lCount{0};

    for (const auto& angle : angles)
    {
        const double norm{aveAngleNormalise(angle)};

        if (norm >= 180)
        {
            lTotal += norm;
            lCount++;
        }
        else
        {
            rTotal += norm;
            rCount++;
        }
    }

    double rAvg{rTotal / max(rCount, 1)};
    double lAvg{lTotal / max(lCount, 1)};

    if (rAvg > lAvg + 180)
        lAvg += 360;

    if (lAvg > rAvg + 180)
        rAvg += 360;

    const double rPart{rAvg * rCount / static_cast<double>(rCount + lCount)};
    const double lPart{lAvg * lCount / static_cast<double>(lCount + rCount)};

    return aveAngleNormalise(rPart + lPart);
}
*/
