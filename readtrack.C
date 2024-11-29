#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <string>
#include <limits>
#include <unordered_set>

#include "TMath.h"
#include "TCanvas.h"
#include "TPolyLine3D.h"
#include "TObjArray.h"
#include "TView.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TString.h"


using namespace std;

struct ParticleInfo
{
    string suffix;
    string name;
    int color;
};

struct Point
{
    float x, y, z;

    static const float epsilon;

    Point(float x_val, float y_val, float z_val) : x(x_val), y(y_val), z(z_val) {}

    struct Hash
    {
        size_t operator()(const Point &p) const
        {
            size_t h1 = std::hash<int>()(static_cast<int>(p.x / epsilon));
            size_t h2 = std::hash<int>()(static_cast<int>(p.y / epsilon));
            size_t h3 = std::hash<int>()(static_cast<int>(p.z / epsilon));
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };

    struct Equal
    {
        bool operator()(const Point &lhs, const Point &rhs) const
        {
            return std::abs(lhs.x - rhs.x) < epsilon && std::abs(lhs.y - rhs.y) < epsilon && std::abs(lhs.z - rhs.z) < epsilon;
        }
    };
};
const float Point::epsilon = 1e-2;

Point getCorePosition(const string &fileName, size_t &particleID, float& energy, float &theta, float &phi);

int main(int argc, char *argv[])
{
    gErrorIgnoreLevel = kWarning;

    size_t randSeed = 0;
    srand(randSeed);

    if (argc < 2)
    {
        cout << "USAGE: " << argv[0] << " DATnnnnnn" << endl;
        cout << "e.g. " << argv[0] << " Data/DAT000001" << endl;
        return -1;
    }

    string baseFileName = argv[1]; // Base name of the input files
    string pngPath = baseFileName.substr(0, baseFileName.find_last_of("/")) + "/PNG/" + baseFileName.substr(baseFileName.find_last_of("/") + 1) + "/";
    string PNGPath = baseFileName.substr(0, baseFileName.find_last_of("/")) + "/PNG/";
    mkdir(PNGPath.c_str(), 0755);
    mkdir(pngPath.c_str(), 0755);

    const size_t nvar = 10; // Number of variables per record
    vector<ParticleInfo> particleTypes = {
        {"track_em", "electromagnetic", kRed},
        {"track_mu", "muonic", kBlue},
        {"track_hd", "hadronic", kGreen},
    };

    // Prepare file names
    vector<string> fileNames;
    for (const auto &particle : particleTypes)
    {
        fileNames.push_back(baseFileName + "." + particle.suffix);
    }

    // Load all data into a 3D array: [file][record_index][variable]
    vector<vector<vector<float>>> allData(fileNames.size());

    for (size_t i = 0; i < fileNames.size(); ++i)
    {
        const auto &fileName = fileNames[i];
        ifstream inputFile(fileName, ios::binary);
        if (!inputFile.is_open())
        {
            cerr << "Error opening file: " << fileName << endl;
            return -1;
        }

        // Read all records from the file
        size_t recordSize = 0;
        vector<float> record(nvar);
        while (inputFile.read(reinterpret_cast<char *>(&recordSize), sizeof(int)))
        {
            if (recordSize != nvar * sizeof(float))
            {
                cerr << "Expected record size: " << nvar * sizeof(float) << ", but got " << recordSize << endl;
                cerr << "Invalid record size in file: " << fileName << endl;
                return -1;
            }

            inputFile.read(reinterpret_cast<char *>(record.data()), recordSize);
            allData[i].push_back(record);

            inputFile.read(reinterpret_cast<char *>(&recordSize), sizeof(int));
        }

        inputFile.close();
    }

    // Sort data by end time
    for (auto &data : allData)
    {
        sort(data.begin(), data.end(), [](const std::vector<float> &a, const std::vector<float> &b)
             { return a[9] < b[9]; });
    }

    // set the time range
    double tMin = 1e64, tMax = 0;
    for (const auto &data : allData)
        tMin = min(tMin, (double)data[0][9]), tMax = max(tMax, (double)data.back().at(9));
    double tStart = 1e64;
    for (const auto &data : allData) // Find the starting time: just before the first mother particle -> secondary particle
        if (data[0][9] == tMin)
            for (size_t i = 0; i < data.size(); ++i)
                if (data[i][9] == data[i + 1][9])
                {
                    if (i - 1 >= 0)
                        tStart = min(tStart, (double)data[i - 1][9]);
                    else
                        tStart = min(tStart, (double)data[i][9]);
                    break;
                }
    double deltaT = (tMax - tStart) / 45; // Time window
    tStart -= 4 * deltaT;                 // Start 4 time windows before the first particle

    // Find the starting index for each file
    vector<size_t> allIdx(fileNames.size(), 0);
    for (size_t fileIdx = 0; fileIdx < fileNames.size(); ++fileIdx)
    {
        const auto &data = allData[fileIdx];
        auto it = lower_bound(data.begin(), data.end(), tStart, [](const vector<float> &vec, float value)
                              { return vec[9] < value; });
        allIdx[fileIdx] = it - data.begin();
        cout << fileNames[fileIdx] << " start index: " << allIdx[fileIdx] << endl;
    }
    // set the spatial range
    size_t particleID = 0;
    float energy = 0, theta = 0, phi = 0;
    Point corePosition = getCorePosition(baseFileName, particleID, energy, theta, phi);
    // set z range
    double zMin = corePosition.z, zMax = 0, xTop = 0, yTop = 0;
    for (size_t fileIdx = 0; fileIdx < fileNames.size(); ++fileIdx)
        if (allIdx[fileIdx])
        {
            const auto particleInfo = allData[fileIdx][allIdx[fileIdx]];
            xTop = particleInfo[2], yTop = particleInfo[3], zMax = particleInfo[4];
        }
    if (!zMax)
        for (const auto &data : allData)
            zMax = max(zMax, (double)data[0][4]);
    xTop *= 1e-2, yTop *= 1e-2, zMax *= 1e-2;

    // set x, y range
    long double xAvg = 0, yAvg = 0;
    unsigned long pcnt = 0;
    for (const auto &data : allData)
        for (const auto &particleInfo : data)
        {
            const auto &xEnd = 1e-2 * particleInfo[6], &yEnd = 1e-2 * particleInfo[7], &zEnd = 1e-2 * particleInfo[8];
            if (zEnd - zMin < 0.1 * zMax)
            {
                xAvg += xEnd;
                yAvg += yEnd;
                ++pcnt;
            }
        }
    xAvg /= pcnt, yAvg /= pcnt;

    double xMin = xAvg - 2e3, xMax = xAvg + 2e3, yMin = yAvg - 2e3, yMax = yAvg + 2e3;
    if (xTop < xMin)
        xMin = xTop;
    else if (xTop > xMax)
        xMax = xTop;
    if (yTop < yMin)
        yMin = yTop;
    else if (yTop > yMax)
        yMax = yTop;
    // Slightly enlarge range
    xMin -= 1 + 1e-3 * xMin, yMin -= 1 + 1e-3 * yMin, zMin -= 1 + 1e-3 * zMin;
    xMax += 1 + 1e-3 * xMax, yMax += 1 + 1e-3 * yMax, zMax += 1 + 1e-3 * zMax;

    printf("Range:\n x: %.2f ~ %.2f\n y: %.2f ~ %.2f\n z: %.2f ~ %.2f\n", xMin, xMax, yMin, yMax, zMin, zMax);

    // Create the canvas and views for plotting
    TCanvas *canvas = new TCanvas("canvas", "Particle Tracks", 1414, 1000);
    canvas->Divide(2, 2);

    TLegend *legend = new TLegend(0.05, 0.8, 0.6, 0.9);
    legend->SetFillColor(0);
    legend->SetBorderSize(0);
    legend->SetNColumns(3);
    legend->SetTextSize(0.05);
    for (size_t i = 0; i < particleTypes.size(); ++i) {
        TMarker *dummyMarker = new TMarker(0.5, 0.5, 20);
        dummyMarker->SetMarkerColor(particleTypes[i].color);
        legend->AddEntry(dummyMarker, particleTypes[i].name.c_str(), "p");
    }

    vector<TView *> views(4);
    Double_t viewAngles[4][3] = {
        {-90, 90, 0}, // x-z plane
        {0, 90, 0},   // y-z plane
        {0, 0, -90},  // x-y plane
        {45, 45, 0},  // 3D view
    };
    string titles[4] = {
        "x-z Plane",
        "y-z Plane",
        "x-y Plane",
        "3D View",
    };
    for (size_t i = 0; i < 4; ++i)
    {
        canvas->cd(i + 1);
        views[i] = TView::CreateView(1);
        views[i]->SetRange(xMin, yMin, zMin, xMax, yMax, zMax);

        Int_t irep = 0;
        views[i]->SetView(viewAngles[i][0], viewAngles[i][1], viewAngles[i][2], irep);
        views[i]->ShowAxis();

        TPaveText *title = new TPaveText(0.1, 0.9, 0.9, 1, "NDC");
        title->SetTextAlign(22);
        title->SetTextSize(0.06); 
        title->SetFillColor(0);
        title->SetLineColor(0);
        title->SetShadowColor(0);
        title->AddText(titles[i].c_str());
        title->Draw();

        legend->Draw();
    }

    TPaveText *time = new TPaveText(0.7, 0.8, 0.9, 0.9, "NDC");
    time->SetTextAlign(12);
    time->SetTextSize(0.05);
    time->SetTextFont(42);
    time->SetFillColor(0);
    time->SetLineColor(0);
    time->SetShadowColor(0);

    TObjArray *plot = new TObjArray();

    // Record which line is used
    unordered_set<Point, Point::Hash, Point::Equal> usedPoints;
    // add the first point
    for (size_t fileIdx = 0; fileIdx < fileNames.size(); ++fileIdx)
    {
        const auto &particleInfo = allData[fileIdx][allIdx[fileIdx]];
        usedPoints.emplace(Point(1e-2 * particleInfo[2], 1e-2 * particleInfo[3], 1e-2 * particleInfo[4]));
    }

    // Plot particle tracks
    cout << "Generating PNG files ..." << endl;
    size_t frame = 0;
    bool hasDrawed = false;
    double currentT = tStart;
    while (currentT <= tMax)
    {
        size_t cnt = 0;
        plot->Clear();
        // Iterate over each file
        for (size_t fileIdx = 0; fileIdx < fileNames.size(); ++fileIdx)
        {
            const auto &data = allData[fileIdx];
            for (auto &idx = allIdx[fileIdx]; idx < data.size(); ++idx)
            {
                const auto &particleInfo = data[idx];
                if (particleInfo[9] < currentT + deltaT)
                {
                    const auto &xStart = 1e-2 * particleInfo[2], &yStart = 1e-2 * particleInfo[3], &zStart = 1e-2 * particleInfo[4];
                    const auto &xEnd = 1e-2 * particleInfo[6], &yEnd = 1e-2 * particleInfo[7], &zEnd = 1e-2 * particleInfo[8];
                    // Clean
                    auto cleanLines = [&]()
                    {
                        if (!usedPoints.count(Point(xStart, yStart, zStart)))
                            return true;
                        if (zStart - zEnd < -1) // Filter out the particles that is moving upwards
                            return true;
                        if (xEnd < xMin || xEnd > xMax)
                            return true;
                        if (yEnd < yMin || yEnd > yMax)
                            return true;
                        if (++cnt > 2e3 && rand() % 100 < 75) // Randomly remove some lines
                            return true;
                        return false;
                    };

                    if (hasDrawed && cleanLines())
                        continue;

                    if (!hasDrawed)
                        if (!cleanLines())
                            hasDrawed = true;

                    usedPoints.emplace(Point(xEnd, yEnd, zEnd));

                    if (!hasDrawed && cleanLines())
                        continue;

                    // Add the particle track to the plot
                    TPolyLine3D *line = new TPolyLine3D(2);
                    line->SetPoint(0, xStart, yStart, zStart); // Start point
                    line->SetPoint(1, xEnd, yEnd, zEnd);       // End point
                    line->SetLineColor(particleTypes[fileIdx].color);
                    line->SetLineWidth(1);
                    plot->Add(line);
                }
                else
                    break;
            }
        }

        // Draw the accumulated data

        size_t entries = plot->GetEntries();
        printf("Frame %d (time: %.0fns): add %d line%s;\n", ++frame, currentT * 1e9, entries, entries < 2 ? "" : "s");

        time->Clear();
        time->AddText(TString::Format("Time: %.2f #mus", currentT * 1e6));
        for (size_t i = 0; i < 4; ++i)
        {
            canvas->cd(i + 1);
            plot->Draw("same");
            time->Draw();
        }
        canvas->Update();
        canvas->Print(TString::Format("%s%.10f.png", pngPath.c_str(), currentT).Data(), "png");

        currentT += deltaT; // Move to the next time window
    }
    canvas->Print(TString::Format("%sPID%d_E%.1f_T%.1f_P%.1f.png", PNGPath.c_str(), particleID, energy, theta * 57.3, phi * 57.3).Data(), "png");
    cout << "PNG generated successfully!" << endl;
    cout << "Generating GIF file ..." << endl;
    string gifName = TString::Format("%sPID%d_E%.1f_T%.1f_P%.1f.gif", PNGPath.c_str(), particleID, energy, theta * 57.3, phi * 57.3).Data();
    system(("convert " + pngPath + "*.png " + gifName).c_str());
    cout << "GIF generated successfully: " << gifName << endl;

    // Cleanup
    cout << "Cleaning ..." << endl;
    system(("rm -r " + pngPath).c_str());
    delete plot;
    delete canvas;
    cout << "Cleaned" << endl;

    return 0;
}

Point getCorePosition(const string &fileName, size_t &particleID, float& energy, float &theta, float &phi)
{
    Point corePosition(0, 0, 0);

    ifstream inputFile(fileName, ios::binary);
    if (!inputFile.is_open())
    {
        cerr << "Error opening file: " << fileName << endl;
        return corePosition;
    }

    // Read records
    size_t recordSize = 0;
    if (!inputFile.read(reinterpret_cast<char *>(&recordSize), sizeof(int))) // Read head record size
    {
        cerr << "Error reading record size from file: " << fileName << endl;
        return corePosition;
    }

    size_t nword = recordSize / 4;
    bool isLong = nword > 5733;
    vector<float> record(nword);
    if (!inputFile.read(reinterpret_cast<char *>(record.data()), recordSize)) // Read the first record
    {
        cerr << "Error reading record from file: " << fileName << endl;
        return corePosition;
    }
    size_t type = 0;
    size_t subword = nword / 21;
    auto &xCore = corePosition.x, &yCore = corePosition.y, &zCore = corePosition.z;
    size_t cnt = 0;
    for (int i = 0; i < nword; i += subword)
    {
        type = 0;
        if (record[i] >= 211284.0 && record[i] <= 211286.0) // RUNH
            type = 1;
        if (record[i] >= 217432.0 && record[i] <= 217434.0) // EVTH
            type = 2;
        if (record[i] >= 52814.0 && record[i] <= 52816.0) // LONG
            type = 3;
        if (record[i] >= 3396.0 && record[i] <= 3398.0) // EVTE
            type = 4;
        if (record[i] >= 3300.0 && record[i] <= 3302.0) // RUNE
            type = 5;

        switch (type)
        {
        case 1: // RUN HEADER
        {
            zCore = record[i + 5]; // height of observation as z core position
            break;
        }
        case 2: // EVENT HEADER
        {
            particleID = record[i + 2];
            energy = record[i + 3];
            theta = record[i + 10];
            phi = record[i + 11];
            if (particleID)
                printf("Particle ID: %d, Energy: %.1f, Theta: %.2f, Phi: %.2f\n", particleID, energy, theta, phi);
            break;
        }
        case 3: // LONG
        {
            // ignore
            break;
        }
        case 4: // EVENT END
        {
            // ignore
            break;
        }
        case 5: // RUN END
        {
            // ignore
            break;
        }
        default: // Particle Block
        {
            // for (int j = i; j < i + subword; j += 7 + isLong)
            // {
            //     size_t id = record[j] / 1000;
            //     if (id < 99e5) // Secondary particles
            //     {
            //         float x = record[j + 4];
            //         float y = record[j + 5];
            //         xCore = (xCore * cnt + x) / (cnt + 1);
            //         yCore = (yCore * cnt + y) / (cnt + 1);
            //         ++ cnt;
            //     }
            // }
        }
        }
    }

    if (!inputFile.read(reinterpret_cast<char *>(&recordSize), sizeof(int))) // Read tail record size
    {
        cerr << "Error reading record size from file: " << fileName << endl;
        return Point(-1, -1, -1);
    }

    // This file should only contain one record, so we can close it now
    inputFile.close();

    // Convert to meters
    xCore *= 1e-2;
    yCore *= 1e-2;
    zCore *= 1e-2;

    return corePosition;
}
