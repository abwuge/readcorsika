#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <vector>
#include <string>
#include <limits>

#include "TMath.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
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
};
const float Point::epsilon = 1e-2;

Point getCorePosition(const string &fileName, size_t &particleID, float &energy, float &theta, float &phi);

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
    system(("rm -r " + pngPath).c_str());
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
    vector<vector<vector<float>>> allData(particleTypes.size());

    for (size_t i = 0; i < particleTypes.size(); ++i)
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
    vector<size_t> allIdx(particleTypes.size(), 0);
    for (size_t fileIdx = 0; fileIdx < particleTypes.size(); ++fileIdx)
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
    for (size_t fileIdx = 0; fileIdx < particleTypes.size(); ++fileIdx)
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
    { // Calculate the average x, y position of the particles
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
    }

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

    TLegend *legend = new TLegend(0.15, 0.8, 0.72, 0.9);
    { // Set the legend
        legend->SetFillColor(0);
        legend->SetBorderSize(0);
        legend->SetNColumns(3);
        legend->SetTextSize(0.05);
    }
    for (size_t i = 0; i < particleTypes.size(); ++i)
    {
        TMarker *dummyMarker = new TMarker(0.5, 0.5, 20);
        dummyMarker->SetMarkerColor(particleTypes[i].color);
        legend->AddEntry(dummyMarker, particleTypes[i].name.c_str(), "p");
    }

    string titleNames[4] = {
        "x-z Plane",
        "y-z Plane",
        "x-y Plane",
        "3D View",
    };
    for (size_t i = 0; i < 4; ++i)
    {
        canvas->cd(i + 1);
        TView *view = TView::CreateView(1);
        view->SetRange(xMin, yMin, zMin, xMax, yMax, zMax);
        TPaveText *title = new TPaveText(0.1, 0.9, 0.9, 1, "NDC");
        title->SetTextAlign(22);
        title->SetTextSize(0.06);
        title->SetFillColor(0);
        title->SetLineColor(0);
        title->SetShadowColor(0);
        title->AddText(titleNames[i].c_str());
        title->Draw();

        legend->Draw();
    }

    TPaveText *time = new TPaveText(0.75, 0.8, 0.95, 0.9, "NDC");
    { // Set the time display
        time->SetTextAlign(12);
        time->SetTextSize(0.05);
        time->SetTextFont(42);
        time->SetFillColor(0);
        time->SetLineColor(0);
        time->SetShadowColor(0);
    }

    // Plot particles
    TGraph2D *blankGraph = new TGraph2D();
    blankGraph->SetPoint(0, 0, 0, 0);
    blankGraph->SetMarkerSize(0);
    cout << "Generating PNG files ..." << endl;
    size_t frame = 0;
    for (double currentT = tStart; currentT <= tMax; currentT += deltaT)
    {
        vector<TGraph2D *> graphs2D(particleTypes.size(), nullptr);
        vector<TGraph *> graphsXZ(particleTypes.size(), nullptr);
        vector<TGraph *> graphsYZ(particleTypes.size(), nullptr);
        vector<TGraph *> graphsXY(particleTypes.size(), nullptr);

        for (size_t i = 0; i < 4; ++i)
        {
            canvas->cd(i + 1);
            gPad->Clear();
        }

        // Iterate over each file
        for (size_t fileIdx = 0; fileIdx < particleTypes.size(); ++fileIdx)
        {
            size_t cnt = 0;

            const auto &data = allData[fileIdx];

            auto graphInit = [](TGraph *graph, int color)
            {
                graph->SetMarkerStyle(20);
                graph->SetMarkerSize(0.5);
                graph->SetMarkerColor(color);
            };
            auto graph2DInit = [](TGraph2D *graph, int color)
            {
                graph->SetMarkerStyle(20);
                graph->SetMarkerSize(0.5);
                graph->SetMarkerColor(color);
            };

            graphsXY[fileIdx] = new TGraph();
            graphInit(graphsXY[fileIdx], particleTypes[fileIdx].color);
            graphsXZ[fileIdx] = new TGraph();
            graphInit(graphsXZ[fileIdx], particleTypes[fileIdx].color);
            graphsYZ[fileIdx] = new TGraph();
            graphInit(graphsYZ[fileIdx], particleTypes[fileIdx].color);
            graphs2D[fileIdx] = new TGraph2D();
            graph2DInit(graphs2D[fileIdx], particleTypes[fileIdx].color);

            for (auto &idx = allIdx[fileIdx]; idx < data.size(); ++idx)
            {
                const auto &particleInfo = data[idx];
                if (particleInfo[9] < currentT + deltaT)
                {
                    const auto &xEnd = 1e-2 * particleInfo[6], &yEnd = 1e-2 * particleInfo[7], &zEnd = 1e-2 * particleInfo[8];

                    if (xEnd < xMin || xEnd > xMax || yEnd < yMin || yEnd > yMax || zEnd < zMin || zEnd > zMax || cnt > 1e3)
                        continue;

                    // Add the particle to the graph
                    graphsXZ[fileIdx]->SetPoint(cnt, xEnd, zEnd);
                    graphsYZ[fileIdx]->SetPoint(cnt, yEnd, zEnd);
                    graphsXY[fileIdx]->SetPoint(cnt, xEnd, yEnd);
                    graphs2D[fileIdx]->SetPoint(cnt++, xEnd, yEnd, zEnd);
                }
                else
                    break;
            }
        }

        // Draw the accumulated data
        time->Clear();
        time->AddText(TString::Format("Time: %.2f #mus", currentT * 1e6));

        auto drawGraphs = [&](const vector<TGraph *> &graphs, size_t i)
        {
            canvas->cd(i);

            if (i == 1) // x-z plane
                gPad->DrawFrame(xMin, zMin, xMax, zMax, ";X;Z");
            else if (i == 2) // y-z plane
                gPad->DrawFrame(yMin, zMin, yMax, zMax, ";Y;Z");
            else if (i == 3) // x-y plane
                gPad->DrawFrame(xMin, yMin, xMax, yMax, ";X;Y");

            TPaveText *title = new TPaveText(0.1, 0.9, 0.9, 1, "NDC");
            { // Set the title
                title->SetTextAlign(22);
                title->SetTextSize(0.06);
                title->SetFillColor(0);
                title->SetLineColor(0);
                title->SetShadowColor(0);
                title->AddText(titleNames[i - 1].c_str());
            }
            title->Draw();

            legend->Draw();

            for (size_t fileIdx = 0; fileIdx < particleTypes.size(); ++fileIdx)
                if (graphs[fileIdx]->GetN())
                    graphs[fileIdx]->Draw("P same");

            time->Draw();
        };
        auto drawGraphs2D = [&](const vector<TGraph2D *> &graphs, size_t i)
        {
            canvas->cd(i);

            TView *view = TView::CreateView(1);
            view->SetRange(xMin, yMin, zMin, xMax, yMax, zMax);
            view->ShowAxis();

            TPaveText *title = new TPaveText(0.1, 0.9, 0.9, 1, "NDC");
            { // Set the title
                title->SetTextAlign(22);
                title->SetTextSize(0.06);
                title->SetFillColor(0);
                title->SetLineColor(0);
                title->SetShadowColor(0);
                title->AddText(titleNames[i - 1].c_str());
            }
            title->Draw();

            legend->Draw();

            blankGraph->Draw("P same");
            for (size_t fileIdx = 0; fileIdx < particleTypes.size(); ++fileIdx)
                if (graphs[fileIdx]->GetN())
                    graphs[fileIdx]->Draw("P same");
            time->Draw();
        };

        drawGraphs(graphsXZ, 1);
        drawGraphs(graphsYZ, 2);
        drawGraphs(graphsXY, 3);
        drawGraphs2D(graphs2D, 4);

        size_t pointNum = 0;
        for (size_t fileIdx = 0; fileIdx < graphs2D.size(); ++fileIdx)
            pointNum += graphs2D[fileIdx]->GetN();

        printf("Frame %2d (time: %6.0fns): add %3d point%s;\n", ++frame, currentT * 1e9, pointNum, pointNum < 2 ? "" : "s");

        canvas->Update();
        canvas->Print(TString::Format("%s%.10f.png", pngPath.c_str(), currentT).Data(), "png");

        for (auto graph : graphs2D)
            delete graph;
    }
    cout << "PNG generated successfully!" << endl;
    cout << "Generating GIF file ..." << endl;
    string gifName = TString::Format("%sPID%d_E%.1f_T%.1f_P%.1fp.gif", PNGPath.c_str(), particleID, energy, theta * 57.3, phi * 57.3).Data();
    system(("convert " + pngPath + "*.png " + gifName).c_str());
    cout << "GIF generated successfully: " << gifName << endl;

    // Cleanup
    cout << "Cleaning ..." << endl;
    system(("rm -r " + pngPath).c_str());

    delete canvas;
    cout << "Cleaned" << endl;

    return 0;
}

Point getCorePosition(const string &fileName, size_t &particleID, float &energy, float &theta, float &phi)
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
