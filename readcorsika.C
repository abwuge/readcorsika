#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include "TFile.h"
#include <fstream>
#include "TNtupleD.h"
#include <string.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TF1.h"
#include "TTree.h"
#include "TH3F.h"
#include "TCanvas.h"

using namespace std;

int main(int argc, char *argv[])
{
   int i, j;
   int test;
   int id;
   int gen;
   int record_size;
   static float pdata[6552], hobs;

   FILE *input;
   int ievn = -1, numsec;
   double Tolenergy, energy, theta, phi, x, y;
   int ev_n = 0;
   int ibinmax, ibin1, ibin2;
   double xmax;

   int irun, ipart;
   float AltStart, AltObs, AltFirstInt, Momentum[3], momentum[3];
   int counts[5];

   if (argc < 3)
   {
      printf("USAGE:%s infilelist ourroot first last maxevt ........\n", argv[0]);
      return 0;
   }
   char *runlist = argv[1];
   int ffirst = argc > 3 ? atoi(argv[3]) : 0;
   int flast = argc > 4 ? atoi(argv[4]) : -1;
   int nmaxevt = argc > 5 ? atoi(argv[5]) : -1;
   int firstevt = argc > 6 ? atoi(argv[6]) : 0;

   TFile *file = TFile::Open(argv[2], "recreate");

   TCanvas *cc = new TCanvas("cc", "", 2100, 2970);
   cc->Divide(1, 2);
   cc->cd(1);
   gPad->SetLogz(1);
   cc->cd(2);
   gPad->SetLogz(1);
   cc->Print("test.pdf[", "pdf");

   TH2F *hxyg = new TH2F("hxy#gamma", "#gamma;x [m];y [m]", 1000, -2000, 2000, 1000, -2000, 2000);
   TH2F *hxye = new TH2F("hxye", "e;x [m];y [m]", 1000, -2000, 2000, 1000, -2000, 2000);
   TH2F *hxyu = new TH2F("hxy#mu", "#mu;x [m];y [m]", 1000, -2000, 2000, 1000, -2000, 2000);
   TH2F *hxyn = new TH2F("hxyu", "n;x [m];y [m]", 1000, -2000, 2000, 1000, -2000, 2000);
   TH1F *hlongg = new TH1F("hlong#gamma", "#gamma-rays;depth [g/cm^{2}];N#gamma", 120, 0, 1200);
   TH1F *hlonge = new TH1F("hlonge", "e;depth [g/cm^{2}];Ne", 120, 0, 1200);
   TH1F *hlongu = new TH1F("hlongu", "#mu;depth [g/cm^{2}];N#mu", 120, 0, 1200);
   TH1F *hlongh = new TH1F("hlongh", "hadronic;depth [g/cm^{2}];Nh", 120, 0, 1200);
   TH1F *hlongq = new TH1F("hlongq", "all charged particles;depth [g/cm^{2}];Nq", 120, 0, 1200);
   TH1F *hlongce = new TH1F("hlongce", "Cherenkov photons;depth [g/cm^{2}];Nce", 120, 0, 1200);

   double rbins[1000];
   int nrbin = 100;
   rbins[0] = pow(10, -0.5);
   for (int ii = 1; ii <= nrbin; ii++)
   {
      rbins[ii] = pow(10, -0.5 + (3.5 - (-0.5)) / nrbin * ii);
   }
   TH1F *hlaterg = new TH1F("hlater#gamma", "#gamma;R [m];#rho_{#gamma} [m^{-2}]", nrbin, rbins);
   TH1F *hlatere = new TH1F("hlatere", "e;R [m];#rho_{e} [m^{-2}]", nrbin, rbins);
   TH1F *hlateru = new TH1F("hlater#mu", "#mu;R [m];#rho_{#mu} [m^{-2}]", nrbin, rbins);
   TH1F *hlatern = new TH1F("hlatern", "n;R [m];#rho_{n} [m^{-2}]", nrbin, rbins);

   char filename[200];
   int nfile = 0;
   bool islist = strstr(runlist, ".txt");
   ifstream fin;
   if (islist)
   {
      fin.open(runlist, std::ios::in);
      if (!fin.is_open())
         return 0;
   }

   while ((!islist) || (islist && (!fin.eof())))
   {
      if (islist)
      {
         fin.getline(filename, 200);
         nfile++;
         if (nfile < ffirst + 1)
            continue;
         if (flast >= 0 && nfile > flast + 1)
            break;
         if (strlen(filename) < 2)
            break;
      }
      else
      {
         if (nfile > 0)
            break;
         strcpy(filename, runlist);
         nfile++;
      }
      printf("filename = %s\n", filename);
      if ((input = fopen(filename, "rb")) == NULL)
         continue;

      // Get file size
      fseek(input, 0, SEEK_END);
      long int loc0 = ftell(input);
      fseek(input, 0, SEEK_SET);

      int nshift = 0;
      bool ignore = true;
      while (!feof(input))
      {
         fread(&record_size, sizeof(int), 1, input);
         int nword = record_size / 4;
         bool islong = (nword > 5733);


         if (ievn < 0 || (ievn >= 0 && ievn < firstevt))
         {
            int unit = nword / 21; // Each block consists of 21 subblocks
            for (int ii = 0; ii < 21; ii++)
            {
               float var1, var2;
               fread(&var1, 4, 1, input);
               if (ferror(input))
                  goto loop;
               fread(&var2, 4, 1, input);
               fseek(input, -8, SEEK_CUR);
               fseek(input, unit * 4, SEEK_CUR);
               if (var1 >= 211284.0 && var1 <= 211286.0) // RUNH
               {
                  irun = (int)var2;
               }
               if (var1 >= 217432.0 && var1 <= 217434.0) // EVTH
               {
                  ievn = (int)var2;
                  printf("ignore event %d of run %d\n", ievn, irun); // This event may actually not be ignored
               }
               if (ftell(input) > loc0)
                  goto loop;
               // printf("nshift=%d ii=%d var1=%.1f var2=%d word=%d\n",nshift,ii,var1,(int)var2,nword);
            }

            if (ievn >= 0 && ievn < firstevt) // Skip the event
            {
               fseek(input, sizeof(int), SEEK_CUR);
               nshift++;
               continue;
            }
            else // Rewind to the beginning of the event
            {
               fseek(input, -(4 * nword), SEEK_CUR);
            }
         }

         fread(pdata, 4, nword, input);
         if (ftell(input) > loc0)
            goto loop;
         fread(&record_size, sizeof(int), 1, input);
         // printf("ii=%d run=%d evt=%d evtrange=%d,%d\n",nshift,irun,ievn,firstevt,nmaxevt);
         nshift++;
         for (i = 0; i < nword; i = i + (nword / 21))
         {
            test = 0;
            if (pdata[i] >= 211284.0 && pdata[i] <= 211286.0) // RUNH
            {
               test = 1;
            }
            if (pdata[i] >= 217432.0 && pdata[i] <= 217434.0) // EVTH
            {
               test = 2;
            }
            if (pdata[i] >= 52814.0 && pdata[i] <= 52816.0) // LONG
            {
               test = 3;
            }
            if (pdata[i] >= 3396.0 && pdata[i] <= 3398.0) // EVTE
            {
               test = 4;
            }
            if (pdata[i] >= 3300.0 && pdata[i] <= 3302.0) // RUNE
            {
               test = 5;
            }
            // printf("i=%d pdata=%.2f test=%d\n",i,pdata[i],test);
            switch (test)
            {
            case 1: // Run Header
            {
               printf("observation level %lf (cm)\n", pdata[i + 5]);
               irun = pdata[i + 1];
               hobs = pdata[i + 5];
            }
            break;
            case 2: // Event Header
               ievn = (int)pdata[i + 1];
               if (nmaxevt >= 0 && ievn > nmaxevt)
                  goto loop;
               ipart = (int)pdata[i + 2];
               Tolenergy = pdata[i + 3];
               AltStart = pdata[i + 157];
               AltFirstInt = pdata[i + 6];
               Momentum[0] = pdata[i + 7];
               Momentum[1] = pdata[i + 8];
               Momentum[2] = pdata[i + 9];
               theta = pdata[i + 10];
               phi = pdata[i + 11];
               printf("run%d evt%d pid=%d eng=%.1f hobs=%.2f AltStart=%.2f AltInt=%.2f dir={%.2f,%.2f} Momentum={%.2f,%.2f,%.2f}\n", irun, ievn, ipart, Tolenergy, hobs, AltStart, fabs(AltFirstInt), theta * 57.3, phi * 57.3, Momentum[0], Momentum[1], Momentum[2]);

               numsec = 0;
               if (ievn < firstevt)
                  ignore = true;
               else if (nmaxevt > 0 && ievn > nmaxevt)
                  ignore = true;
               else
                  ignore = false;
               if (ignore)
                  break;
               hxyg->Reset();
               hxye->Reset();
               hxyu->Reset();
               hxyn->Reset();
               hlongg->Reset();
               hlonge->Reset();
               hlongu->Reset();
               hlongh->Reset();
               hlongq->Reset();
               hlongce->Reset();
               hlaterg->Reset();
               hlatere->Reset();
               hlateru->Reset();
               hlatern->Reset();
               break;
            case 3: // Long
               if (ignore)
                  break;
               for (int iblk = 0; iblk < 26; iblk++)
               {
                  float depth = pdata[i + iblk * 10 + 13];
                  // if(depth<=0) continue;
                  // float subpid=pdata[iblk];
                  float nph = pdata[i + iblk * 10 + 14];
                  float nep = pdata[i + iblk * 10 + 15];
                  float nem = pdata[i + iblk * 10 + 16];
                  float nup = pdata[i + iblk * 10 + 17];
                  float num = pdata[i + iblk * 10 + 18];
                  float nhr = pdata[i + iblk * 10 + 19];
                  float nq = pdata[i + iblk * 10 + 20];
                  float nce = pdata[i + iblk * 10 + 22];
                  hlongg->Fill(depth, nph);
                  hlonge->Fill(depth, nep + nem);
                  hlongu->Fill(depth, nup + num);
                  hlongh->Fill(depth, nhr);
                  hlongq->Fill(depth, nq);
                  hlongce->Fill(depth, nce);
                  // int whichstep=(int)(depth/1.+0.5)-1;
                  // printf("idep=%d depth=%.2f n={%f,%f,%f,%f,%f,%f,%f}\n",iblk,depth,nph,nep,nem,nup,num,nq,nce);
               }
               break;
            case 4: // Event End
               if (ignore)
                  break;
               counts[0] = pdata[i + 2]; // photon
               counts[1] = pdata[i + 3]; // electrons
               counts[2] = pdata[i + 4]; // hadrons
               counts[3] = pdata[i + 5]; // muons
               counts[4] = pdata[i + 6]; // all particles
               ev_n++;
               // printf("counts={%d,%d,%d,%d,%d}\n",counts[0],counts[1],counts[2],counts[3],counts[4]);
               // process the histograms

               cc->cd(1);
               gPad->SetLogy(0);
               hxyg->Draw("colz");
               cc->cd(2);
               gPad->SetLogy(0);
               hxye->Draw("colz");
               cc->Print("test.pdf");

               cc->cd(1);
               hxyu->Draw("colz");
               cc->cd(2);
               hxyn->Draw("colz");
               cc->Print("test.pdf");

               for (int ibin = 1; ibin <= hlatere->GetNbinsX(); ibin++)
               {
                  double r1 = hlatere->GetXaxis()->GetBinLowEdge(ibin);
                  double r2 = hlatere->GetXaxis()->GetBinUpEdge(ibin);
                  double area = TMath::Pi() * (r2 * r2 - r1 * r1);

                  double ng = hlaterg->GetBinContent(ibin);
                  if (ng > 0)
                  {
                     hlaterg->SetBinContent(ibin, ng / area);
                     hlaterg->SetBinError(ibin, hlaterg->GetBinError(ibin) / area);
                  }

                  double ne = hlatere->GetBinContent(ibin);
                  if (ne > 0)
                  {
                     hlatere->SetBinContent(ibin, ne / area);
                     hlatere->SetBinError(ibin, hlatere->GetBinError(ibin) / area);
                  }

                  double nu = hlateru->GetBinContent(ibin);
                  if (nu > 0)
                  {
                     hlateru->SetBinContent(ibin, nu / area);
                     hlateru->SetBinError(ibin, hlateru->GetBinError(ibin) / area);
                  }

                  double nn = hlatern->GetBinContent(ibin);
                  if (nn > 0)
                  {
                     hlatern->SetBinContent(ibin, nn / area);
                     hlatern->SetBinError(ibin, hlatern->GetBinError(ibin) / area);
                  }
               }

               cc->cd(1);
               gPad->SetLogx(1);
               gPad->SetLogy(1);
               hlaterg->Draw("hist");
               cc->cd(2);
               gPad->SetLogx(1);
               gPad->SetLogy(1);
               hlatere->Draw("hist");
               cc->Print("test.pdf");

               cc->cd(1);
               hlateru->Draw("hist");
               cc->cd(2);
               hlatern->Draw("hist");
               cc->Print("test.pdf");

               cc->cd(1);
               gPad->SetLogx(0);
               hlongg->Draw("hist");
               cc->cd(2);
               gPad->SetLogx(0);
               hlonge->Draw("hist");
               cc->Print("test.pdf");

               cc->cd(1);
               hlongu->Draw("hist");
               cc->cd(2);
               hlongh->Draw("hist");
               cc->Print("test.pdf");

               cc->cd(1);
               hlongq->Draw("hist");
               cc->cd(2);
               hlongce->Draw("hist");
               cc->Print("test.pdf");

               if (nmaxevt >= 0 && ievn >= nmaxevt)
                  goto loop;

               break;
            case 5: // Run End
            {
               printf("RUN End\n");
               goto loop;
            }
            break;
            default: // Particle Block
               if (ignore)
                  break;
               for (j = i; j < i + (islong ? 311 : 272); j = j + (islong ? 8 : 7))
               {
                  if (theta <= 90. / (180. / TMath::Pi()))
                  {
                     int index0 = (int)(pdata[j]);
                     id = int(pdata[j] / 1000.); // sub_particle ID
                     double wavelength = 0;
                     double weight = 1;
                     if (id >= 9900)
                     {                    // Cerenkov light
                        x = pdata[j + 1]; // cm
                        y = pdata[j + 2];
                        wavelength = fabs((float)(pdata[j + 7]));
                        weight = ((index0 % 100000) - 1) / 10.;

                        double hgen = pdata[j + 6];
                        double dirc[3] = {pdata[j + 3], pdata[j + 4], -sqrt(1 - pow(pdata[j + 3], 2) - pow(pdata[j + 4], 2))};
                        double zenith = acos(-dirc[2] / sqrt(dirc[0] * dirc[0] + dirc[1] * dirc[1] + dirc[2] * dirc[2]));
                     }
                     else
                     { // Secondary particles
                        momentum[0] = pdata[j + 1];
                        momentum[1] = pdata[j + 2];
                        momentum[2] = pdata[j + 3];
                        energy = sqrt(momentum[0] * momentum[0] + momentum[1] * momentum[1] + momentum[2] * momentum[2]); // energy in  GeV
                        x = pdata[j + 4];                                                                                 // cm
                        y = pdata[j + 5];
                        weight = islong ? pdata[j + 7] : 1;
                        numsec++;
                        // printf("id=%d xy={%.2lf,%.2lf} momentum={%.2lf,%.2lf,%.2lf} energy=%.2le weight=%.2lf\n",id,x,y,momentum[0],momentum[1],momentum[2],energy,weight);

                        double rr = sqrt(x * x + y * y);
                        if (id == 1)
                        { // gamma
                           hxyg->Fill(x / 100, y / 100, weight);
                           hlaterg->Fill(rr / 100, weight);
                        }
                        if (id == 2 || id == 3)
                        { // e+ & e-
                           hxye->Fill(x / 100, y / 100, weight);
                           hlatere->Fill(rr / 100, weight);
                        }
                        if (id == 5 || id == 6)
                        { // mu+ & mu-
                           hxyu->Fill(x / 100, y / 100, weight);
                           hlateru->Fill(rr / 100, weight);
                        }
                        if (id == 13)
                        { // n
                           hxyn->Fill(x / 100, y / 100, weight);
                           hlatern->Fill(rr / 100, weight);
                        }
                     }
                  }
               }
               break;
            }
         }
      }
   loop:
      fclose(input);
      printf("..........file %s was closed ...........\n", filename);
      printf("........      %d events         ........\n\n", ev_n);
   }

   file->cd();
   file->Close();

   cc->Print("test.pdf]", "pdf");

   return 0;
}
