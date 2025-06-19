#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <chrono>
#include <cstdio>     // remove, rename

#include "TTree.h"
#include "TTreeIndex.h"
#include "TFile.h"
//#define WRITE_CSV

using namespace std;

struct CascadeData {
  vector<int32_t> primaryEnergies; //Holds the primary particle enrgy
  vector<float> primary_xLocs;
  vector<float> primary_yLocs;
  vector<float> primary_zLocs;
  vector<vector<float>> recoil_xLocs;
  vector<vector<float>> recoil_yLocs;
  vector<vector<float>> recoil_zLocs;
  vector<vector<int16_t>> recoil_nVacs;
  
};

//Functions
void ProcessThrow(CascadeData& data, TTree* tree, int32_t& energy, vector<float>& xs, vector<float>& ys, vector<float>& zs, vector<int16_t>& nVacs, float clusteringDistance_nm, int32_t binSize_eV, int maxEntriesPerBin, int32_t maxEnergy_eV, int32_t minEnergy_eV);

//Takes dx,dy,dz unit vector and computes the matrix to rotate arbitrary points in that space to 1,0,0
array<array<float, 3>, 3> ComputeRotationMatrix(float dx, float dy, float dz);

std::map<int, int> energyBinCounter;

int main(int argc, char* argv[]) {

  //Srim settings  
  int32_t minEnergy_eV = 0;
  int32_t maxEnergy_eV = 5000e3;  //We start our SRIM sim slightly above this, at least 100 keV to allow for "burn in"


  int maxEntriesPerBin = 5e2;  //Avoid filling up too much data at lower energys
  int32_t binSize_eV = 1e3;       //Bin size
  
  //
  float clusteringDistance_nm = 4.026*0.1*1.5;  //4.026 Angstroms is the lattice distance, but we're using 1.5 lattice spacings
  vector<string> atomsToTrack = {"09"};         //Don't care about Li vacancies, not optical
  //--------------------//
  //Parse cmd line input//
  //--------------------//
  if ((argc != 2)&&(argc != 3)) {
      cout << "Usage: " << argv[0] << " <COLLISON.txt> <optiinal output filename>" << endl;
      return 1;
  }
  string filename = argv[1];

  string outputFilename = "trimTracks.root";
  if (argc == 3) {
    outputFilename = argv[2];
  }


  cout << "Input file: " << filename << endl;
  //--------------------//
  //Make the output file//
  //--------------------//
  TFile* outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  TTree* trimTree = new TTree("trimTree", "Clustered recoil tracks");
  
  int32_t energy = 0;
  vector<float> xs;
  vector<float> ys;
  vector<float> zs;
  vector<int16_t> nVacs;

  trimTree->Branch("energy_eV", &energy, "energy_eV/i");
  trimTree->Branch("xs_nm", &xs);
  trimTree->Branch("ys_nm", &ys);
  trimTree->Branch("zs_nm", &zs);
  trimTree->Branch("nVacs", &nVacs);

  trimTree->SetAutoFlush(0);  // disable autosizing
  trimTree->SetBasketSize("*", 1000000); //trying to debug...

  //-------------//
  //Open the file//
  //-------------// 
  ifstream infile(filename);
  if (!infile.is_open()) {
      cerr << "Error: could not open file " << filename << endl;
      return 1;
  }
  
  //--------------------------------//
  //Load collision data into vectors//
  //--------------------------------//
  string line;
  int lineCounter = 0;

  // Data vectors
  CascadeData data;
  vector<float> tmpRecoil_xLocs;
  vector<float> tmpRecoil_yLocs;
  vector<float> tmpRecoil_zLocs;
  vector<int16_t> tmpRecoil_nVacs;

  int throwNum = 0; //Only for tracking how far in the processing we are
  bool lookingForRecoils = false;
  bool ignoreReplacedVacancy = false;
  string dummy;
  string energyStr, xStr, yStr, zStr, atomStr, vacancyStr, replacementStr;

  //---------------------------//
  //Load through collision file//
  //---------------------------//
  auto t_start = std::chrono::high_resolution_clock::now();
  while (getline(infile, line)) {
    //Clean weird ascii chars
    for (char& c : line) {
      if (!isprint(static_cast<unsigned char>(c))) {
          c = ' ';
      }
    }
    istringstream iss(line); 

    //If "==" or "--" in line, skip
    if (line.find("===") != string::npos || line.find("--") != string::npos) continue;     

    //elif "New Cascade" in line, start of a new primary--we need to store the energy and x,y,z positions
    else if (line.find("New Cascade") != string::npos) {
      iss >> dummy >> energyStr >> xStr >> yStr >> zStr;
      float eFloat = std::stof(energyStr);
      int32_t eRounded = static_cast<int32_t>(std::round(eFloat*1000.0f)); //convert to eV
      data.primaryEnergies.push_back(eRounded);
      data.primary_xLocs.push_back(stof(xStr)*0.1);
      data.primary_yLocs.push_back(stof(yStr)*0.1);
      data.primary_zLocs.push_back(stof(zStr)*0.1);
    }

    //elif "Prime Recoil" in line, start of a new cascade. Enable tracking.
    else if (line.find("Prime Recoil") != string::npos) {
      lookingForRecoils = true;
    }
    //elif "Summary" in line, end of a new cascade. Push back the cascade and clear vectors
    else if (line.find("Summary") != string::npos) {      
      lookingForRecoils = false;

      data.recoil_xLocs.push_back(tmpRecoil_xLocs);
      data.recoil_yLocs.push_back(tmpRecoil_yLocs);
      data.recoil_zLocs.push_back(tmpRecoil_zLocs);
      data.recoil_nVacs.push_back(tmpRecoil_nVacs);

      //Clean vectors of the previous cascade data
      tmpRecoil_xLocs.clear();
      tmpRecoil_yLocs.clear();
      tmpRecoil_zLocs.clear();      
      tmpRecoil_nVacs.clear();
    }

    //elif "For ion" in line--end of a primary. Process the throw and then clear vectors
    else if (line.find("For Ion") != string::npos) {
      if (throwNum%100==0) {
        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;

        cout << "Processed " << throwNum << " throws in "
            << elapsed.count() << " seconds." << endl;
        t_start = std::chrono::high_resolution_clock::now();
      }
      
      throwNum++; 
      ProcessThrow(data,trimTree,energy,xs,ys,zs,nVacs,clusteringDistance_nm,binSize_eV,maxEntriesPerBin,maxEnergy_eV,minEnergy_eV);

      //Clear primary info
      data.primaryEnergies.clear();
      data.primary_xLocs.clear();
      data.primary_yLocs.clear();
      data.primary_zLocs.clear();
      //Clear recoil info
      data.recoil_xLocs.clear();
      data.recoil_yLocs.clear();
      data.recoil_zLocs.clear();
      data.recoil_nVacs.clear();
    }

    //If enable recoil tracking, store the energy, xloc, yloc, zloc if nVacanies >0 and atom is in atomsToTrack
    if (lookingForRecoils) {
      iss >> dummy >> atomStr >> energyStr >> xStr >> yStr >> zStr >> vacancyStr >> replacementStr;
      int16_t nVacs=static_cast<int16_t>(stoi(vacancyStr));
      if (vacancyStr=="1") {
        if (ignoreReplacedVacancy==true) {
          ignoreReplacedVacancy = false;
          nVacs=0;
        }
        if (find(atomsToTrack.begin(), atomsToTrack.end(), atomStr) == atomsToTrack.end()) {
          nVacs=0;
        }
        tmpRecoil_xLocs.push_back(stof(xStr)*0.1);
        tmpRecoil_yLocs.push_back(stof(yStr)*0.1);
        tmpRecoil_zLocs.push_back(stof(zStr)*0.1);
        tmpRecoil_nVacs.push_back(nVacs);
      }
      else if (replacementStr=="1") {
        ignoreReplacedVacancy=true;
      }
    }
  }

  //Having some memory issues sorting, so we're going to try writing to disk, closing, re-opening, and sorting.
  trimTree->Write("trimTree",TObject::kOverwrite);
  outputFile->Close();

  // Reopen input file
  TFile* fin = new TFile(outputFilename.c_str(), "READ");
  TTree* tree = (TTree*)fin->Get("trimTree");
  tree->BuildIndex("energy_eV");

  // Sort and write to temp file
  TFile* fout = new TFile("trimTree_sorted.root", "RECREATE");
  TTree* sorted = tree->CopyTree("");
  sorted->SetBasketSize("*", 1000000);
  sorted->SetAutoFlush(0);
  sorted->Write("trimTree", TObject::kOverwrite);
  fout->Close();
  fin->Close();

  // Replace original file
  if (std::remove(outputFilename.c_str()) != 0)
      std::cerr << "Failed to remove old file\n";
  if (std::rename("trimTree_sorted.root", outputFilename.c_str()) != 0)
      std::cerr << "Failed to rename sorted file\n";

}



///////////////////
//THROW PROCESSOR//
///////////////////
void ProcessThrow(CascadeData& data, TTree* tree, int32_t& energy, vector<float>& xs, vector<float>& ys, vector<float>& zs, vector<int16_t>&nVacs,float clusteringDistance_nm, int32_t binSize_eV, int maxEntriesPerBin,int32_t maxEnergy_eV,int32_t minEnergy_eV) {
  float startX = 0;
  float startY = 0;
  float startZ = 0;

  float prevX = startX;
  float prevY = startY;
  float prevZ = startZ;

  for (size_t i = 0; i < data.primaryEnergies.size(); i++) {
    //Shouldn't be empty but we'll check
    if (data.recoil_xLocs[i].empty()) continue;

    //If energy is below max and above min, process
    if ((data.primaryEnergies[i]<=maxEnergy_eV)&&(data.primaryEnergies[i]>minEnergy_eV)) {

      //If we have too much data in the bin, skip
      int binIndex = static_cast<int>(data.primaryEnergies[i] / binSize_eV);
      if (energyBinCounter[binIndex] <= maxEntriesPerBin) {

        //1. Get direction of primary
        //Get direction primary took from previous 
        float dx = data.primary_xLocs[i] - prevX;
        float dy = data.primary_yLocs[i] - prevY;
        float dz = data.primary_zLocs[i] - prevZ;
        //unit normalize
        float norm = sqrt(dx*dx + dy*dy + dz*dz);
        if (norm == 0) continue;  //shouldn't happen
        dx /= norm;
        dy /= norm;
        dz /= norm;

        //2. Compute rotation matrix
        array<array<float, 3>, 3> rotMat = ComputeRotationMatrix(dx,dy,dz);

        //3. Calculate the offset to apply to all points
        float offsetX = data.primary_xLocs[i] + startX;
        float offsetY = data.primary_yLocs[i] + startY;
        float offsetZ = data.primary_zLocs[i] + startZ;

        //Flatten the data, apply offset
        vector<float> flatX;
        vector<float> flatY; 
        vector<float> flatZ;
        vector<int16_t> flatNVacs;
        for (size_t j = i; j < data.primaryEnergies.size(); j++ ) {
          for (size_t k = 0; k < data.recoil_xLocs[j].size(); k++) {
            //Flatten future data while applying the offset 
            flatX.push_back(data.recoil_xLocs[j][k]-offsetX);
            flatY.push_back(data.recoil_yLocs[j][k]-offsetY);
            flatZ.push_back(data.recoil_zLocs[j][k]-offsetZ);
            flatNVacs.push_back(data.recoil_nVacs[j][k]);
          }
        }

        //Rotate & cluster
        std::map<std::tuple<int, int, int>, int> clusters;   
        for (size_t j=0; j < flatX.size(); j++ ) {
          //Rotate all subsequent recoils such that they correspond to the incident particle traveling in the (1,0,0 direction)
          float x = rotMat[0][0]*flatX[j] + rotMat[0][1]*flatY[j] + rotMat[0][2]*flatZ[j];
          float y = rotMat[1][0]*flatX[j] + rotMat[1][1]*flatY[j] + rotMat[1][2]*flatZ[j];
          float z = rotMat[2][0]*flatX[j] + rotMat[2][1]*flatY[j] + rotMat[2][2]*flatZ[j];
        
          int ix = round(x / clusteringDistance_nm);
          int iy = round(y / clusteringDistance_nm);
          int iz = round(z / clusteringDistance_nm);
          auto key = make_tuple(ix, iy, iz);

          clusters[key] += flatNVacs[j];  // nVacs
        }

        //5. Clear branch vectors
        xs.clear();
        ys.clear();
        zs.clear();
        nVacs.clear();

        //6. Fill
        energy = data.primaryEnergies[i];
        for (const auto& [key, count] : clusters) {
          auto [ix, iy, iz] = key;
      
          // Convert back to physical lattice site
          float x = ix * clusteringDistance_nm;
          float y = iy * clusteringDistance_nm;
          float z = iz * clusteringDistance_nm;
      
          if (count>0) {
            xs.push_back(x);
            ys.push_back(y);
            zs.push_back(z);
            nVacs.push_back(count);
          }
        }
        if (xs.size()>0) {
          tree->Fill();
          energyBinCounter[binIndex]++;
        }
      }
    }

    //7. Update prevX, prevY, prevZ for the next iteration
    prevX = data.primary_xLocs[i];
    prevY = data.primary_yLocs[i];
    prevZ = data.primary_zLocs[i];
  }
}

//Uses Rodriguez formulat to calculate the rotation matrix from (dx,dy,dz) to (1,0,0) frame.
//We can then easily apply this to all recoil points
array<array<float, 3>, 3> ComputeRotationMatrix(float dx, float dy, float dz) {
  const float sinThetaCutoff = 0.000000174532; //sin(0.0000001 deg) — treat as aligned if below this

  //Normalize direction vector
  float norm = sqrt(dx*dx + dy*dy + dz*dz);
  if (norm == 0.0f) {
    //Degenerate case
    return {{
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.}
    }};
  }
  dx /= norm;
  dy /= norm;
  dz /= norm;

  //Rotation axis = (dx,dy,dz) × (1,0,0) = (0, dz, -dy)
  float ax = 0.f;
  float ay = dz;
  float az = -dy;

  //sin(theta) = |axis| (since both vectors are unit length)
  float sinTheta = sqrt(ay * ay + az * az);

  //If already aligned, skip rotation
  if (sinTheta < sinThetaCutoff) {
    return {{
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.}
    }};
  }

  //Normalize rotation axis
  ax /= sinTheta;
  ay /= sinTheta;
  az /= sinTheta;

  float cosTheta = dx; //cos(theta) = dot((dx,dy,dz), (1,0,0)) = dx
  float t = 1.0f - cosTheta;
  float s = sinTheta;

  //Rodrigues rotation matrix
  return {{
    {cosTheta + t*ax*ax,     t*ax*ay - s*az,    t*ax*az + s*ay},
    {t*ax*ay + s*az,         cosTheta + t*ay*ay, t*ay*az - s*ax},
    {t*ax*az - s*ay,         t*ay*az + s*ax,    cosTheta + t*az*az}
  }};
}
