#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <chrono>

#include "TTree.h"
#include "TTreeIndex.h"
#include "TFile.h"
//#define WRITE_CSV

using namespace std;

struct CascadeData {
  vector<float> primaryEnergies;
  vector<float> primary_xLocs;
  vector<float> primary_yLocs;
  vector<float> primary_zLocs;
  vector<vector<float>> recoilEnergies;
  vector<vector<float>> recoil_xLocs;
  vector<vector<float>> recoil_yLocs;
  vector<vector<float>> recoil_zLocs;
};

//Functions
void ProcessThrow(CascadeData& data, TTree* tree, float& energy, vector<float>& xs, vector<float>& ys, vector<float>& zs, vector<int>& nVacs, float startOffset_nm, float clusteringDistance_nm, float binSize_keV, int maxEntriesPerBin, float maxEnergy_keV);

//Takes dx,dy,dz unit vector and computes the matrix to rotate arbitrary points in that space to 1,0,0
array<array<float, 3>, 3> ComputeRotationMatrix(float dx, float dy, float dz);

std::map<int, int> energyBinCounter;

int main(int argc, char* argv[]) {

  //Srim settings
  float startOffset_nm = 1.0;  //From SRIM input, 10 A = 1nm, we start partway in to allow backscattered particles to show
  float maxEnergy_keV = 5000;  //We start our SRIM sim slightly above this, at least 100 keV to allow for "burn in"

  int maxEntriesPerBin = 5e4;  //Avoid filling up too much data at lower energys
  float binSize_keV = 5;       //Bin size
  
  //
  float clusteringDistance_nm = 4.026*0.1*0.5;  //4.026 Angstroms is the lattice distance, but we're using 1.5 lattice spacings
  vector<string> atomsToTrack = {"19"};         //Don't care about Li vacancies, not optical
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
  //outputFile->SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kLZ4);
  //outputFile->SetCompressionLevel(4);  
  TTree* unsortedTree = new TTree("unsortedTree", "Clustered recoil tracks");
  unsortedTree->SetDirectory(0);
  
  float energy = 0.0;
  vector<float> xs;
  vector<float> ys;
  vector<float> zs;
  vector<int> nVacs;

  unsortedTree->Branch("energy_keV", &energy, "energy_keV/F");
  unsortedTree->Branch("xs_nm", &xs);
  unsortedTree->Branch("ys_nm", &ys);
  unsortedTree->Branch("zs_nm", &zs);
  unsortedTree->Branch("nVacs", &nVacs);
  
  TTree* trimTree = (TTree*)unsortedTree->CloneTree(0);
  trimTree->SetName("trimTree" );
  trimTree->SetDirectory(outputFile);

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
  vector<float> tmpRecoilEnergies;
  vector<float> tmpRecoil_xLocs;
  vector<float> tmpRecoil_yLocs;
  vector<float> tmpRecoil_zLocs;

  int throwNum = 0; //Only for tracking how far in the processing we are
  bool lookingForRecoils = false;
  string dummy;
  string energyStr, xStr, yStr, zStr, atomStr, vacancyStr;

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

    //elif "New Cascade" in line, then we need to store the energy and x,y,z positions
    else if (line.find("New Cascade") != string::npos) {
      iss >> dummy >> energyStr >> xStr >> yStr >> zStr;
      data.primaryEnergies.push_back(stof(energyStr));
      data.primary_xLocs.push_back(stof(xStr)*0.1);
      data.primary_yLocs.push_back(stof(yStr)*0.1);
      data.primary_zLocs.push_back(stof(zStr)*0.1);
    }
    //elif "Prime Recoil" in line, enable tracking. If "Summary" in line, disable tracking and process the event
    else if (line.find("Prime Recoil") != string::npos) {
      lookingForRecoils = true;
    }
    else if (line.find("Summary") != string::npos) {

      lookingForRecoils = false;

      data.recoilEnergies.push_back(tmpRecoilEnergies);
      data.recoil_xLocs.push_back(tmpRecoil_xLocs);
      data.recoil_yLocs.push_back(tmpRecoil_yLocs);
      data.recoil_zLocs.push_back(tmpRecoil_zLocs);

      //Clean vectors of the previous cascade data
      tmpRecoilEnergies.clear();
      tmpRecoil_xLocs.clear();
      tmpRecoil_yLocs.clear();
      tmpRecoil_zLocs.clear();
      
    }
    else if (line.find("For Ion") != string::npos) {
      if (throwNum%100==0) {
        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;

        cout << "Processed " << throwNum << " throws in "
            << elapsed.count() << " seconds." << endl;
        t_start = std::chrono::high_resolution_clock::now();

        if (throwNum>500) {
          break;
        }
      }
      
      throwNum++; 
      ProcessThrow(data,unsortedTree,energy,xs,ys,zs,nVacs,startOffset_nm,clusteringDistance_nm,binSize_keV,maxEntriesPerBin,maxEnergy_keV);

      //Clear primary info
      data.primaryEnergies.clear();
      data.primary_xLocs.clear();
      data.primary_yLocs.clear();
      data.primary_zLocs.clear();
      //Clear recoil ingo
      data.recoilEnergies.clear();
      data.recoil_xLocs.clear();
      data.recoil_yLocs.clear();
      data.recoil_zLocs.clear();
    }

    //If enable recoil tracking, store the energy, xloc, yloc, zloc if nVacanies >0 and atom is in atomsToTrack
    if (lookingForRecoils) {
      iss >> dummy >> atomStr >> energyStr >> xStr >> yStr >> zStr >> vacancyStr;
      if ((vacancyStr=="1") && (find(atomsToTrack.begin(), atomsToTrack.end(), atomStr) != atomsToTrack.end())) {

        tmpRecoilEnergies.push_back(stof(energyStr));
        tmpRecoil_xLocs.push_back(stof(xStr)*0.1);
        tmpRecoil_yLocs.push_back(stof(yStr)*0.1);
        tmpRecoil_zLocs.push_back(stof(zStr)*0.1);
      }
    }
  }

  //Sort
  cout<<"Sorting tree...";
  unsortedTree->BuildIndex("energy_keV");
  TTreeIndex* treeIndex = (TTreeIndex*)unsortedTree->GetTreeIndex();
  for( size_t i = 0; i < treeIndex->GetN(); i++ ) {
    unsortedTree->GetEntry(treeIndex->GetIndex()[i]);
    trimTree->Fill();      
  }
  cout<<"sorted!"<<endl;

  //Write CSV
  #ifdef WRITE_CSV
  cout<<"Writing CSV..."<<endl;
  ofstream csvOut("trimTracks.csv");
  csvOut << "energy_keV,xs_nm,ys_nm,zs_nm,nVacs\n";

  for (Long64_t i = 0; i < trimTree->GetEntries(); i++) {
      trimTree->GetEntry(i);
      for (size_t j = 0; j < xs.size(); j++) {
          csvOut<<energy<< ","<<xs[j]<<","<<ys[j]<< ","<<zs[j]<<","<<(int)nVacs[j]<<"\n";
      }
  }
  csvOut.close();
  cout<<"Done!"<<endl;
  #endif

  //Write root file
  trimTree->Write("trimTree",TObject::kOverwrite);
  unsortedTree->SetDirectory(0);
  delete unsortedTree;
  outputFile->Close();
}


void ProcessThrow(CascadeData& data, TTree* tree, float& energy, vector<float>& xs, vector<float>& ys, vector<float>& zs, vector<int>&nVacs, float startOffset_nm, float clusteringDistance_nm, float binSize_keV, int maxEntriesPerBin,float maxEnergy_keV) {
  float startX = startOffset_nm;
  float startY = 0;
  float startZ = 0;

  float prevX = startX;
  float prevY = startY;
  float prevZ = startZ;

  for (size_t i = 0; i < data.primaryEnergies.size() - 1; i++) {
    //Shouldn't be empty but we'll check
    if (data.recoil_xLocs[i].empty()) continue;

    //If energy is above our max, skip
    if (data.primaryEnergies[i]<=maxEnergy_keV) {

      //If we have too much data in the bin, skip
      int binIndex = static_cast<int>(data.primaryEnergies[i] / binSize_keV);
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
        for (size_t j = i; j < data.primaryEnergies.size(); j++ ) {
          for (size_t k = 0; k < data.recoil_xLocs[j].size(); k++) {
            //Flatten future data while applying the offset 
            flatX.push_back(data.recoil_xLocs[j][k]-offsetX);
            flatY.push_back(data.recoil_yLocs[j][k]-offsetY);
            flatZ.push_back(data.recoil_zLocs[j][k]-offsetZ);
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
          clusters[key] += 1;
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
      
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
          nVacs.push_back(count);
        }
        tree->Fill();
        energyBinCounter[binIndex]++;
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
  //If
  const float sinThetaCutoff = 0.00000174532; //sin(0.00001 degrees). If a vector is within this of the x-axis, we will skip the rotation.

  //Compute rotation axis: (dx,dy,dz) x (1,0,0)
  float ax = 0; //dy * uz - dz * uy;
  float ay = dz; //dz * ux - dx * uz;
  float az = -dy; //dx * uy - dy * ux;

  // Magnitude of cross product = |a||b|sin(theta) = sin(theta)
  float sinTheta = sqrt(ay * ay + az * az); //sqrt(ax * ax + ay * ay + az * az)

  // If vectors are nearly aligned, return identity
  if (sinTheta < sinThetaCutoff) {
    return {{
        {1., 0., 0.},
        {0., 1., 0.},
        {0., 0., 1.}
    }};
}
  // Normalize rotation axisaxis
  ax /= sinTheta;
  ay /= sinTheta;
  az /= sinTheta;

  //(dx,dy,dz) dot (1,0,0) = |a||b|cos(theta) = cos(theta)
  float cosTheta = dx; //vx * ux + vy * uy + vz * uz;

  // Rodrigues formula
  float t = 1 - cosTheta;

  return {{
      {cosTheta, -sinTheta*az, sinTheta*ay},      //{t*ax*ax + c    , t*ax*ay - s*az, t*ax*az + s*ay},
      {sinTheta*az, t*ay*ay + cosTheta, t*ay*az}, //{t*ax*ay + s*az, t*ay*ay + c    , t*ay*az - s*ax},
      {-sinTheta*ay, t*ay*az, t*az*az + cosTheta} //{t*ax*az - s*ay, t*ay*az + s*ax, t*az*az + c    }
  }};
}