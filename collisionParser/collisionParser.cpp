#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <map>
#include <chrono>
#include <cstdio> 
#include <filesystem>
#include <zlib.h>

#include "gzstream.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TFile.h"

using namespace std;

struct CascadeData {
  vector<int32_t> primaryEnergies; //Holds the primary particle enrgy
  vector<float> primary_xLocs;
  vector<float> primary_yLocs;
  vector<float> primary_zLocs;
  vector<vector<float>> recoil_xLocs;
  vector<vector<float>> recoil_yLocs;
  vector<vector<float>> recoil_zLocs;
  vector<vector<uint8_t>> recoil_nVacs;
  vector<vector<uint8_t>> recoil_displacedZs;
  vector<vector<int32_t>> recoil_energies_eV;
};


//Functions
void ProcessThrow(CascadeData& data, 
  TTree* tree, int32_t& energy, vector<float>& xs, vector<float>& ys, vector<float>& zs, vector<uint8_t>& nVacs, 
  vector<uint8_t>& recoil_displacedZs, vector<int32_t>& recoilEnergies_eV,vector<int16_t>& recoilNums,
  int32_t binSize_eV, int maxEntriesPerBin, int32_t maxEnergy_eV, int32_t minEnergy_eV);

//Takes dx,dy,dz unit vector and computes the matrix to rotate arbitrary points in that space to 1,0,0
array<array<float, 3>, 3> ComputeRotationMatrix(float dx, float dy, float dz);

std::map<int, int> energyBinCounter;

int main(int argc, char* argv[]) {

  //Srim settings  
  int32_t minEnergy_eV = 0;
  int32_t maxEnergy_eV = 200e3;  //We start our SRIM sim slightly above this, at least 100 keV to allow for "burn in"

  int maxEntriesPerBin = 100;     //Avoid filling up too much data at lower energys
  int32_t binSize_eV = 1e3;       //Bin size, in eV. 1 keV right now
 
  //--------------------//
  //Parse cmd line input//
  //--------------------//
  if ((argc != 2)&&(argc != 3)) {
      cout << "Usage: " << argv[0] << " <COLLISON.txt> <optional output filename>" << endl;
      return 1;
  }

  //Input filename
  string inputFilename = argv[1];

  //Make output filename 
  std::string outputFilename;
  if (argc == 3) {
    outputFilename = argv[2];
  }
  else {
    std::filesystem::path p(inputFilename);
    std::string prefix = p.stem().string(); 
    outputFilename = prefix + ".root";
  }

  std::cout << "Input file: " << inputFilename << std::endl;
  std::cout << "Output file: " << outputFilename << std::endl;

  //--------------------//
  //Make the output file//
  //--------------------//
  TFile* outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  TTree* trimTree = new TTree("trimTree", "Clustered recoil tracks");
  
  //Incident ion
  int32_t ionEnergy_eV = 0; 
  //Values per site
  vector<float> xs;
  vector<float> ys;
  vector<float> zs;
  vector<uint8_t> nVacs;
  vector<uint8_t> displacedZs;
  vector<int32_t> recoilEnergies_eV;
  vector<int16_t> recoilNums;

  trimTree->Branch("ionEnergy_eV", &ionEnergy_eV, "ionEnergy_eV/i");
  trimTree->Branch("recoilNums", &recoilNums);
  trimTree->Branch("xs_nm", &xs);
  trimTree->Branch("ys_nm", &ys);
  trimTree->Branch("zs_nm", &zs);
  trimTree->Branch("nVacs", &nVacs);
  trimTree->Branch("displacedAtoms_Z", &displacedZs);
  trimTree->Branch("recoilEnergies_eV", &recoilEnergies_eV);
  
  trimTree->SetMaxTreeSize(200000000000LL); //200GB
  trimTree->SetAutoFlush(10000);
  trimTree->SetAutoSave(500000000);  // Every 500 MB
  //trimTree->SetBasketSize("*", 1000000); //trying to debug...

  //-------------//
  //Open the file//
  //-------------// 
  //Support for either zipped or txt files
  std::ifstream txtfile;
  igzstream gzfile;      // Declare gzstream object
  std::istream* infile = nullptr;

  if (inputFilename.size() >= 3 && inputFilename.substr(inputFilename.size() - 3) == ".gz") {
    gzfile.open(inputFilename.c_str(), std::ios::in);
    if (!gzfile) {
      std::cerr << "Error opening gz file: " << inputFilename << std::endl;
      return 1;
    }
    infile = &gzfile;
  } else {
    txtfile.open(inputFilename);
    if (!txtfile) {
      std::cerr << "Error opening text file: " << inputFilename << std::endl;
      return 1;
    }
    infile = &txtfile;
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
  vector<uint8_t> tmpRecoil_nVacs;
  vector<uint8_t> tmpRecoil_displacedZs;
  vector<int32_t> tmpRecoil_energies_eV;

  int throwNum = 0; //Only for tracking how far in the processing we are
  bool lookingForRecoils = false;
  bool ignoreReplacedVacancy = false;
  string dummy;
  string energyStr, recoilEnergyStr, xStr, yStr, zStr, atomStr, vacancyStr, replacementStr;
  
  double initialThrownEnergy_eV = 0;
  double thrownEnergy_eV = 0;

  //---------------------------//
  //Load through collision file//
  //---------------------------//
  auto t_start = std::chrono::high_resolution_clock::now();
  while (std::getline(*infile, line)) {
    //Clean weird ascii chars
    for (char& c : line) {
      if (!isprint(static_cast<unsigned char>(c))) {
          c = ' ';
      }
    }
    istringstream iss(line); 

    //Parse initial energy line
    if (line.find("Ion Energy =") != string::npos) {
        std::string digits;
        for (char c : line) {
            if ((c >= '0' && c <= '9') || c == '.') {
                digits.push_back(c);
            }
        }
        if (!digits.empty()) {
            float E0_keV = std::stof(digits);
            initialThrownEnergy_eV = static_cast<int32_t>(std::round(E0_keV * 1000.0f));
            thrownEnergy_eV = initialThrownEnergy_eV;
        }
        continue;
    }

    //If "==" or "--" in line, skip
    if (line.find("===") != string::npos || line.find("--") != string::npos) continue;     

    //elif "New Cascade" in line, start of a new primary--we need to store the energy and x,y,z positions
    else if (line.find("New Cascade") != string::npos) {
      iss >> dummy >> energyStr >> xStr >> yStr >> zStr >> dummy >> dummy >> recoilEnergyStr;
      
      float preCollisionEnergy_keV = std::stof(energyStr);
      int32_t preCollisionEnergy_eV = static_cast<int32_t>(std::round(preCollisionEnergy_keV*1000.0f)); //convert to eV

      int32_t collisionRecoilEnergy_eV = static_cast<int32_t>(std::round(std::stof(recoilEnergyStr)));

      data.primaryEnergies.push_back(thrownEnergy_eV);
      data.primary_xLocs.push_back(stof(xStr)*0.1);
      data.primary_yLocs.push_back(stof(yStr)*0.1);
      data.primary_zLocs.push_back(stof(zStr)*0.1);

      thrownEnergy_eV = preCollisionEnergy_eV - collisionRecoilEnergy_eV;
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
      data.recoil_displacedZs.push_back(tmpRecoil_displacedZs);
      data.recoil_energies_eV.push_back(tmpRecoil_energies_eV);

      //Clean vectors of the previous cascade data
      tmpRecoil_xLocs.clear();
      tmpRecoil_yLocs.clear();
      tmpRecoil_zLocs.clear();      
      tmpRecoil_nVacs.clear();
      tmpRecoil_displacedZs.clear();
      tmpRecoil_energies_eV.clear();
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
      ProcessThrow(data, //Contains data to load
        trimTree,ionEnergy_eV,xs,ys,zs,nVacs,displacedZs,recoilEnergies_eV,recoilNums, //Branches to fill
        binSize_eV,maxEntriesPerBin,maxEnergy_eV,minEnergy_eV); //Flags

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
      data.recoil_displacedZs.clear();
      data.recoil_energies_eV.clear();
      //Reset exit energy
      thrownEnergy_eV = initialThrownEnergy_eV;
    }

    //If enable recoil tracking, store the energy, xloc, yloc, zloc if nVacanies >0 and atom is in atomsToTrack
    if (lookingForRecoils) {
      iss >> dummy >> atomStr >> energyStr >> xStr >> yStr >> zStr >> vacancyStr >> replacementStr;
      uint8_t nVacs=static_cast<int8_t>(std::stoul(vacancyStr));
      if (vacancyStr=="1") {
        if (ignoreReplacedVacancy==true) {
          ignoreReplacedVacancy = false;
          nVacs=0;
        }
        tmpRecoil_xLocs.push_back(stof(xStr)*0.1);
        tmpRecoil_yLocs.push_back(stof(yStr)*0.1);
        tmpRecoil_zLocs.push_back(stof(zStr)*0.1);
        tmpRecoil_nVacs.push_back(nVacs);
        tmpRecoil_displacedZs.push_back(static_cast<uint8_t>(std::stoul(atomStr)));
        tmpRecoil_energies_eV.push_back(static_cast<int32_t>(std::round(stof(energyStr))));
      }
      else if (replacementStr=="1") {
        ignoreReplacedVacancy=true;
      }
    }
  }

  //Write
  trimTree->Write("trimTree",TObject::kOverwrite);
  outputFile->Close();
}

///////////////////
//THROW PROCESSOR//
///////////////////
void ProcessThrow(CascadeData& data, 
  TTree* tree, int32_t& energy, vector<float>& xs, vector<float>& ys, vector<float>& zs, vector<uint8_t>&nVacs, 
  vector<uint8_t>& displacedZs, vector<int32_t>& recoilEnergies_eV,vector<int16_t>& recoilNums,
  int32_t binSize_eV, int maxEntriesPerBin,int32_t maxEnergy_eV,int32_t minEnergy_eV) {

  float startX = 0;
  float startY = 0;
  float startZ = 0;

  float prevX = startX;
  float prevY = startY;
  float prevZ = startZ;

  for (size_t i = 0; i < data.primaryEnergies.size(); i++) {

    if (data.recoil_xLocs[i].empty()) {
      std::cout<<"Found recoil without any vacancies"<<std::endl;
    }
    
    //If energy is below max and above min, process
    if ((data.primaryEnergies[i]<=maxEnergy_eV)&&(data.primaryEnergies[i]>minEnergy_eV)) {

      //If we have too much data in the bin, skip
      int binIndex = static_cast<int>(data.primaryEnergies[i] / binSize_eV);
      if (energyBinCounter[binIndex] <= maxEntriesPerBin) {

        // If no vacancies, still record an empty entry
        if (data.recoil_xLocs[i].empty()) {
            xs.clear();
            ys.clear();
            zs.clear();
            nVacs.clear();
            displacedZs.clear();
            recoilEnergies_eV.clear();
            recoilNums.clear();
            energy = data.primaryEnergies[i];
            tree->Fill();
            energyBinCounter[binIndex]++;
            continue;
        }

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
        float offsetX = prevX;
        float offsetY = prevY;
        float offsetZ = prevZ;

        //Flatten the data, apply offset to all downstream depositions
        vector<float> flatX;
        vector<float> flatY; 
        vector<float> flatZ;
        vector<uint8_t> flatNVacs;
        vector<uint8_t> flatDisplacedZs;
        vector<float> flatRecoilEnergies_eV;
        vector<int16_t> flatRecoilNums;
        for (size_t j = i; j < data.primaryEnergies.size(); j++ ) {
          for (size_t k = 0; k < data.recoil_xLocs[j].size(); k++) {
            //Flatten future data while applying the offset 
            flatX.push_back(data.recoil_xLocs[j][k]-offsetX);
            flatY.push_back(data.recoil_yLocs[j][k]-offsetY);
            flatZ.push_back(data.recoil_zLocs[j][k]-offsetZ);
            flatNVacs.push_back(data.recoil_nVacs[j][k]);
            flatDisplacedZs.push_back(data.recoil_displacedZs[j][k]);
            flatRecoilEnergies_eV.push_back(data.recoil_energies_eV[j][k]);
            flatRecoilNums.push_back(static_cast<int16_t>(j-i));
          }
        }

        //Clear branch vectors
        xs.clear();
        ys.clear();
        zs.clear();
        nVacs.clear();
        displacedZs.clear();
        recoilEnergies_eV.clear();
        recoilNums.clear();

        //Rotate & Fill
        std::map<std::tuple<int, int, int>, int> clusters;   
        for (size_t j=0; j < flatX.size(); j++ ) {
          //Rotate all subsequent recoils such that they correspond to the incident particle traveling in the (1,0,0 direction)
          float x = rotMat[0][0]*flatX[j] + rotMat[0][1]*flatY[j] + rotMat[0][2]*flatZ[j];
          float y = rotMat[1][0]*flatX[j] + rotMat[1][1]*flatY[j] + rotMat[1][2]*flatZ[j];
          float z = rotMat[2][0]*flatX[j] + rotMat[2][1]*flatY[j] + rotMat[2][2]*flatZ[j];
        
          xs.push_back(x);
          ys.push_back(y);
          zs.push_back(z);
          nVacs.push_back(flatNVacs.at(j));
          displacedZs.push_back(flatDisplacedZs.at(j));
          recoilEnergies_eV.push_back(flatRecoilEnergies_eV.at(j));
          recoilNums.push_back(flatRecoilNums.at(j));
        }

        //6. Fill
        energy = data.primaryEnergies[i];
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

//Uses Rodriguez formula to calculate the rotation matrix from (dx,dy,dz) to (1,0,0) frame.
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
