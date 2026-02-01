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
  vector<int32_t> primaryEnergies;
  vector<float> xLocs;
  vector<float> yLocs;
  vector<float> zLocs;
  vector<float> nVacs;
  vector<uint8_t> displacedZs;
  vector<float> recoilEnergies;
};

//Functions
void ProcessThrow(CascadeData& data, 
  TTree* tree, int32_t& energy, vector<float>& xs, vector<float>& ys, vector<float>& zs, vector<float>& nVacs, 
  vector<uint8_t>& recoil_displacedZs, vector<float>& recoilEnergies_eV,vector<int16_t>& recoilNums,
  int32_t binSize_eV, int maxEntriesPerBin, int32_t maxEnergy_eV, int32_t minEnergy_eV);
int ElementSymbolToZ(const std::string& s);

//Takes dx,dy,dz unit vector and computes the matrix to rotate arbitrary points in that space to 1,0,0
array<array<float, 3>, 3> ComputeRotationMatrix(float dx, float dy, float dz);

std::map<int, int> energyBinCounter;



int main(int argc, char* argv[]) {

  //Srim settings  
  int32_t minEnergy_eV = 0;
  int32_t maxEnergy_eV = 300e3;  //We start our SRIM sim slightly above this, at least 1% or ~100 keV to allow for "burn in"

  int maxEntriesPerBin = 1;     //Avoid filling up too much data at lower energys
  int32_t binSize_eV = 1;       //Bin size, in eV
 
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
  vector<float> nVacs;
  vector<uint8_t> displacedZs;
  vector<float> recoilEnergies_eV;
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

  int throwNum = 0; //Only for tracking how far in the processing we are
  bool lookingForRecoils = false;
  string dummy;
  string energyStr, recoilEnergyStr, xStr, yStr, zStr, atomStr, vacancyStr;
  
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
            std::cout<<"Input energy is "<<thrownEnergy_eV<<" eV"<<std::endl;
        }
        continue;
    }
    //If "==" or "--" in line, skip
    else if (line.find("===") != string::npos || line.find("--") != string::npos) {
      continue;     
    }
    //elif "REPLAC INTER" in line, start of a new primary--we need to store the energy and x,y,z positions
    else if (line.find("REPLAC INTER") != string::npos) {
      lookingForRecoils = true;
    }
    else if (line.find("For Ion") != string::npos) {
      lookingForRecoils = false; 

      //Print
      if (throwNum%100==0) {
        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = t_end - t_start;

        cout << "Processed " << throwNum << " throws in "
            << elapsed.count() << " seconds." << endl;
        t_start = std::chrono::high_resolution_clock::now();
      }
      throwNum++; 

      //Process data
      ProcessThrow(data, 
        trimTree,ionEnergy_eV,xs,ys,zs,nVacs,displacedZs,recoilEnergies_eV,recoilNums, //Branches to fill
        binSize_eV,maxEntriesPerBin,maxEnergy_eV,minEnergy_eV); //Flags

      //Clear vectors
      data.primaryEnergies.clear();
      data.xLocs.clear();
      data.yLocs.clear();
      data.zLocs.clear();
      data.displacedZs.clear();
      data.nVacs.clear();
      data.recoilEnergies.clear();
      
      //Reset exit energy
      thrownEnergy_eV = initialThrownEnergy_eV;
    }
    else if (lookingForRecoils==true) {
      //Parse string, load into vectors
      iss >> dummy >> energyStr >> xStr >> yStr >> zStr >> dummy >> atomStr >> recoilEnergyStr >> vacancyStr;
      data.primaryEnergies.push_back(thrownEnergy_eV);
      data.xLocs.push_back(stof(xStr)*0.1);
      data.yLocs.push_back(stof(yStr)*0.1);
      data.zLocs.push_back(stof(zStr)*0.1);
      data.nVacs.push_back(stof(vacancyStr)*0.1);
      data.displacedZs.push_back(ElementSymbolToZ(atomStr)); 
      data.recoilEnergies.push_back(stof(recoilEnergyStr));
      
      //Calculate the energy of this ion before the collision--we will use this as the throw energy for a subsequent primary
      float preCollisionEnergy_keV = std::stof(energyStr);
      float preCollisionEnergy_eV = preCollisionEnergy_keV*1000.0f;
      float collisionRecoilEnergy_eV = static_cast<int32_t>(std::round(std::stof(recoilEnergyStr)));
      thrownEnergy_eV = preCollisionEnergy_eV - collisionRecoilEnergy_eV;
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
  TTree* tree, int32_t& energy, vector<float>& xs, vector<float>& ys, vector<float>& zs, vector<float>&nVacs, 
  vector<uint8_t>& displacedZs, vector<float>& recoilEnergies_eV,vector<int16_t>& recoilNums,
  int32_t binSize_eV, int maxEntriesPerBin,int32_t maxEnergy_eV,int32_t minEnergy_eV) {

  float startX = 0;
  float startY = 0;
  float startZ = 0;

  float prevX = startX;
  float prevY = startY;
  float prevZ = startZ;
  
  for (size_t i = 0; i < data.primaryEnergies.size(); i++) {

    //Should only trigger if we throw sub-damage threshold primaries
    if (data.xLocs.empty()) {
      std::cout<<"Found recoil without any vacancies"<<std::endl;
    }
    
    //If energy is below max and above min, process
    if ((data.primaryEnergies[i]<=maxEnergy_eV)&&(data.primaryEnergies[i]>minEnergy_eV)) {
      
      //If we have too much data in the bin, skip
      int binIndex = static_cast<int>(data.primaryEnergies[i] / binSize_eV);
      if (energyBinCounter[binIndex] < maxEntriesPerBin) {

        // If no vacancies, still record an empty entry
        if (data.xLocs.empty()) {
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
        float dx = data.xLocs[i] - prevX;
        float dy = data.yLocs[i] - prevY;
        float dz = data.zLocs[i] - prevZ;
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
        vector<float> flatNVacs;
        vector<uint8_t> flatDisplacedZs;
        vector<float> flatRecoilEnergies_eV;
        vector<int16_t> flatRecoilNums;
        for (size_t j = i; j < data.primaryEnergies.size(); j++ ) {
          //Flatten future data while applying the offset 
          flatX.push_back(data.xLocs[j]-offsetX);
          flatY.push_back(data.yLocs[j]-offsetY);
          flatZ.push_back(data.zLocs[j]-offsetZ);
          flatNVacs.push_back(data.nVacs[j]);
          flatDisplacedZs.push_back(data.displacedZs[j]);
          flatRecoilEnergies_eV.push_back(data.recoilEnergies[j]);
          flatRecoilNums.push_back(static_cast<int16_t>(j-i));
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
    prevX = data.xLocs[i];
    prevY = data.yLocs[i];
    prevZ = data.zLocs[i];
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

int ElementSymbolToZ(const std::string& s) {
  if (s == "H")  return 1;
  if (s == "He") return 2;
  if (s == "Li") return 3;
  if (s == "Be") return 4;
  if (s == "B")  return 5;
  if (s == "C")  return 6;
  if (s == "N")  return 7;
  if (s == "O")  return 8;
  if (s == "F")  return 9;
  if (s == "Ne") return 10;
  if (s == "Na") return 11;
  if (s == "Mg") return 12;
  if (s == "Al") return 13;
  if (s == "Si") return 14;
  if (s == "P")  return 15;
  if (s == "S")  return 16;
  if (s == "Cl") return 17;
  if (s == "Ar") return 18;
  if (s == "K")  return 19;
  if (s == "Ca") return 20;
  if (s == "Sc") return 21;
  if (s == "Ti") return 22;
  if (s == "V")  return 23;
  if (s == "Cr") return 24;
  if (s == "Mn") return 25;
  if (s == "Fe") return 26;
  if (s == "Co") return 27;
  if (s == "Ni") return 28;
  if (s == "Cu") return 29;
  if (s == "Zn") return 30;
  if (s == "Ga") return 31;
  if (s == "Ge") return 32;
  if (s == "As") return 33;
  if (s == "Se") return 34;
  if (s == "Br") return 35;
  if (s == "Kr") return 36;
  if (s == "Rb") return 37;
  if (s == "Sr") return 38;
  if (s == "Y")  return 39;
  if (s == "Zr") return 40;
  if (s == "Nb") return 41;
  if (s == "Mo") return 42;
  if (s == "Tc") return 43;
  if (s == "Ru") return 44;
  if (s == "Rh") return 45;
  if (s == "Pd") return 46;
  if (s == "Ag") return 47;
  if (s == "Cd") return 48;
  if (s == "In") return 49;
  if (s == "Sn") return 50;
  if (s == "Sb") return 51;
  if (s == "Te") return 52;
  if (s == "I")  return 53;
  if (s == "Xe") return 54;
  if (s == "Cs") return 55;
  if (s == "Ba") return 56;
  if (s == "La") return 57;
  if (s == "Ce") return 58;
  if (s == "Pr") return 59;
  if (s == "Nd") return 60;
  if (s == "Pm") return 61;
  if (s == "Sm") return 62;
  if (s == "Eu") return 63;
  if (s == "Gd") return 64;
  if (s == "Tb") return 65;
  if (s == "Dy") return 66;
  if (s == "Ho") return 67;
  if (s == "Er") return 68;
  if (s == "Tm") return 69;
  if (s == "Yb") return 70;
  if (s == "Lu") return 71;
  if (s == "Hf") return 72;
  if (s == "Ta") return 73;
  if (s == "W")  return 74;
  if (s == "Re") return 75;
  if (s == "Os") return 76;
  if (s == "Ir") return 77;
  if (s == "Pt") return 78;
  if (s == "Au") return 79;
  if (s == "Hg") return 80;
  if (s == "Tl") return 81;
  if (s == "Pb") return 82;
  if (s == "Bi") return 83;
  if (s == "Po") return 84;
  if (s == "At") return 85;
  if (s == "Rn") return 86;
  if (s == "Fr") return 87;
  if (s == "Ra") return 88;
  if (s == "Ac") return 89;
  if (s == "Th") return 90;
  if (s == "Pa") return 91;
  if (s == "U")  return 92;
  if (s == "Np") return 93;
  if (s == "Pu") return 94;
  if (s == "Am") return 95;
  if (s == "Cm") return 96;
  if (s == "Bk") return 97;
  if (s == "Cf") return 98;
  if (s == "Es") return 99;
  if (s == "Fm") return 100;

  return -1;
}
