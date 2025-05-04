#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include "TTree.h"

using namespace std;

struct CascadeData {
  std::vector<float> primaryEnergies;
  std::vector<float> primary_xLocs;
  std::vector<float> primary_yLocs;
  std::vector<float> primary_zLocs;
  std::vector<std::vector<int>> recoilIons;
  std::vector<std::vector<float>> recoilEnergies;
  std::vector<std::vector<float>> recoil_xLocs;
  std::vector<std::vector<float>> recoil_yLocs;
  std::vector<std::vector<float>> recoil_zLocs;
};

//Global variables
//For tree
float energy = 0.0;
std::vector<float>* xs = new std::vector<float>();
std::vector<float>* ys = new std::vector<float>();
std::vector<float>* zs = new std::vector<float>();
std::vector<int>* nVac = new std::vector<int>();

//Hard-cded srim settings
vector<string> atomsToTrack = {"19"};
float startOffset_nm = 1.0; //From SRIM input, 10 A = 1nm, we start partway in to allow backscattered particles to show
//Hard coded globa variables
float clusteringDistance_nm = 4.026*0.1*1.5; //4.026 Angstroms is the lattice distance, but we're using 1.5 lattice spacings

//Functions
void ProcessThrow(CascadeData& data, TTree* tree);

//Takes dx,dy,dz unit vector and computes the matrix to rotate arbitrary points in that space to 1,0,0
std::array<std::array<float, 3>, 3> ComputeRotationMatrix(float dx, float dy, float dz);

int main(int argc, char* argv[]) {

  //--------------------//
  //Parse cmd line input//
  //--------------------//
  if (argc != 3) {
      cout << "Usage: " << argv[0] << " <COLLISON.txt> <ionName>" << endl;
      return 1;
  }

  string filename = argv[1];
  string ionName = argv[2];

  cout << "Input file: " << filename << endl;
  cout << "Ion species: " << ionName << endl;

  //--------------------//
  //Make the output file//
  //--------------------//
  TFile* outputFile = new TFile("trimTracks.root", "RECREATE");
  TTree* recoilTree = new TTree("trimTree", "Clustered recoil tracks");

  recoilTree->Branch("energy", &energy, "energy/F");
  recoilTree->Branch("xs", &xs);
  recoilTree->Branch("ys", &ys);
  recoilTree->Branch("zs", &zs);
  recoilTree->Branch("nVac", &nVac);


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
  vector<int> tmpRecoilIons;
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
      
      throwNum++; 
      iss >> dummy >> energyStr >> xStr >> yStr >> zStr;
      data.primaryEnergies.push_back(stof(energyStr));
      data.primary_xLocs.push_back(stof(xStr)*0.1);
      data.primary_yLocs.push_back(stof(yStr)*0.1);
      data.primary_zLocs.push_back(stof(zStr)*0.1);
    }
    //elif "Prime Recoil" in line, enable tracking. If "Summary" in line, disable tracking and process the event
    else if (line.find("Prime Recoil") != string::npos) {
      lookingForRecoils = true;
      tmpRecoilIons.clear();
      tmpRecoilEnergies.clear();
      tmpRecoil_xLocs.clear();
      tmpRecoil_yLocs.clear();
      tmpRecoil_zLocs.clear();
    }
    else if (line.find("Summary") != string::npos) {
      lookingForRecoils = false;

      data.recoilIons.push_back(tmpRecoilIons);
      data.recoilEnergies.push_back(tmpRecoilEnergies);
      data.recoil_xLocs.push_back(tmpRecoil_xLocs);
      data.recoil_yLocs.push_back(tmpRecoil_yLocs);
      data.recoil_zLocs.push_back(tmpRecoil_zLocs);
    }
    else if (line.find("For ion") != string::npos) {
      ProcessThrow(data,tree);

      //Clear cascade data
      data.primaryEnergies.clear();
      data.primary_xLocs.clear();
      data.primary_yLocs.clear();
      data.primary_zLocs.clear();
      data.recoilIons.clear();
      data.recoilEnergies.clear();
      data.recoil_xLocs.clear();
      data.recoil_yLocs.clear();
      data.recoil_zLocs.clear();
    }

    //If enable recoil tracking, store the energy, xloc, yloc, zloc if nVacanies >0 and atom is in atomsToTrack
    if (lookingForRecoils) {
      iss >> dummy >> atomStr >> energyStr >> xStr >> yStr >> zStr >> vacancyStr;
      if ((vacancyStr=="1") && (find(atomsToTrack.begin(), atomsToTrack.end(), atomStr) != atomsToTrack.end())) {

        tmpRecoilIons.push_back(stoi(atomStr));
        tmpRecoilEnergies.push_back(stof(energyStr));
        tmpRecoil_xLocs.push_back(stof(xStr)*0.1);
        tmpRecoil_yLocs.push_back(stof(yStr)*0.1);
        tmpRecoil_zLocs.push_back(stof(zStr)*0.1);
      }
    }
  }
  outputFile->cd();
  tree->Write("trimTree","TObject::kOverwrite");
  outputFile->Close();
}


void ProcessThrow(CascadeData& data, TTree* tree) {
  float startX = startOffset_nm;
  float startY = 0;
  float startZ = 0;

  float prevX = startX;
  float prevY = startY;
  float prevZ = startZ;

  for (size_t i = 0; i < data.primaryEnergies.size() - 1; i++) {
    //Shouldn't be empty but we'll check
    if (data.recoil_xLocs[i].empty()) continue;

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
    std::array<std::array<float, 3>, 3> rotMat = ComputeRotationMatrix(dx,dy,dz);

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
    std::unordered_map<std::tuple<int, int, int>, int> clusters;
    for (size_t j=0; j < flatX.size(); j++ ) {
      //Rotate all subsequent recoils such that they correspond to the incident particle traveling in the (1,0,0 direction)
      float x = rotMat[0][0]*flatX[j] + rotMat[0][1]*flatY[j] + rotMat[0][2]*flatZ[j];
      float y = rotMat[1][0]*flatX[j] + rotMat[1][1]*flatY[j] + rotMat[1][2]*flatZ[j];
      float z = rotMat[2][0]*flatX[j] + rotMat[2][1]*flatY[j] + rotMat[2][2]*flatZ[j];
    
      int ix = std::round(x / clusteringDistance_nm);
      int iy = std::round(y / clusteringDistance_nm);
      int iz = std::round(z / clusteringDistance_nm);
      auto key = std::make_tuple(ix, iy, iz);
      clusters[key] += 1;
    }

    //5. Clear branch vectors
    xs->clear();
    ys->clear();
    zs->clear();
    nVac->clear();

    //6. Fill
    energy = data.primaryEnergies[i];
    for (const auto& [key, count] : clusters) {
      auto [ix, iy, iz] = key;
  
      // Convert back to physical lattice site
      float x = ix * clusteringDistance_nm;
      float y = iy * clusteringDistance_nm;
      float z = iz * clusteringDistance_nm;
  
      xs->push_back(x);
      ys->push_back(y);
      zs->push_back(z);
      nVac->push_back(count);
    }
    tree->Fill();

    //7. Update prevX, prevY, prevZ for the next iteration
    prevX = data.primary_xLocs[i];
    prevY = data.primary_yLocs[i];
    prevZ = data.primary_zLocs[i];
  }
}

//Uses Rodriguez formulat to calculate the rotation matrix from (dx,dy,dz) to (1,0,0) frame.
//We can then easily apply this to all recoil points
std::array<std::array<float, 3>, 3> ComputeRotationMatrix(float dx, float dy, float dz) {
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