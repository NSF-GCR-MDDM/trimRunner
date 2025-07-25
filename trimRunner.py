import os
import numpy as np
import pickle
import uuid
import shutil
import gzip
import tarfile

#Incident ion
ionName = "16O"
energy = 5500 #keV
nps = 250000

#Mass evaluation from https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
#subtracting off # of protons * 0.000548579905 amu to remove electron contribution
massDict = {
  "58Fe": 57.91901049747,
  "57Fe": 56.92112887247,
  "56Fe": 55.92067245947,  #91.8%
  "54Fe": 53.92534511147,
  "30Si": 29.96609001833,
  "29Si": 28.96881454567,
  "28Si": 27.96924641575,   #92.2%
  "26Mg": 25.97601001314,
  "25Mg": 24.97925400714,
  "24Mg": 23.97845873014,   #79%
  "19F": 18.993465942925,
  "16O": 15.99052598002,    #99.8%
  "7Li": 7.014357694545,
  "6Li": 6.013477147705,
  "4He": 4.00150609432,
  "3He": 3.01493216216,
  "3H": 3.015500701415,
  "2H": 2.013553197939,
  "1H": 1.007276451993
}
zDict = {
  "58Fe": 26,
  "57Fe": 26,
  "56Fe": 26,
  "30Si": 14,
  "29Si": 14,
  "28Si": 14,
  "26Mg": 12,
  "25Mg": 12,
  "24Mg": 12,
  "19F": 9,
  "16O": 8,
  "7Li": 3,
  "6Li": 3,
  "4He": 2,
  "3He": 2,
  "3H": 1,
  "2H": 1,
  "1H": 1,
}
Z = zDict[ionName]
mass = massDict[ionName]
runMode = 2

#Target description
target_name = "Olivine"
if target_name == "LiF":
  target_nElements = 2
elif target_name == "Olivine":
  target_nElements = 4
target_nLayers = 1
target_start_offset = 0 #Angstroms, we offset the start point to save backscattered ions

#Note if you want the same element in different layers, you can optionally include it twice to track
#recoils separately and set separate displacement energies
if target_name=="LiF":
  target_element_names = ["Li","F"]
  target_element_Zs = [3,9]
  target_element_masses = [6.941,18.998] #Natural abundances, atomic masses
  target_element_TDEs = [25,25] #eV, displacement energy for each element.
  target_element_LBEs = [3., 3.] #eV, lattice binding energies
  target_element_SBEs = [1.67, 2] #eV, surface binding energies
  layer_names = ["LiF"]
  layer_densities  = [2.635] # g/cm^3]
  #stoichs correspond to target_elements
  layer_stoichs = [[0.5,0.5]]
elif target_name == "Olivine":
  target_element_names = ["Fe","Mg","Si","O"]
  target_element_Zs = [26,12,14,8]
  target_element_masses = [55.845,24.305,28.085,15.999]
  #From https://theses.hal.science/tel-01126887v1/file/2014PA112184.pdf
  target_element_TDEs = [25,25,25,25] #eV, displacement energy for each element.
  target_element_LBEs = [1.9,1.9,4,15.5] #eV, lattice binding energies
  target_element_SBEs = [6.47,3.7,6.8,4.75] #eV, surface binding energies
  layer_names = ["Olivine"]
  #https://webmineral.com/data/Olivine.shtml
  layer_densities  = [3.32] # g/cm^3]
  #stoichs correspond to target_elements
  layer_stoichs = [[0.4,1.6,1.0,4.0]]
layer_depths = [20000000] #Angstroms. (2mm) - we want to be sure ions won't range out
target_depth = sum(layer_depths)
layer_phases = [1]

def makeTrimInputString(energy):
  lines=[]
  headerLine = "==> SRIM-2013.00 This file controls TRIM Calculations."
  ionHeaderLine = "Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number."
  ionLine = "     {0}   {1:.3f}         {2}       0  {3} 0    {4}".format(Z,mass,energy,nps,nps+1)
  lines.append(headerLine)
  lines.append(ionHeaderLine)
  lines.append(ionLine)

  cascadeHeaderLine = "Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders"
  cascadeLine = f"     {runMode}     0     0"
  lines.append(cascadeHeaderLine)
  lines.append(cascadeLine)

  #Note, 1 = new file, 2 = extend
  diskHeaderLine = "Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file"
  diskLine = "     0     0     0     0     2     0"
  lines.append(diskHeaderLine)
  lines.append(diskLine)

  targetHeaderLine = "Target material : Number of Elements & Layers"
  targetLine = '"{0} ({1}) into {2}"     {3}     {4}'.format(ionName,energy,target_name,target_nElements,target_nLayers)
  lines.append(targetHeaderLine)
  lines.append(targetLine)
  
  plotHeaderLine = "PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]"
  plotLine = "     5     {0}     {1}".format(target_start_offset,target_depth)
  lines.append(plotHeaderLine)
  lines.append(plotLine)

  targetElementsHeaderLine = "Target Elements:    Z   Mass(amu)"
  lines.append(targetElementsHeaderLine)
  for i,elemName in enumerate(target_element_names):
     targetElementLine = "Atom {0} = {1} =".format(i,elemName)
     nSpaces = 15-len(targetElementLine)
     targetElementLine += " "*nSpaces + "{0}  {1:.3f}".format(target_element_Zs[i],target_element_masses[i])
     lines.append(targetElementLine)

  layerMetaHeaderLine = "Layer   Layer Name /               Width Density"
  for i,elemName in enumerate(target_element_names):
    layerMetaHeaderLine +="     {0}({1})".format(elemName,target_element_Zs[i])
  lines.append(layerMetaHeaderLine)

  layerHeaderLine = "Numb.   Description                (Ang) (g/cm3)"
  for i in range(0,len(layer_names)):
     layerHeaderLine+="    Stoich"
  lines.append(layerHeaderLine)
  for i,layer_name in enumerate(layer_names):
    layerLine = '{0}      "{1}"           {2}  {3:.3f}'.format(i,layer_name,layer_depths[i],layer_densities[i])
    for stoich in layer_stoichs[i]:
      layerLine+="    {0}".format(stoich)
    lines.append(layerLine)

  layerPhaseHeaderLine = "Target layer phases (0=Solid, 1=Gas)"
  layerPhaseLine=""
  for layer_phase in layer_phases:
    layerPhaseLine+="{0} ".format(layer_phase)
  lines.append(layerPhaseHeaderLine)
  lines.append(layerPhaseLine)
  
  layerBraggCorrectionHeaderLine = "Target Compound Corrections (Bragg)"
  layerBraggCorrectionLine=""
  for layer_phase in layer_phases:
    layerBraggCorrectionLine+="1 "
  lines.append(layerBraggCorrectionHeaderLine)
  lines.append(layerBraggCorrectionLine)
  
  targetElementTDEHeaderLine = "Individual target atom displacement energies (eV)"
  lines.append(targetElementTDEHeaderLine)
  targetElementTDELine = ""
  for TDE in target_element_TDEs:
    targetElementTDELine+="{0} ".format(TDE)  
  lines.append(targetElementTDELine)

  targetElementLBEHeaderLine = "Individual target atom lattice binding energies (eV)"
  lines.append(targetElementLBEHeaderLine)
  targetElementLBELine = ""
  for LBE in target_element_LBEs:
    targetElementLBELine+="{0} ".format(LBE)
  lines.append(targetElementLBELine)

  targetElementSBEHeaderLine = "Individual target atom surface binding energies (eV)"
  lines.append(targetElementSBEHeaderLine)
  targetElementSBELine = ""
  for SBE in target_element_SBEs:
    targetElementSBELine+="{0} ".format(SBE)
  lines.append(targetElementSBELine)

  stoppingPowerHeaderline="Stopping Power Version (1=2008, 0=2008)"
  stoppingPowerLine=" 1"
  lines.append(stoppingPowerHeaderline)
  lines.append(stoppingPowerLine)
  return lines

if __name__ == "__main__":
  # Create a unique working directory for this ion run, cp srim to that folder
  unique_id = uuid.uuid4().hex[:8]
  temp_workdir = os.path.join("C:/Users/Sam/Desktop/SRIM_tmp", f"{ionName}_{unique_id}")
  shutil.copytree("C:/Users/Sam/Desktop/SRIM_exe", temp_workdir)
  srimFolder = temp_workdir + "/"

  # Make output folder, output file names
  outputFolder = "C:/Users/Sam/Documents/code/trimRunner/data/txt/"
  if not os.path.exists(outputFolder):
    os.mkdir(outputFolder)

  txt_name = "{0}_{1}.txt".format(target_name, ionName)
  txt_path = os.path.join(outputFolder, txt_name)

  # Overwrite if exists
  if os.path.exists(txt_path):
    os.remove(txt_path)
  tar_path =  "{0}_{1}.tar.gz".format(target_name, ionName)
  if os.path.exists(tar_path):
    os.remove(tar_path)

  # Create new SRIM input
  os.chdir(srimFolder)
  lines = makeTrimInputString(energy)
  with open(f"{srimFolder}TRIM.in", "w") as trimFile:
    for line in lines:
      if not line.endswith("\n"):
        line += "\n"
      trimFile.write(line)
    trimFile.close()

  # Run SRIM  
  exit_code = os.system("TRIM.exe")

  # Upon successful completion, mv collison.txt to file name, compress, remove temp folder
  if exit_code == 0:
    print("TRIM completed successfully.")

    # Move and compress output
    src = os.path.join(srimFolder, "SRIM Outputs", "COLLISON.txt")
    os.replace(src, txt_path)

    # Create .tar.gz
    with tarfile.open(tar_path, "w:gz") as tar:
      tar.add(txt_path, arcname=os.path.basename(txt_path))
      tar.close()
