import os
import numpy as np
import pickle
import lzma

#Incident ion
ionName = "4He"
energies = [2150]
nps = 50000

#Mass evaluation from https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
#subtracting off # of protons * 0.000548579905 amu to remove electron contribution
massDict = {
  "19F": 18.993465942925,
  "7Li": 7.014357694545,
  "6Li": 6.013477147705,
  "4He": 4.00150609432,
  "3He": 3.01493216216,
  "3H": 3.015500701415,
  "2H": 2.013553197939,
  "1H": 1.007276451993
}
zDict = {
  "19F": 9,
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

#Target description
target_name = "LiF"
target_nElements = 2
target_nLayers = 1
target_start_offset = 0 #Angstroms, we offset the start point to save backscattered ions

#Note if you want the same element in different layers, you can optionally include it twice to track
#recoils separately and set separate displacement energies
target_element_names = ["Li","F"]
target_element_Zs = [3,9]
target_element_masses = [6.941,18.998]
target_element_TDEs = [25,25] #eV, displacement energy for each element.
target_element_LBEs = [3., 3.] #eV, lattice binding energies
target_element_SBEs = [1.67, 2] #eV, surface binding energies

layer_names = ["LiF"]
layer_depths = [2000000] #Angstroms. (2mm) - we want to be sure ions won't range out
target_depth = sum(layer_depths)
layer_densities  = [2.635] # g/cm^3]
#stoichs correspond to target_elements
layer_stoichs = [[0.5,0.5]]
layer_phases = [1]

#Location of SRIM folder
srimFolder = "C:/Users/Sam/Desktop/SRIM_exe/"
if not srimFolder.endswith("/"):
  srimFolder+="/"

outputFolder = "C:/Users/Sam/Documents/code/trimRunner/outputs/{0}/{1}/".format(target_name,ionName)
if not os.path.exists(outputFolder):
  os.mkdir(outputFolder)

def makeTrimInputString(energy):
  lines=[]
  headerLine = "==> SRIM-2013.00 This file controls TRIM Calculations."
  ionHeaderLine = "Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number."
  ionLine = "     {0}   {1:.3f}         {2}       0  {3} 0    {4}".format(Z,mass,energy,nps,nps+1)
  lines.append(headerLine)
  lines.append(ionHeaderLine)
  lines.append(ionLine)

  cascadeHeaderLine = "Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders"
  cascadeLine = "     2     0     0"
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
  plotLine = "     0     {0}     {1}".format(target_start_offset,target_depth)
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

def writeLines(lines,srimFolder):
  trimFile = open("{0}TRIM.in".format(srimFolder),"w")
  for line in lines:
    if not line.endswith("\n"):
      line+="\n"
      trimFile.write(line)
  trimFile.close()

# ========== MAIN ==========
for energy in energies:
  #Go to SRIM folder, make TRIMIN.txt
  os.chdir(srimFolder)
  lines = makeTrimInputString(energy)
  writeLines(lines,srimFolder)

  #Run SRIM
  os.system("TRIM.exe")  # This runs SRIM in batch mode 

  #Mv output
  os.rename("C:\\Users\\Sam\\Desktop\\SRIM_exe\\SRIM Outputs\\COLLISON.txt", "C:\\Users\\Sam\\Desktop\\SRIM_exe\\SRIM Outputs\\{0}_{1}.txt".format(target_name,ionName))
