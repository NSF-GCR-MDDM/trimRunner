import os
import numpy as np
import pickle
import lzma

ionName = "7Li"
Z = 3
mass = 7.01600343 #amu
energies = range(25,2050,50) #keV
nps = 500

target_name = "LiF"
target_nElements = 2
target_nLayers = 1
target_start_offset = 10 #Angstroms, we offset the start point to save backscattered ions

#Note if you want the same element in different layers, you can optionally include it twice to track
#recoils separately and set separate displacement energies
target_element_names = ["Li","F"]
target_element_Zs = [3,9]
target_element_masses = [6.941,18.998]
target_element_TDEs = [25,25] #eV, displacement energy for each element.
target_element_LBEs = [5., 6.] #eV, lattice binding energies
target_element_SBEs = [4.0, 4.5] #eV, surface binding energies

layer_names = ["LiF"]
layer_depths = [20000000] #Angstroms. (2mm) - we want to be sure ions won't range out
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
     targetElementLine += " "*nSpaces + "{0}  {1:.3f}".format(target_element_masses[i],target_element_Zs[i])
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

#Set up to parse collision files of type C (Full recoil cascades)
def parseCollisionC(filename):
  startLine=10
  inpFile = open(filename,"r")

  lookingForRecoils=0

  ionNum = 0
  ionNums=[]
  recoilAtoms=[]
  recoilEnergies=[]
  recoilXs=[]
  recoilYs=[]
  recoilZs=[]
  recoilVacancies=[]
  recoilReplacements=[]

  for iline,line in enumerate(inpFile):
    
    #Skip header
    if iline<28:
      continue
    elif "For Ion" in line:      
      ionNum+=1
      lookingForRecoils=0

    elif "Recoil Atom" in line:
      lookingForRecoils=1 
    elif "===" in line:
      lookingForRecoils=0
    elif lookingForRecoils==1:
      line=line.strip("\n")
      line = line.replace('³', ' ')         # Replace weird delimiters with space
      line = line.replace('Û', '')          # Remove the garbage Û characters
      lineParts = line.strip().split()      # Split on whitespace into clean fields
      if lookingForRecoils==1:
        if int(lineParts[6])>0:
          ionNums.append(ionNum)
          recoilAtoms.append(int(lineParts[1]))
          recoilEnergies.append(float(lineParts[2])/1.e-6) #MeV
          recoilXs.append(float(lineParts[3])*0.1 - target_start_offset*0.1)
          recoilYs.append(float(lineParts[4])*0.1 - 0)
          recoilZs.append(float(lineParts[5])*0.1 - 0)
          recoilVacancies.append(int(lineParts[6]))
          recoilReplacements.append(int(lineParts[7]))
    else:
      continue
  
  events = {
      "historyNum": np.array(ionNums, dtype=np.uint16),
      "secondaryAtom": np.array(recoilAtoms, dtype=np.uint8),
      "energy": np.array(recoilEnergies, dtype=np.float32),
      "x": np.array(recoilXs, dtype=np.float32),
      "y": np.array(recoilYs, dtype=np.float32),
      "z": np.array(recoilZs, dtype=np.float32),
      "nVacancies": np.array(recoilVacancies, dtype=np.uint16),
  }
  return events

# ========== MAIN ==========
for energy in energies:
  #Go to SRIM folder, make TRIMIN.txt
  os.chdir(srimFolder)
  lines = makeTrimInputString(energy)
  writeLines(lines,srimFolder)

  #Run SRIM
  os.system("TRIM.exe")  # This runs SRIM in batch mode 

  #Make output file
  events = parseCollisionC(srimFolder+"SRIM Outputs/COLLISON.txt")

  outFilename = outputFolder + "{01:06.3f}_keV.npz".format(energy/1000.)
  np.savez_compressed(outFilename, **events)
