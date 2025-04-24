import os
import shutil
import subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

# USER CONFIGURATION
ionName = "7Li"
Z = 3
mass = 7.01600343  # amu
nps = 500
energies = list(range(1050, 0, -100))  # in keV
target_name = "LiF"
target_nElements = 2
target_nLayers = 1
target_start_offset = 10  # Angstroms
target_element_names = ["Li", "F"]
target_element_Zs = [3, 9]
target_element_masses = [6.941, 18.998]
target_element_TDEs = [25, 25]  # eV
target_element_LBEs = [5., 6.]  # eV
target_element_SBEs = [4.0, 4.5]  # eV
layer_names = ["LiF"]
layer_depths = [20000000]
target_depth = sum(layer_depths)
layer_densities = [2.635]
layer_stoichs = [[0.5, 0.5]]
layer_phases = [1]

srim_base = "C:/Users/Sam/Desktop/SRIM_exe/"
out_base = "C:/Users/Sam/Documents/code/trimRunner/outputs/LiF/{0}/".format(ionName)
temp_base = "C:/Users/Sam/Documents/code/trimRunner/tmp_runs/"

# Ensure output dirs exist
os.makedirs(out_base, exist_ok=True)
os.makedirs(temp_base, exist_ok=True)

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

def writeLines(lines, srimFolder):
    trimFile = open(os.path.join(srimFolder, "TRIM.in"), "w")
    for line in lines:
        if not line.endswith("\n"):
            line += "\n"
        trimFile.write(line)
    trimFile.close()
    
def run_single_energy(energy):
  global Z, mass, nps, target_name, target_start_offset
  global target_element_names, target_element_Zs, target_element_masses
  global target_element_TDEs, target_element_LBEs, target_element_SBEs
  global layer_names, layer_depths, target_depth, layer_densities
  global layer_stoichs, layer_phases
  # Set up unique temp folder for this energy
  workdir = os.path.join(temp_base, f"{ionName}_{energy}keV")
  if os.path.exists(workdir):
      shutil.rmtree(workdir)
  shutil.copytree(srim_base, workdir)

  # Generate TRIM.IN and run TRIM
  lines = makeTrimInputString(energy)
  writeLines(lines, workdir)
  subprocess.run(["TRIM.exe"], cwd=workdir, shell=True)

  # Parse output
  collision_file = os.path.join(workdir, "SRIM Outputs", "COLLISON.txt")
  events = parseCollisionC(collision_file)

  # Save output compressed
  out_filename = os.path.join(out_base, f"{energy/1000.:06.3f}_keV.npz")
  np.savez_compressed(out_filename, **events)

  return f"{energy} keV done"

if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=multiprocessing.cpu_count()-1) as pool:
        results = list(pool.map(run_single_energy, energies))
    for r in results:
        print(r)
