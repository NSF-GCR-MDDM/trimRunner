#Mass evaluation from https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt
#subtracting off # of protons * 0.000548579905 amu to remove electron contribution
import pandas as pd
import json
import sys
import uuid
import os
import shutil

def createMassDict(massFile="mass_1.mas20.txt"):
  m_e_amu = 0.000548579905 
  widths = [1,3,5,5,5,1,3,4,1,14,12,13,1, #a1,i3,i5,i5,i5,1x,a3,a4,1x,f14.6,f12.6,f13.5,1x,
            10,1,2,13,11,1,3,1,13,12]     #f10.5,1x,a2,f13.5,f11.5,1x,i3,1x,f13.6,f12.6
  names = ["skip1","N-Z","N","Z","A","skip2","El","O","skip3","mass_excess_keV","mass_excess_unc","binding_energy_keV","skip4",
           "binding_energy_unc","skip5","skip6","beta_decay_energy_keV","beta_decay_energy_unc","skip7","mass_1_uamu","skip8","mass_2_uamu","mass_unc"]
  df = pd.read_fwf(massFile, widths=widths, names=names, skiprows=37,dtype=str)

  df = df[["A","Z","El","mass_1_uamu","mass_2_uamu"]].copy()
  df["A"] = df["A"].str.strip()
  df["El"] = df["El"].str.strip()
  df["symbol"] = df["A"]+df["El"]
  df["Z"] = df["Z"].astype(int)
  df["mass_1_uamu"] = df["mass_1_uamu"].str.replace("#", "")
  df["mass_2_uamu"] = df["mass_2_uamu"].str.replace("#", "")
  df["mass_amu"] = (df["mass_1_uamu"].astype(float)*1e6 + df["mass_2_uamu"].astype(float))*1e-6 - df["Z"]*m_e_amu
  df = df[["symbol","Z","mass_amu"]].dropna()

  mass_dict = {}
  for _, row in df.iterrows():
    mass_dict[row["symbol"]] = {"Z": row["Z"], "mass_amu": row["mass_amu"]}

  return mass_dict

def parseConfig(inpFileName):
  with open(inpFileName, "r") as f:
    lines = f.readlines()
  lines = [line.split("#", 1)[0] for line in lines]  # drop inline comments
  text = "".join(lines)
  data = json.loads(text)
  return data

def checkTrimArgs(data,mass_dict,materials_dict):
  required_keys = ["material","ionSymbol","runMode","calcMode","nps","energy_keV","outputPath"]
  #Check all top level keys are there
  for key in required_keys:
    if key not in data:
      print(f"Error! Config missing required key: '{key}'")
      sys.exit()

  #Check material is valid
  if not data["material"] in materials_dict:
    print(f"Error! Material {data['material']} not supplied in materials_dict!")
    print("Choices are:")
    for key in materials_dict.keys():
      print(f"\t{key}")
    sys.exit()
  
  #Check ion is valid
  if not data["ionSymbol"] in mass_dict:
    print(f"Error! Ion {data['ionSymbol']} not valid!")
    sys.exit()

  #Check calc mode is valid
  if not data["calcMode"] in ["quick","full"]:
    print(f"Error! calcMode {data['calcMode']} not valid!")
    print(f"Valid choices are: {['quick','full']}")
    sys.exit()

  #Check run mode is valid
  if not data["runMode"] in ["damage","efficiency"]:
    print(f"Error! runMode {data['runMode']} not valid!")
    print(f"Valid choices are: {['damage','efficiency']}")
    sys.exit()

  #Create output path if needed
  outdir = data["outputPath"]
  os.makedirs(outdir, exist_ok=True)

def makeTrimInputString(data,mass_dict,material_dict):
  calcMode_str = data["calcMode"]
  if calcMode_str == "full":
    calcMode = 2
  elif calcMode_str == "quick":
    calcMode = 1

  nps = data["nps"]
  energy = data["energy_keV"]
  ion_name = data["ionSymbol"]
  target_name = data["material"]

  Z = mass_dict[ion_name]["Z"]
  mass = mass_dict[ion_name]["mass_amu"]
  target_nElements = material_dict[target_name]["nElements"]
  target_element_names = material_dict[target_name]["element_names"]
  target_element_Zs = material_dict[target_name]["Zs"]
  target_element_masses = material_dict[target_name]["element_masses"]
  layer_names = [target_name]
  layer_stoichs = material_dict[target_name]["stoichs"]
  layer_densities = material_dict[target_name]["densities"]
  target_element_TDEs = material_dict[target_name]["TDEs"]
  target_element_LBEs = material_dict[target_name]["LBEs"]
  target_element_SBEs = material_dict[target_name]["SBEs"]

  target_nLayers = 1
  target_start_offset = 0

  layer_depths = [20000000] #Angstroms. (2mm) - we want to be sure ions won't range out
  target_depth = sum(layer_depths)
  layer_phases = [0]

  lines=[]
  headerLine = "==> SRIM-2013.00 This file controls TRIM Calculations."
  ionHeaderLine = "Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number."
  ionLine = "     {0}   {1:.3f}         {2}       0  {3} 0    {4}".format(Z,mass,energy,nps,nps+1)
  lines.append(headerLine)
  lines.append(ionHeaderLine)
  lines.append(ionLine)

  cascadeHeaderLine = "Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders"
  cascadeLine = f"     {calcMode}     0     0"
  lines.append(cascadeHeaderLine)
  lines.append(cascadeLine)

  #Note, 1 = new file, 2 = extend
  diskHeaderLine = "Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file"
  diskLine = "     0     0     0     0     2     0"
  lines.append(diskHeaderLine)
  lines.append(diskLine)

  targetHeaderLine = "Target material : Number of Elements & Layers"
  targetLine = '"{0} ({1}) into {2}"     {3}     {4}'.format(ion_name,energy,target_name,target_nElements,target_nLayers)
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

def makeTempSRIMFolder(tag,tmp_path,srim_path):
  
  srim_path = os.path.abspath(os.path.normpath(srim_path))
  tmp_path = os.path.abspath(os.path.normpath(tmp_path))

  print("SRIM_EXE_PATH:", srim_path)
  print("SRIM_EXE_PATH isdir:", os.path.isdir(srim_path))
  print("SRIM_EXE_PATH contents count:", len(os.listdir(srim_path)))
  
  unique_id = uuid.uuid4().hex[:8]
  temp_workdir = os.path.join(tmp_path, f"{tag}_{unique_id}")
  temp_workdir = os.path.abspath(os.path.normpath(temp_workdir))
  print("TEMP_WORKDIR:", temp_workdir)

  # Ensure destination does NOT exist. If it does, delete it.
  if os.path.exists(temp_workdir):
    shutil.rmtree(temp_workdir)

  # Do not swallow errors. Let it crash with a traceback if it fails.
  shutil.copytree(srim_path, temp_workdir, dirs_exist_ok=False)

  print("TEMP_WORKDIR contents count:", len(os.listdir(temp_workdir)))
  print("TEMP_WORKDIR first entries:", os.listdir(temp_workdir)[:20])

  import subprocess
  subprocess.run(["explorer", os.path.abspath(temp_workdir)])

  return temp_workdir

def runSRIM(srimFolder,input_data,mass_dict,materials_dict):
  os.chdir(srimFolder)
  print(srimFolder)
  lines = makeTrimInputString(input_data, mass_dict, materials_dict)
  with open(os.path.join(srimFolder, "TRIM.in"), "w", encoding="utf-8", newline="\n") as f:
    for line in lines:
      if not line.endswith("\n"):
        line += "\n"
      f.write(line)

  sys.exit()
  exit_code = os.system("TRIM.exe")
  if exit_code != 0:
    raise RuntimeError(f"TRIM.exe failed with exit_code={exit_code}")

  collison_path = os.path.join(srimFolder, "SRIM Outputs", "COLLISON.txt")
  if not os.path.exists(collison_path):
    raise RuntimeError("TRIM completed but COLLISON.txt not found")

  return collison_path

def createMaterialsDict():
  materials_dict = {}
  #LiF
  materials_dict["LiF"] = {
    "nElements": 2,
    "element_names": ["Li","F"],
    "element_masses": [6.941,18.998], #Natural abundances, atomic masses
    "Zs": [3,9],
    #From TRIM Compound Library
    "TDEs": [25.,25.],
    "LBEs": [3.,3.],
    "SBEs": [1.67,2.],
    "densities": [2.635],
    "stoichs": [[0.5,0.5]]}
  #Olivine - (Fe0.4Mg1.6)SiO4
  materials_dict["Olivine"] = {
    "nElements": 4,
    "element_names": ["Fe","Mg","Si","O"],
    "element_masses": [55.845,24.305,28.085,15.999], #Natural abundances, atomic masses
    "Zs": [26,12,14,8],
    #From https://theses.hal.science/tel-01126887v1/file/2014PA112184.pdf
    "TDEs": [25,25,25,25], #eV, displacement energy for each element.
    "LBEs": [1.9,1.9,4,15.5], #eV, lattice binding energies
    "SBEs": [6.47,3.7,6.8,4.75], #eV, surface binding energies
    #https://webmineral.com/data/Olivine.shtml
    "densities": [3.32],
    #stoichs correspond to target_elements
    "stoichs": [[0.4,1.6,1.0,4.0]]}
  #Diamond
  materials_dict["Diamond"] = {
    "nElements": 1,
    "element_names": ["C"],
    "element_masses": [12.01], #Natural abundances, atomic masses
    "Zs": [6],
    "TDEs": [31.], #eV, displacement energy for each element. From https://arxiv.org/pdf/2206.06772
    "LBEs": [3.], #eV, lattice binding energies. From https://www.sciencedirect.com/science/article/abs/pii/S0168583X21001890
    "SBEs": [7.4], #eV, surface binding energies. From https://www.sciencedirect.com/science/article/abs/pii/S0168583X21001890
    #https://www.sciencedirect.com/science/article/abs/pii/S0168583X21001890 
    "densities": [3.51],
    #stoichs correspond to target_elements
    "stoichs": [[1.0]]}
  return materials_dict