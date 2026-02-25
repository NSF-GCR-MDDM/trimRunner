import sys
import gzip
import multiprocessing as mp
import tqdm
import numpy as np
import os
import array
import h5py
import parserUtils

#Usage: python3 <input filename> <output filename> <mode='fast' or 'full'>

def main():
  ################
  ##RUN SETTINGS##
  ################
  #Only save recoils with energies in this range. In general you should start your TRIM run at least a few hundred keV
  #above the max range to allow for 'burn in'
  energyRange_eV = [1,500e3] 
  #[3.7e7,9e6] #for alphas

  #We impose artificial limits in the number of recoils we store within given energy ranges. This keeps us from
  #being overwhelmed by low-energy recoils
  spacing = "log"
  nBins = 1000
  #nBins = 5300 #for alphas
  maxEntriesPerBin = 100

  #RAM/CPU settings
  #maxThrowstoProcessAtOnce = 1000 #for alphas
  nCores = 8
  maxCascadesToHold = 1000000 #, only relevant in 'full' mode. 1e6 cascades ~5GB peak RAM usage

  #################
  ##INPUT PARSING##
  #################
  if len(sys.argv)<4:
    print("Error! Usage is:\n\n\tpython <input file .gz or .txt> <output file name .h5 or .root> <fast or full>\n\nExiting!")
    sys.exit()

  inpFileName = sys.argv[1]
  if inpFileName.endswith(".gz"):
    inpFile = gzip.open(inpFileName, "rt", encoding="utf-8", errors="replace")
  elif inpFileName.endswith(".txt"):
    inpFile = open(inpFileName, "rt", encoding="utf-8", errors="replace")
  else:
    print(f"Error! Acceptable input file types are .gz or .txt.\nYou passed in {inpFileName}!\n")
    sys.exit()

  outputFileName = sys.argv[2]
  if outputFileName.endswith(".root"):
    outputFileType = "root"
  elif outputFileName.endswith(".h5"):
    outputFileType = "h5"
  else:
    print(f"Error! Acceptable output file types are .root or .h5.\nYou passed in {outputFileName}!\n")
  #Get base path of output name, make if doesn't exist
  outDir = os.path.dirname(os.path.abspath(outputFileName))
  if outDir and not os.path.isdir(outDir):
    os.makedirs(outDir, exist_ok=True)

  mode = sys.argv[3]
  if not mode in ["fast","full"]:
    print(f"Error! Third input arg must be 'fast' or 'full'.\nYou pass in {mode}!\n")
    sys.exit
  
  ######################
  ##GET INITIAL ENERGY##
  ######################
  initialEnergy_eV = parserUtils.getInitialEnergy(inpFile)

  ######################################################
  ##DETERMINE MAX ENTRIES TO BE WRITTEN TO OUTPUT FILE##
  ######################################################
  if spacing=="linear":
    bin_edges = np.linspace(energyRange_eV[0],energyRange_eV[1],nBins+1)
  else:
    bin_edges = np.logspace(np.log10(energyRange_eV[0]), np.log10(energyRange_eV[1]), nBins + 1)

  counts_per_bin = np.zeros(len(bin_edges) - 1)

  ##################
  ##PREPARE OUTPUT##
  ##################
  if outputFileType=="root":
    import ROOT
    outputFile = ROOT.TFile(outputFileName,"RECREATE")
    tree = ROOT.TTree("trimTree","")

    ionEnergy_eV = array.array("f",[0.])
    pka_endpoint_x = array.array("f",[0.])
    pka_endpoint_y = array.array("f",[0.])
    pka_endpoint_z = array.array("f",[0.])
    xs = ROOT.std.vector("float")()
    ys = ROOT.std.vector("float")()
    zs = ROOT.std.vector("float")()
    nVacs = ROOT.std.vector("float")()
    displacedZs = ROOT.std.vector("unsigned char")()
    recoilEnergies_eV = ROOT.std.vector("float")()
    recoilNums = ROOT.std.vector("short")()

    tree.Branch("ionEnergy_eV", ionEnergy_eV, "ionEnergy_eV/F")
    tree.Branch("pka_endpoint_x_nm", pka_endpoint_x, "pka_endpoint_x_nm/F")
    tree.Branch("pka_endpoint_y_nm", pka_endpoint_y, "pka_endpoint_y_nm/F")
    tree.Branch("pka_endpoint_z_nm", pka_endpoint_z, "pka_endpoint_z_nm/F")
    tree.Branch("recoilNums", recoilNums)
    tree.Branch("xs_nm", xs)
    tree.Branch("ys_nm", ys)
    tree.Branch("zs_nm", zs)
    tree.Branch("nVacs", nVacs)
    tree.Branch("displacedAtoms_Z", displacedZs)
    tree.Branch("recoilEnergies_eV", recoilEnergies_eV)

    branches = [ionEnergy_eV,pka_endpoint_x,pka_endpoint_y,pka_endpoint_z,xs,ys,zs,nVacs,displacedZs,recoilEnergies_eV,recoilNums]
  else:
    maxEntries = nBins*maxEntriesPerBin
    outputFile = h5py.File(outputFileName, "w")

    dt_float = h5py.vlen_dtype(np.float32)
    dt_uint8 = h5py.vlen_dtype(np.uint8)
    dt_int16 = h5py.vlen_dtype(np.int16)

    d_ionE = outputFile.create_dataset("ionEnergy_eV", shape=(maxEntries,), maxshape=(maxEntries,), dtype=np.float32, chunks=True)
    d_pka_endpoint_x = outputFile.create_dataset("pka_endpoint_x_nm", shape=(maxEntries,), maxshape=(maxEntries,), dtype=np.float32, chunks=True)
    d_pka_endpoint_y = outputFile.create_dataset("pka_endpoint_y_nm", shape=(maxEntries,), maxshape=(maxEntries,), dtype=np.float32, chunks=True)
    d_pka_endpoint_z = outputFile.create_dataset("pka_endpoint_z_nm", shape=(maxEntries,), maxshape=(maxEntries,), dtype=np.float32, chunks=True)

    d_xs = outputFile.create_dataset("xs_nm", shape=(maxEntries,), maxshape=(maxEntries,), dtype=dt_float, chunks=True)
    d_ys = outputFile.create_dataset("ys_nm", shape=(maxEntries,), maxshape=(maxEntries,), dtype=dt_float, chunks=True)
    d_zs = outputFile.create_dataset("zs_nm", shape=(maxEntries,), maxshape=(maxEntries,), dtype=dt_float, chunks=True)

    d_nv = outputFile.create_dataset("nVacs", shape=(maxEntries,), maxshape=(maxEntries,), dtype=dt_float, chunks=True)
    d_disZ = outputFile.create_dataset("displacedAtoms_Z", shape=(maxEntries,), maxshape=(maxEntries,), dtype=dt_uint8, chunks=True)
    d_recoilE = outputFile.create_dataset("recoilEnergies_eV", shape=(maxEntries,), maxshape=(maxEntries,), dtype=dt_float, chunks=True)
    d_recoilNums = outputFile.create_dataset("recoilNums", shape=(maxEntries,), maxshape=(maxEntries,), dtype=dt_int16, chunks=True)

    branches = (d_ionE,d_pka_endpoint_x,d_pka_endpoint_y,d_pka_endpoint_z,d_xs, d_ys, d_zs, d_nv, d_disZ, d_recoilE, d_recoilNums)

  #######################################
  ##PARSE, PROCESS, AND WRITE IN CHUNKS##
  #######################################
  nCascadesInBuffer = 0
  nThrowsProcessed = 0
  nEventsWritten = 0
  throwBuffer = []
  EOF = False
  while not EOF:
    #Read one throw
    if mode=="full":
      primary_steps_arr,cascades_arr,EOF = parserUtils.read_throw_from_csv(inpFile)
    else:
      primary_steps_arr,cascades_arr,EOF = parserUtils.read_throw_from_csv_fast(inpFile)

    nThrowsProcessed+=1
    if nThrowsProcessed%500==0:
      print(f"Processed {nThrowsProcessed} throws")

    #Check if primary_steps_arr is empty, break
    if primary_steps_arr.shape[0]==0:
      break
    
    #Add event to our buffer of events to process
    throwBuffer.append((primary_steps_arr,cascades_arr,initialEnergy_eV))      
    nCascadesInBuffer += cascades_arr.shape[0]

    #If we reached the threshold of events to process, process them
    if nCascadesInBuffer >= maxCascadesToHold or EOF:
      with mp.Pool(processes=nCores) as pool:
        primaries_results = tqdm.tqdm(pool.map(parserUtils.processThrow,throwBuffer,chunksize=10))

      #Write output
      if outputFileType=="root":
        for res in primaries_results:
          parserUtils.fillTree(tree,branches,res,bin_edges,counts_per_bin,maxEntriesPerBin)
      else:
        for res in primaries_results:
          nEventsWritten = parserUtils.fillh5(branches,res,bin_edges,counts_per_bin,maxEntriesPerBin,nEventsWritten) 

      #Clear buffers
      nCascadesInBuffer=0
      throwBuffer.clear()
      
  print(f"Output contains {int(np.sum(counts_per_bin))} out of max of {int(nBins*maxEntriesPerBin)} events\n")

  inpFile.close()

  if outputFileType=="root":
    tree.Write("trimTree",ROOT.TObject.kOverwrite)
    outputFile.Close()
  else:
    #Resize
    for d in branches:
      d.resize((nEventsWritten,))
    outputFile.close()

if __name__ == "__main__":
  main()