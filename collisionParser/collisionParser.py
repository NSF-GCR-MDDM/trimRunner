import sys
import gzip
import multiprocessing as mp
import tqdm
import numpy as np
import ROOT
import os
import array
import h5py
import parserUtils

def _processThrow_star(args):
  return parserUtils.processThrow(*args)

def _processFastThrow_start(args):
  return parserUtils.processFastThrow(*args)

def main():
  ################
  ##RUN SETTINGS##
  ################
  #Only save recoils with energies in this range. In general you should start your TRIM run at least a few hundred keV
  #above the max range to allow for 'burn in'
  energyRange_eV = [20,500e3] 

  #We impose artificial limits in the number of recoils we store within given energy ranges. This keeps us from
  #being overwhelmed by low-energy recoils
  spacing = "log"
  nBins = 1000
  maxEntriesPerBin = 100

  #RAM/CPU settings
  maxThrowsToProcessAtOnce=2000 #Reduce to decrease ram usage, only relevant for 'full' mode. Recommend starting at 100 and tuning
  nCores = 8

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

  ##################
  ##PREPARE OUTPUT##
  ##################
  if outputFileType=="root":
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
    outputFile = h5py.File(outputFileName, "w")

    # vlen types (now floats for nVacs and recoilEnergies per your change)
    dt_float = h5py.vlen_dtype(np.float32)
    dt_uint8 = h5py.vlen_dtype(np.uint8)
    dt_int16 = h5py.vlen_dtype(np.int16)

    d_ionE = outputFile.create_dataset("ionEnergy_eV", shape=(0,), maxshape=(None,), dtype=np.float32)
    d_pka_endpoint_x = outputFile.create_dataset("pka_endpoint_x_nm", shape=(0,), maxshape=(None,), dtype=np.float32)
    d_pka_endpoint_y = outputFile.create_dataset("pka_endpoint_y_nm", shape=(0,), maxshape=(None,), dtype=np.float32)
    d_pka_endpoint_z = outputFile.create_dataset("pka_endpoint_z_nm", shape=(0,), maxshape=(None,), dtype=np.float32)
    d_xs = outputFile.create_dataset("xs_nm", shape=(0,), maxshape=(None,), dtype=dt_float)
    d_ys = outputFile.create_dataset("ys_nm", shape=(0,), maxshape=(None,), dtype=dt_float)
    d_zs = outputFile.create_dataset("zs_nm", shape=(0,), maxshape=(None,), dtype=dt_float)

    d_nv = outputFile.create_dataset("nVacs", shape=(0,), maxshape=(None,), dtype=dt_float)            # float32
    d_disZ = outputFile.create_dataset("displacedAtoms_Z", shape=(0,), maxshape=(None,), dtype=dt_uint8)
    d_recoilE = outputFile.create_dataset("recoilEnergies_eV", shape=(0,), maxshape=(None,), dtype=dt_float)  # float32
    d_recoilNums = outputFile.create_dataset("recoilNums", shape=(0,), maxshape=(None,), dtype=dt_int16)

    branches = (d_ionE,d_pka_endpoint_x,d_pka_endpoint_y,d_pka_endpoint_z,d_xs, d_ys, d_zs, d_nv, d_disZ, d_recoilE, d_recoilNums)


  #######################################
  ##PARSE, PROCESS, AND WRITE IN CHUNKS##
  #######################################
  if spacing=="linear":
    bin_edges = np.linspace(energyRange_eV[0],energyRange_eV[1],nBins+1)
  else:
    bin_edges = np.logspace(np.log10(energyRange_eV[0]), np.log10(energyRange_eV[1]), nBins + 1)

  counts_per_bin = np.zeros(len(bin_edges) - 1)
  EOF = False
  nThrows=0

  if mode=="full":
    while not EOF:
      #Parse
      primaries_df,cascades_df,EOF = parserUtils.read_csv(inpFile,maxThrowsToProcessAtOnce)

      #Check if at end, give user feedback
      if len(primaries_df)==0:
        break
      nThrows += primaries_df["ionNum"].nunique()
      print(f"Parsed {nThrows} throws")

      #Processs
      # Split dfs into groups by ionNum
      primaries_list = {ionNum: df for ionNum, df in primaries_df.groupby("ionNum", sort=False)}
      cascades_list  = {ionNum: df for ionNum, df in cascades_df.groupby("ionNum", sort=False)}

      items = []
      for ionNum, p_df in primaries_list.items():
        c_df = cascades_list.get(ionNum)
        if c_df is None:
          continue
        items.append((p_df, c_df, initialEnergy_eV, bin_edges, counts_per_bin, maxEntriesPerBin))

      #Multiprocessing call
      with mp.Pool(processes=nCores) as pool:
        results = list(tqdm.tqdm(
        pool.imap(_processThrow_star, items, chunksize=100),
        total=len(items),
        desc="Processing throws"
      ))
        
      #Write
      if outputFileType=="root":
        parserUtils.fillTree(tree,branches,results,bin_edges,counts_per_bin,maxEntriesPerBin)
      else:
        events = []
        for df in results:
          if df is not None and len(df) > 0:
            events.extend(df.to_dict("records"))
        parserUtils.fillh5(branches, events, bin_edges, counts_per_bin, maxEntriesPerBin) 
  else:
    primaries_df = parserUtils.read_csv_fast_mode(inpFile)
    primaries_list = [df for _, df in primaries_df.groupby("ionNum", sort=False)]
    items = [(p, initialEnergy_eV, energyRange_eV[0],energyRange_eV[1]) for p in primaries_list]
    #Multiprocessing call
    with mp.Pool(processes=nCores) as pool:
      results = list(tqdm.tqdm(
      pool.imap(_processFastThrow_start, items, chunksize=100),
      total=len(items),
      desc="Processing throws"
    ))
    #Write
    if outputFileType=="root":
      parserUtils.fillTree(tree,branches,results,bin_edges,counts_per_bin,maxEntriesPerBin)
    else:
      events = []
      for df in results:
        if df is not None and len(df) > 0:
          events.extend(df.to_dict("records"))
      parserUtils.fillh5(branches, events, bin_edges, counts_per_bin, maxEntriesPerBin) 

  print(f"Output contains {int(np.sum(counts_per_bin))} out of max of {int(nBins*maxEntriesPerBin)} events\n")

  inpFile.close()

  if outputFileType=="root":
    tree.Write("trimTree",ROOT.TObject.kOverwrite)
    outputFile.Close()
  else:
    outputFile.close()

if __name__ == "__main__":
  main()