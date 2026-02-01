import uproot
import h5py
import numpy as np
import sys
import tqdm
import gc
# Usage: python convert_root_to_hdf5.py input.root output.h5
root_file = sys.argv[1]
output_file = sys.argv[2]
treeName = "trimTree"
print(f"Opening ROOT file: {root_file}")
with uproot.open(root_file) as f:
  tree = f[treeName]
  print("Reading branches...")
  nEntries = tree.num_entries
  ionEnergy_eV = np.empty(nEntries, dtype=np.uint32)
  xs_nm = [None] * nEntries
  ys_nm = [None] * nEntries
  zs_nm = [None] * nEntries
  nVacs = [None] * nEntries
  displacedAtoms_Z = [None] * nEntries
  recoilEnergies_eV = [None] * nEntries

  batch_size = 10000
  for start in tqdm.trange(0, nEntries, batch_size, desc="Reading in batches"):
      stop = min(start + batch_size, nEntries)
      arr = tree.arrays(["ionEnergy_eV", "xs_nm", "ys_nm", "zs_nm", "nVacs","displacedAtoms_Z","recoilEnergies_eV"],
                        entry_start=start, entry_stop=stop,
                        library="np")
      ionEnergy_eV[start:stop] = arr["ionEnergy_eV"]
      xs_nm[start:stop] = arr["xs_nm"]
      ys_nm[start:stop] = arr["ys_nm"]
      zs_nm[start:stop] = arr["zs_nm"]
      nVacs[start:stop] = arr["nVacs"]
      displacedAtoms_Z[start:stop] = arr["displacedAtoms_Z"]
      recoilEnergies_eV[start:stop] = arr["recoilEnergies_eV"]
      
print("Writing HDF5 file...")
with h5py.File(output_file, "w") as h5f:
  h5f.create_dataset("ionEnergy_eV", data=ionEnergy_eV)
  # Define variable-length datatypes
  dt_float = h5py.vlen_dtype(np.float32)
  dt_uint8 = h5py.vlen_dtype(np.uint8)
  dt_int32 = h5py.vlen_dtype(np.int32)
  # Create empty datasets
  dset_xs = h5f.create_dataset("xs_nm", shape=(len(xs_nm),), dtype=dt_float)
  dset_ys = h5f.create_dataset("ys_nm", shape=(len(ys_nm),), dtype=dt_float)
  dset_zs = h5f.create_dataset("zs_nm", shape=(len(zs_nm),), dtype=dt_float)
  dset_nv = h5f.create_dataset("nVacs", shape=(len(nVacs),), dtype=dt_uint8)
  dset_displacedZ = h5f.create_dataset("displacedAtoms_Z", shape=(len(displacedAtoms_Z),), dtype=dt_uint8)
  dset_recoilEs = h5f.create_dataset("recoilEnergies_eV", shape=(len(recoilEnergies_eV),), dtype=dt_int32)
  for i in tqdm.tqdm(range(len(xs_nm)),total=len(xs_nm)):
    dset_xs[i] = xs_nm[i]
    dset_ys[i] = ys_nm[i]
    dset_zs[i] = zs_nm[i]
    dset_nv[i] = nVacs[i]
    dset_displacedZ[i] = displacedAtoms_Z[i]
    dset_recoilEs[i] = recoilEnergies_eV[i]
print(f"Saved HDF5 to {output_file}")
