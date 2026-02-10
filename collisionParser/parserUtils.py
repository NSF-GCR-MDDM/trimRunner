import pandas as pd
import numba
import numpy as np
import io

def getInitialEnergy(inpFile):
  for line in inpFile:
    if "Ion Energy =" in line:
      start = line.find("=")+1
      end = line.find("keV")
      throwEnergy_eV = float(line[start:end].strip())*1e3
      return throwEnergy_eV

def read_csv(inpFile,maxThrowsToProcessAtOnce):
  nThrowsRead=0
  cascadeNum=0
  lookingForRecoils=False
  primeRecoils = []
  cascades = []

  while True:
    line = inpFile.readline()
    #Check for EOF
    if line == "":
      break
    #Skip comments
    if line.startswith("==") or line.startswith("--"):
      continue
    #Replace weird chars
    line = line.replace("�"," ")

    #One step of the primary recoil
    if "Start of New Cascade" in line:
      lineParts = line.split()      
      primeRecoils.append({
        "ionNum":int(lineParts[0]),
        "energy_eV":float(lineParts[1])*1000,
        "x_nm":float(lineParts[2])*0.1,
        "y_nm":float(lineParts[3])*0.1,
        "z_nm":float(lineParts[4])*0.1,
        "recoilEnergy_eV":float(lineParts[7])
      })
    #Start of a cascade
    elif "Prime Recoil" in line:
      lookingForRecoils=True
      cascadeNum+=1
    #End of a cascade
    elif "Summary of Above Cascade" in line:
      lookingForRecoils=False
    #End of a throw
    elif "For Ion" in line:
      cascadeNum=0
      nThrowsRead+=1
      if nThrowsRead==maxThrowsToProcessAtOnce:
        #Convert to data frame
        primeRecoils_df = pd.DataFrame(primeRecoils)
        cascades_df = pd.DataFrame(cascades)
        return primeRecoils_df,cascades_df,False

    #Parse cascades if this is in a cascades line
    if lookingForRecoils==True:
      lineParts = line.split()
      nVacs = int(lineParts[6])
      nReps = int(lineParts[7])
      if nVacs > 0 and nReps==0:
        cascades.append({
          "ionNum":primeRecoils[-1]["ionNum"],
          "cascadeNum":cascadeNum,
          "atom": int(lineParts[1]),
          "recoilEnergy_eV": float(lineParts[2]),
          "x_nm": float(lineParts[3])*0.1,
          "y_nm": float(lineParts[4])*0.1,
          "z_nm": float(lineParts[5])*0.1,
          "nVacs": int(lineParts[6])
        })

  #Parse final chunk
  primeRecoils_df = pd.DataFrame(primeRecoils)
  cascades_df = pd.DataFrame(cascades)
  return primeRecoils_df,cascades_df,True

def read_csv_fast_mode(inpFile):
  nIons=0
  startedEvents = False
  output = "ionNum energy_keV x_A y_A z_A Se atom_hit recoilEnergy_eV nVacs nReps nInters\n"
  ignore_keys = ["Displacement","Replacements","Vacancies","Interstitials","Sputtered","Transmitted","Ion","Num"]
  while True:
    line = inpFile.readline()
    #Check for EOF
    if line == "":
      break
    #Skip comments
    if line.startswith("==") or line.startswith("--"):
      continue
    #Skip header
    if "REPLAC INTER" in line:
      startedEvents=True
      nIons+=1
      if nIons%1000==0:
        print(nIons)
      continue

    if startedEvents==True:
      if any(k in line for k in ignore_keys):
        continue

      #Replace weird chars
      line = line.replace("�"," ")

      output+=line
      
  df = pd.read_csv(io.StringIO(output), sep=r"\s+",usecols=["ionNum", "energy_keV", "x_A", "y_A", "z_A", "atom_hit", "recoilEnergy_eV", "nVacs"])
  return df

def processThrow(primaries_df,cascades_df,initialEnergy_eV,bin_edges,counts_per_bin,maxEntriesPerBin):
  outputEvents = []

  #Get primaries columns
  primaries_energy_eV = primaries_df["energy_eV"].to_numpy(copy=False)
  primaries_recoilEnergy_eV = primaries_df["recoilEnergy_eV"].to_numpy(copy=False)
  primaries_xs_nm = primaries_df["x_nm"].to_numpy(copy=False)
  primaries_ys_nm = primaries_df["y_nm"].to_numpy(copy=False)
  primaries_zs_nm = primaries_df["z_nm"].to_numpy(copy=False)

  cascade_nums = cascades_df["cascadeNum"].to_numpy(copy=False)
  cascades_atoms = cascades_df["atom"].to_numpy(copy=False)
  cascade_xs_nm = cascades_df["x_nm"].to_numpy(copy=False)
  cascade_ys_nm = cascades_df["y_nm"].to_numpy(copy=False)
  cascade_zs_nm = cascades_df["z_nm"].to_numpy(copy=False)
  cascade_nVacs = cascades_df["nVacs"].to_numpy(copy=False)
  cascade_recoilEnergy_eV = cascades_df["recoilEnergy_eV"].to_numpy(copy=False)

  #Get the start positions of what will become each of our new primaries. Assume the first starts at (0,0,0), and each 
  #subsequent primary starts at the collision site of the previous primary
  primaries_start_xs_nm = np.empty_like(primaries_xs_nm)
  primaries_start_ys_nm = np.empty_like(primaries_ys_nm)
  primaries_start_zs_nm = np.empty_like(primaries_zs_nm)
  primaries_start_xs_nm[0] = 0.0
  primaries_start_ys_nm[0] = 0.0
  primaries_start_zs_nm[0] = 0.0
  if primaries_start_xs_nm.size > 1:
    primaries_start_xs_nm[1:] = primaries_xs_nm[:-1]
    primaries_start_ys_nm[1:] = primaries_ys_nm[:-1]
    primaries_start_zs_nm[1:] = primaries_zs_nm[:-1]
  
  #Calculate the direction it traveled, assuming a straight line between the start point and its collision point. Normalize
  vxs = primaries_xs_nm - primaries_start_xs_nm
  vys = primaries_ys_nm - primaries_start_ys_nm
  vzs = primaries_zs_nm - primaries_start_zs_nm
  #If the primary didn't move between collisions, assign same direction as previous track. This happens near the ends of tracks
  #very seldomly
  if vxs.size > 0:
    for i in range(0, vxs.size):
      if (vxs[i] == 0.0) and (vys[i] == 0.0) and (vzs[i] == 0.0):
        if i>0:
          vxs[i] = vxs[i-1]
          vys[i] = vys[i-1]
          vzs[i] = vzs[i-1]
        else: #Shouldn't ever happen, but if the first primary interaction is at (0,0,0) assign the default SRIM direction of (+1,0,0)
          vxs[i] = 1.
          vys[i] = 0.
          vzs[i] = 0.
  norms = np.sqrt(vxs*vxs + vys*vys + vzs*vzs)
  primaries_vxs = vxs/norms
  primaries_vys = vys/norms
  primaries_vzs = vzs/norms

  #Calculate the energy each primary started with. In SRIM at each collision site we get the energy of the primary just prior
  #to the collision. To get the energy at the initial location it was generated, we need to take the energy of the previous
  #primary and subtract the recoil energy. The first primary starts out with initialEnergy_eV
  primaries_startEnergy_eV = np.empty_like(primaries_energy_eV)
  primaries_startEnergy_eV[0] = float(initialEnergy_eV)
  if primaries_startEnergy_eV.size > 1:
    primaries_startEnergy_eV[1:] = primaries_energy_eV[:-1] - primaries_recoilEnergy_eV[:-1]

  #Now step through primaries. Check if we should process it (ie is the bin already full we skip, or if the primary energy is
  #outside of our bin edges, skip). 
  for i in range(0,len(primaries_df)):
    #Check if this primary is out of bounds or the bin is already filled
    primaryEnergy_eV = float(primaries_startEnergy_eV[i])
    bin_idx = np.searchsorted(bin_edges, primaryEnergy_eV, side="right") - 1
    if bin_idx < 0 or bin_idx >= (len(bin_edges) - 1):
      continue
    if counts_per_bin[bin_idx] >= maxEntriesPerBin:
      continue
    
    #Get all downstream cascades from this recoil. If there are none, skip it
    mask = (cascade_nums >= (i + 1))
    if not np.any(mask):
      continue

    #Compute rotation matrix to align the new lower energy primary in the +X direction
    R = compute_rotation_matrix(primaries_vxs[i],primaries_vys[i],primaries_vzs[i])

    #Get flat arrays of cascade positions, offset them so they start at 0,0,0
    xs = cascade_xs_nm[mask] - primaries_start_xs_nm[i]
    ys = cascade_ys_nm[mask] - primaries_start_ys_nm[i]
    zs = cascade_zs_nm[mask] - primaries_start_zs_nm[i]
    #Rotate all the cascades so they correspond to the +X axis
    rotate_points_inplace(xs, ys, zs, R)

    #Get the PKA endpoint - we reset each time because we are rotating in-place
    pka_end_x = np.array([primaries_xs_nm[-1] - primaries_start_xs_nm[i]])
    pka_end_y = np.array([primaries_ys_nm[-1] - primaries_start_ys_nm[i]])
    pka_end_z = np.array([primaries_zs_nm[-1] - primaries_start_zs_nm[i]])
    #Rotate the pka endpoint
    rotate_points_inplace(pka_end_x, pka_end_y, pka_end_z, R)

    #Get energies, nums
    recoil_energies = cascade_recoilEnergy_eV[mask]
    recoil_nums = cascade_nums[mask] - 1 - i
    recoil_nVacs = cascade_nVacs[mask]
    recoil_atoms = cascades_atoms[mask]

    outputEvents.append({
      "primaryIndex": i,
      "energy_eV": primaryEnergy_eV,
      "pka_endpoint_x": pka_end_x[0],
      "pka_endpoint_y": pka_end_y[0],
      "pka_endpoint_z": pka_end_z[0],
      "xs_nm": xs.tolist(),
      "ys_nm": ys.tolist(),
      "zs_nm": zs.tolist(),
      "nVacs": recoil_nVacs.tolist(),
      "displacedAtoms_Z": recoil_atoms.tolist(),
      "recoilEnergies_eV": recoil_energies.tolist(),
      "recoilNums": recoil_nums.tolist()
    })

  #convert to dataframe
  return pd.DataFrame(outputEvents)

sym_to_Z = {
  "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
  "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
  "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
  "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
  "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
  "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
  "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
  "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
  "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
  "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
}

def processFastThrow(primaries_df,initialEnergy_eV,minEnergy,maxEnergy):
  outputEvents = []

  #Get primaries columns
  primaries_energy_eV = primaries_df["energy_keV"].to_numpy(copy=False)*1000
  primaries_recoilEnergy_eV = primaries_df["recoilEnergy_eV"].to_numpy(copy=False)
  primaries_xs_nm = primaries_df["x_A"].to_numpy(copy=False)*0.1
  primaries_ys_nm = primaries_df["y_A"].to_numpy(copy=False)*0.1
  primaries_zs_nm = primaries_df["z_A"].to_numpy(copy=False)*0.1
  primaries_nVacs = primaries_df["nVacs"].to_numpy(copy=False)
  primaries_atomHit = primaries_df["atom_hit"].map(sym_to_Z).to_numpy(dtype=np.uint8, copy=False)

  #Get the start positions of what will become each of our new primaries. Assume the first starts at (0,0,0), and each 
  #subsequent primary starts at the collision site of the previous primary
  primaries_start_xs_nm = np.empty_like(primaries_xs_nm)
  primaries_start_ys_nm = np.empty_like(primaries_ys_nm)
  primaries_start_zs_nm = np.empty_like(primaries_zs_nm)
  primaries_start_xs_nm[0] = 0.0
  primaries_start_ys_nm[0] = 0.0
  primaries_start_zs_nm[0] = 0.0
  if primaries_start_xs_nm.size > 1:
    primaries_start_xs_nm[1:] = primaries_xs_nm[:-1]
    primaries_start_ys_nm[1:] = primaries_ys_nm[:-1]
    primaries_start_zs_nm[1:] = primaries_zs_nm[:-1]
  
  #Calculate the direction it traveled, assuming a straight line between the start point and its collision point. Normalize
  vxs = primaries_xs_nm - primaries_start_xs_nm
  vys = primaries_ys_nm - primaries_start_ys_nm
  vzs = primaries_zs_nm - primaries_start_zs_nm
  #If the primary didn't move between collisions, assign same direction as previous track. This happens near the ends of tracks
  #very seldomly
  if vxs.size > 0:
    for i in range(0, vxs.size):
      if (vxs[i] == 0.0) and (vys[i] == 0.0) and (vzs[i] == 0.0):
        if i>0:
          vxs[i] = vxs[i-1]
          vys[i] = vys[i-1]
          vzs[i] = vzs[i-1]
        else:
          vxs[i] = 1.
          vys[i] = 0.
          vzs[i] = 0.
  norms = np.sqrt(vxs*vxs + vys*vys + vzs*vzs)
  primaries_vxs = vxs/norms
  primaries_vys = vys/norms
  primaries_vzs = vzs/norms
  
  #Calculate the energy each primary started with. In SRIM at each collision site we get the energy of the primary just prior
  #to the collision. To get the energy at the initial location it was generated, we need to take the energy of the previous
  #primary and subtract the recoil energy. The first primary starts out with initialEnergy_eV
  primaries_startEnergy_eV = np.empty_like(primaries_energy_eV)
  primaries_startEnergy_eV[0] = float(initialEnergy_eV)
  if primaries_startEnergy_eV.size > 1:
    primaries_startEnergy_eV[1:] = primaries_energy_eV[:-1] - primaries_recoilEnergy_eV[:-1]

  #Now step through primaries. Check if we should process it (ie is the bin already full we skip, or if the primary energy is
  #outside of our bin edges, skip). 
  for i in range(0,len(primaries_df)):
    #Check if this primary is out of bounds or the bin is already filled
    primaryEnergy_eV = float(primaries_startEnergy_eV[i])
    if primaryEnergy_eV < minEnergy or primaryEnergy_eV >= maxEnergy:
      continue

    #Get flat arrays of cascade positions, offset them so they start at 0,0,0
    xs = (primaries_xs_nm[i:] - primaries_start_xs_nm[i]).copy()
    ys = (primaries_ys_nm[i:] - primaries_start_ys_nm[i]).copy()
    zs = (primaries_zs_nm[i:] - primaries_start_zs_nm[i]).copy()

    #Compute rotation matrix to align these with the +X axis
    R = compute_rotation_matrix(primaries_vxs[i],primaries_vys[i],primaries_vzs[i])

    #Rotate all the cascades so they correspond to the +X axis
    rotate_points_inplace(xs, ys, zs, R)

    pka_end_x = xs[-1]
    pka_end_y = ys[-1]
    pka_end_z = zs[-1]

    #Get energies, nums
    recoil_energies = primaries_recoilEnergy_eV[i:]
    recoil_nums = np.array([j-i for j in range(i,len(primaries_df))])
    recoil_nVacs = primaries_nVacs[i:]
    recoil_atoms = primaries_atomHit[i:]

    outputEvents.append({
      "primaryIndex": i,
      "energy_eV": primaryEnergy_eV,
      "pka_endpoint_x": pka_end_x,
      "pka_endpoint_y": pka_end_y,
      "pka_endpoint_z": pka_end_z,
      "xs_nm": xs.tolist(),
      "ys_nm": ys.tolist(),
      "zs_nm": zs.tolist(),
      "nVacs": recoil_nVacs.tolist(),
      "displacedAtoms_Z": recoil_atoms.tolist(),
      "recoilEnergies_eV": recoil_energies.tolist(),
      "recoilNums": recoil_nums.tolist()
    })

  #convert to dataframe
  return pd.DataFrame(outputEvents)

@numba.njit(cache=True, fastmath=True)
def compute_rotation_matrix(dx, dy, dz):
  sinThetaCutoff = 1.74532e-7  # sin(1e-7 deg)
  
  # normalize input direction
  norm = (dx*dx + dy*dy + dz*dz) ** 0.5
  if norm == 0.0:
    return np.eye(3, dtype=np.float32)
  dx /= norm
  dy /= norm
  dz /= norm

  # axis = v x (1,0,0) = (0, dz, -dy)
  ax = 0.0
  ay = dz
  az = -dy

  sinTheta = (ay*ay + az*az) ** 0.5
  if sinTheta < sinThetaCutoff:  
    #Closely aligned with +X
    if dx > 0.0:
      return np.eye(3, dtype=np.float32)

    #Closely aligned with -X
    R = np.empty((3, 3), dtype=np.float32)
    R[0, 0] = -1.0; R[0, 1] =  0.0; R[0, 2] =  0.0
    R[1, 0] =  0.0; R[1, 1] =  1.0; R[1, 2] =  0.0
    R[2, 0] =  0.0; R[2, 1] =  0.0; R[2, 2] = -1.0
    return R

  # normalize axis
  ay /= sinTheta
  az /= sinTheta

  cosTheta = dx
  t = 1.0 - cosTheta
  s = sinTheta

  #Rodrigues rotation matrix
  R = np.empty((3, 3), dtype=np.float32)
  R[0, 0] = cosTheta + t*ax*ax
  R[0, 1] = t*ax*ay - s*az
  R[0, 2] = t*ax*az + s*ay

  R[1, 0] = t*ax*ay + s*az
  R[1, 1] = cosTheta + t*ay*ay
  R[1, 2] = t*ay*az - s*ax

  R[2, 0] = t*ax*az - s*ay
  R[2, 1] = t*ay*az + s*ax
  R[2, 2] = cosTheta + t*az*az
  return R

@numba.njit(cache=True, fastmath=True)
def rotate_points_inplace(x, y, z, R):
  # x,y,z are 1D float32 arrays, modified in-place
  for i in range(x.size):
    xi = x[i]
    yi = y[i]
    zi = z[i]
    x[i] = R[0, 0]*xi + R[0, 1]*yi + R[0, 2]*zi
    y[i] = R[1, 0]*xi + R[1, 1]*yi + R[1, 2]*zi
    z[i] = R[2, 0]*xi + R[2, 1]*yi + R[2, 2]*zi

def fillTree(tree,branches,output,bin_edges,counts_per_bin,maxEntriesPerBin):
  ionEnergy_eV, pka_end_x, pka_end_y, pka_end_z, xs, ys, zs, nVacs, displacedZs, recoilEnergies_eV, recoilNums = branches

  def clear_vectors():
    xs.clear()
    ys.clear()
    zs.clear()
    nVacs.clear()
    displacedZs.clear()
    recoilEnergies_eV.clear()
    recoilNums.clear()

  for df in output:
    if df is None or len(df) == 0:
      continue

    for row in df.itertuples(index=False):
      ionEnergy_eV[0] = float(row.energy_eV)

      bin_idx = np.searchsorted(bin_edges, ionEnergy_eV[0], side="right") - 1
      if bin_idx < 0 or bin_idx >= (len(bin_edges) - 1):
        continue
      if counts_per_bin[bin_idx] >= maxEntriesPerBin:
        continue

      pka_end_x[0] = float(row.pka_endpoint_x)
      pka_end_y[0] = float(row.pka_endpoint_y)
      pka_end_z[0] = float(row.pka_endpoint_z)
      
      clear_vectors()

      # positions
      for v in row.xs_nm:
        xs.push_back(float(v))
      for v in row.ys_nm:
        ys.push_back(float(v))
      for v in row.zs_nm:
        zs.push_back(float(v))

      for v in row.nVacs:
        nVacs.push_back(float(v)) 
      for v in row.displacedAtoms_Z:
        displacedZs.push_back(int(v)) 
      for v in row.recoilEnergies_eV:
        recoilEnergies_eV.push_back(float(v))
      for v in row.recoilNums:
        recoilNums.push_back(int(v))

      counts_per_bin[bin_idx]+=1

      tree.Fill()

def fillh5(dsets, events, bin_edges, counts_per_bin, maxEntriesPerBin):
  d_ionE, d_pka_end_x, d_pka_end_y, d_pka_end_z, d_xs, d_ys, d_zs, d_nv, d_disZ, d_recoilE, d_recoilNums = dsets

  to_write = []
  for ev in events:
    E = float(ev["energy_eV"])
    bin_idx = int(np.searchsorted(bin_edges, E, side="right") - 1)
    if bin_idx < 0 or bin_idx >= (len(bin_edges) - 1):
      continue
    if counts_per_bin[bin_idx] >= maxEntriesPerBin:
      continue
    to_write.append((E, ev))
    counts_per_bin[bin_idx] += 1

  num_new_events = len(to_write)
  if num_new_events == 0:
    return

  old_num_events = d_ionE.shape[0]
  num_events_total = old_num_events + num_new_events

  d_ionE.resize((num_events_total,))
  d_pka_end_x.resize((num_events_total,))
  d_pka_end_y.resize((num_events_total,))
  d_pka_end_z.resize((num_events_total,))
  d_xs.resize((num_events_total,))
  d_ys.resize((num_events_total,))
  d_zs.resize((num_events_total,))
  d_nv.resize((num_events_total,))
  d_disZ.resize((num_events_total,))
  d_recoilE.resize((num_events_total,))
  d_recoilNums.resize((num_events_total,))

  j = old_num_events
  for E, ev in to_write:
    d_ionE[j] = float(E)
    d_pka_end_x[j] = float(ev["pka_endpoint_x"])
    d_pka_end_y[j] = float(ev["pka_endpoint_y"])
    d_pka_end_z[j] = float(ev["pka_endpoint_z"])
    d_xs[j] = np.asarray(ev["xs_nm"], dtype=np.float32)
    d_ys[j] = np.asarray(ev["ys_nm"], dtype=np.float32)
    d_zs[j] = np.asarray(ev["zs_nm"], dtype=np.float32)
    d_nv[j] = np.asarray(ev["nVacs"], dtype=np.float32)
    d_disZ[j] = np.asarray(ev["displacedAtoms_Z"], dtype=np.uint8)
    d_recoilE[j] = np.asarray(ev["recoilEnergies_eV"], dtype=np.float32)
    d_recoilNums[j] = np.asarray(ev["recoilNums"], dtype=np.int16)

    j += 1