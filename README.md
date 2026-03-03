# TRIM Python Utilities

Python tools for running **SRIM/TRIM** in batch mode and converting `COLLISON.txt` output into structured `.h5` or `.root` recoil track libraries.

Designed to efficiently generate reusable recoil-track libraries for downstream analysis.

---

# Repository Structure

```
trimRunner/
  trimRunner.py
  trimUtils.py
  mass_1.mas20.txt
  config/
    annotated.jsonc

collisionParser/
  collisionParser.py
  parserUtils.py
```

- `trimRunner` — Generates TRIM simulations
- `collisionParser` — Parses TRIM output into track libraries

---

# Part I — Running TRIM

## trimRunner.py

Runs SRIM/TRIM in batch mode using a JSON config file.

Each run:

- Creates a temporary SRIM working directory
- Writes `TRIM.in`
- Executes `TRIM.exe`
- Compresses `COLLISON.txt`
- Cleans up temporary files automatically

---

## Usage

From inside `trimRunner/`:

```bash
python trimRunner.py config/annotated.jsonc
```

## Required SRIM Setup

Before running:

1. Open `TRIMAUTO` inside your SRIM installation.
2. Change the first line from:

```
0
```

to:

```
1
```

This enables batch execution. Note that launching TRIM via the SRIM interface resets this to 1, so you will need to set it back to 0 before running this script if you have used the SRIM GUI.

---

## Materials

Materials are defined in:

```
trimRunner/trimUtils.py
```

Currently implemented:

- `LiF`
- `Olivine`
- `Diamond`

To add a new material, modify `createMaterialsDict()` in `trimUtils.py`.

---

## Config files

An annotated example configuration file is provided:

```
trimRunner/config/annotated.jsonc
```

This shows all required keys and expected structure.

---

### Required Config Fields

Your config must include:

| Key | Description |
|------|------------|
| `material` | Must exist in materials dictionary |
| `ionSymbol` | e.g. 208Pb, 16O, etc. |
| `runMode` | `"damage"` for tracks or `"efficiency" for vacancy-production efficiency |
| `calcMode` | `"quick"` or `"full"` |
| `nps` | Number of primary ions |
| `energy_keV` | Ion energy |
| `outputPath` | Directory to store output files|
| `SRIM_TEMP_FOLDER` | Temporary working directory |
| `SRIM_FOLDER` | Path to SRIM installation |

---

### Notes
* In `damage` mode, output format is <material>_<ion>_<energy>keV.tar.gz. This is just a zipped ASCII file of the standard TRIM output
* In `efficiency` mode, output format is <material>_<ion>_efficiency.csv.
* In `efficiency` mode, energy_keV should be an energy where efficiency==1, typically a few keV is sufficient.
* In `damage` mode, energy_keV should be several hundred keV greater than the largest recoil you care about, allowing for "burn in"
  
---

# Part II — Parsing TRIM Output

## collisionParser.py

Converts `COLLISON.txt` (or `.gz`) into:

- `.h5`
- `.root`

---

## Usage

From inside `collisionParser/`:

```bash
python collisionParser.py <input.txt|.gz> <output.h5|.root> <fast|full>
```

Example:

```bash
python collisionParser.py LiF_16O_200keV.tar.gz tracks.h5 full
```

---

## Fast vs Full Modes

| Mode | Description |
|-------|------------|
| `full` | Parses full cascade output |
| `fast` | Parses quick/summary format |

The parser mode must match how TRIM was run.

---

## Artificial Energy Binning

TRIM output is not energy-binned.  
The parser imposes artificial binning to:

- Limit file size
- Prevent low-energy recoils from dominating storage

Key parameters (inside `collisionParser.py`):

| Parameter | Example | Description |
|------------|--------------|-------------|
| `energyRange_eV` | `[1, 300e3]` | Energy window (in eV) of recoils to retain |
| `spacing` | `"linear"` | Energy bin spacing (`"linear"` or `"log"`) |
| `nBins` | `100` | Number of artificial energy bins |
| `nCores` | `8` | Number of cores to use for processing tracks |
| `maxCascadesToHold` | `1000000` | only relevant in 'full' mode. 1e6 cascades ~5GB peak RAM usage |
| `monoMode` | `False` | For not doing the downstream ion track generation, i.e. pure mono output|

---

## Output Structure

Each stored event contains:

- `ionEnergy_eV`
- `pka_endpoint_x_nm`
- `pka_endpoint_y_nm`
- `pka_endpoint_z_nm`
- `xs_nm`
- `ys_nm`
- `zs_nm`
- `nVacs`
- `displacedAtoms_Z`
- `recoilEnergies_eV`
- `recoilNums`

Coordinates are rotated so that the primary recoil direction aligns with +X.

---

# Dependencies

- Python 3
- SRIM/TRIM (installed separately)
- numpy
- pandas
- h5py (for `.h5` output)
- ROOT (for `.root` output)
- numba
- tqdm

---

# Design Philosophy

Instead of running thousands of mono-energetic TRIM simulations:

1. Simulate primaries slightly above your maximum energy of interest.
2. Extract lower-energy recoils from secondary cascades.
3. Build a reusable recoil-track library.
4. Apply detector efficiency separately in downstream analysis.

This dramatically reduces required TRIM runtime.

---

# Notes

- Efficiency curves are not embedded in the track library.
