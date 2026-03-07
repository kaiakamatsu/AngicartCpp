# Angicart

Angicart is an efficient software for vascular network segmentation and analysis, designed and maintained by the [Savage Lab](https://savagelab.eeb.ucla.edu/) at UCLA.

## References

The original algorithm for segmentation is described in:

> Newberry MG, Ennis DB, Savage VM (2015) Testing Foundations of Biological Scaling Theory Using Automated Measurements of Vascular Networks. *PLoS Comput Biol* 11(8): e1004455. https://doi.org/10.1371/journal.pcbi.1004455

Applications to PET–CT scan data from cancer patients are described in:

> Brummer AB, Savage VM (2021) Cancer as a Model System for Testing Metabolic Scaling Theory. *Front. Ecol. Evol.* https://doi.org/10.3389/fevo.2021.691830

The method for quantifying curvature and tortuosity is described in:

> Brummer AB, Hunt D, Savage VM (2021) Improving Blood Vessel Tortuosity Measurements via Highly Sampled Numerical Integration of the Frenet–Serret Equations. *IEEE Trans. Med. Imaging* https://doi.org/10.1109/TMI.2020.3025467

---

## Quick installation

> Note that the source code includes an optional GUI build using wxWidgets. For our installation instructions and executable files, we exclude the GUI for faster processing. If you want the GUI, you would have to (1) install wxWidgets, (2) set `TRY_WX` to 1 in the first line of `source/minSurfTests.cpp` and (3) rebuild. 

### Windows

A prebuilt x64 executable can be found in releases: https://github.com/kaiakamatsu/AngicartCpp/releases/tag/v1.0.0windowsmac

**Build from source (Visual Studio, not VS Code):**

1. Install [Visual Studio 2022 Community](https://visualstudio.microsoft.com/downloads/) (with C++ workload).
2. Open this repository in Visual Studio.
3. Set the top toolbar to **x64-Release**.
4. In Solution Explorer, open `angicarttest.sln`.
5. Right-click the solution → **Properties** → **Configuration Properties** (set Configuration to *All configurations*) → **Linker** → **System** → set **Subsystem** to **Console (/SUBSYSTEM:CONSOLE)**.
6. **Build** → **Build angicarttest.sln**. You are now ready to run the executable or directly build and run from Visual Studio.

### macOS

A prebuilt macOS executable can be found in releases: https://github.com/kaiakamatsu/AngicartCpp/releases/tag/v1.0.0windowsmac

**Build from source (Xcode):**

1. Install the Xcode command-line tools:
   ```bash
   xcode-select --install
   ```
2. Open the folder on xcode. You are now ready to build and run from xcode.

> We noticed that on some macOS devices, users encounter privacy/security warnings when trying to run the executable file. Please see `other_resources/Angicart Notes for MacOS Users.pdf` for a simple fix.

---

## Quick run on example images, inputs and outputs

Example input images are in `sample_data/FITC-MCA0_N12_PI001_s1/`. Each `.png` file is a 2D slice of the 3D volume; Angicart reconstructs the volume from these slices before analysis.

**Command-line usage (10 or 11 arguments):**

```
angicart.exe <image_dir> <im_start> <im_end> <output_base> <vox_x> <vox_y> <vox_z> <length_unit> <threshold> [thread_count]
```

| Argument      | Description |
|---------------|-------------|
| `image_dir`   | Directory containing the input slice images WITHOUT the last '/' (e.g. `sample_data/FITC-MCA0_N12_PI001_s1`). |
| `im_start`    | First slice index to include (inclusive). |
| `im_end`      | Last slice index to include (inclusive). |
| `output_base` | Base path for output files (`_vessels.png`, `.tsv`, `_withRoots.tsv`, `.dat`, etc. are appended to this prefix). |
| `vox_x`       | Voxel size in the x direction (same units as `length_unit`). |
| `vox_y`       | Voxel size in the y direction (same units as `length_unit`). |
| `vox_z`       | Voxel size in the z direction (same units as `length_unit`). |
| `length_unit` | Name of the length unit for voxel dimensions and output (e.g. `um` or `mm`). |
| `threshold`   | Normalized intensity threshold (0–1) to separate vascular from non-vascular voxels. |
| `thread_count`| *(Optional)* Number of threads. If omitted, defaults to (CPU cores − 1). |

> CRITICAL: `threshold` is the only hyperparameter that has to be tuned, usually for each tissue-dataset pair. Before running Angicart, it helps to explore different threshold values visually to identify an optimal normalized treshold value for identifying candidate vascular voxels. 

**Example run** (from the repository root, or adjust paths):

```bash
angicart.exe sample_data/FITC-MCA0_N12_PI001_s1 0 61 sample_outputs/FITC-MCA0_N12_PI001_s1 0.412 0.412 0.492 um 0.2
```

Optional 11th argument: thread count (e.g. `8`). If omitted, the program uses a default based on your CPU.

All outputs use the path given as `output_base` with the following suffixes:

| Output file | Description |
|-------------|-------------|
| `output_base.tsv` | Summary table per backbone: name, volume, length, radius from volume/length, observed mean radius per vessel, number of adjacent segments, and adjacency list. |
| `output_base_withRoots.tsv` | Same as above but with a rooted tree: each segment has a parent and children instead of a flat adjacency list. |
| `output_base.dat` | Per-backbone tables: for each segment, vertebra index, voxel coordinates (x, y, z), and per-vertebra radii (r_solid, r_surf). Blocks separated by blank lines; each block has a header with the segment name (same as in the TSV). |
| `output_base_vessels.png` | Sanity-check visualization of the binary vessel volume before backbone extraction. |

Each segment is labeled (name column) such that the index of the backbone within its connected component appears before the decimal point, and the connected component label appears after the  decimal point.

## Algorithm details

Pipeline: intensity threshold → binary volume → sphere coarsening → adjacency and false-tip removal → critical vertebrae → backbone erosion → segmentation and tip extension → segregation by connected component.

**1. Input and preprocessing**
- Load 2D slices into a 3D volume; threshold intensities to obtain a binary vessel mask.
- Keep the largest connected components; fill interior holes.

**2. Sphere coarsening** (`sphereCoarsen`, `critFrac` ≈ 0.16)
- Build disjoint spheres that roughly cover the vasculature. Spheres are grown in voxel space (not physical dimensions).
- For each candidate seed: grow a sphere while the fraction of vascular voxels in the expanding shell stays above `critFrac`. Stop if the shell becomes too non-vascular.
- If the growing sphere hits a voxel already in another sphere, reject the seed; otherwise the sphere is valid.
- Valid sphere centers may be repositioned to the center of mass of their vascular region until the center and radius stabilize.
- Continue until a specified number of consecutive candidate seeds fail.
- Once spheres are determined, adjacency is computed: which intersphere regions connect which intrasphere (sphere) seeds.

**3. Connectivity and false-tip removal**
- **Intrasphere space:** voxels within each sphere's radius. **Intersphere space:** voxels in the complement that form connections between spheres.
- Remove false tips: drop spheres with only one connection and a singly connected intersphere component; drop any intersphere region touching fewer than two spheres (don't worry, tips will be extended later)

**4. Critical vertebrae and raw backbones**
- **Critical vertebrae** = center of mass of each singly connected component at the interface between intrasphere and intersphere space. They define the tips for skeletonization.
- **Backbone erosion:** from these tips, erode the vessel to obtain centerlines (ordered voxel sequences) farthest from the boundary. One raw backbone per sphere or connection.

**5. Segmentation and tip extension**
- Order and stitch raw backbones using critical vertebrae; split at branch points so each segment has its own backbone. Remove empty backbones.
- Assign vessel voxels to the nearest backbone (meat per segment). For tip segments, extend the backbone into its meat (erosion from current tip toward farthest voxel in the segment). Compute per-segment volume and average radius (and per-centerline voxel radius).

**6. Output**
- Backbones are segregated by connected components. Results are written to the TSV, rooted TSV, and `.dat` files.


## Contributors 

This repository is currently being maintained by Kai Akamatsu (kakamatsu@g.ucla.edu) and Dr. Van Savage (vsavage@g.ucla.edu).

- Dr. Van Savage
- Jocelyn J. Shen
- Mitchell Newberry
- Dr. David Hunt
- Dr. Alex Brummer
- Anderson Ju
- Kai Akamatsu

Our software is tested across multiple devices thanks to talented undergraduate researchers:

- Kaili Tam
- David Zhang
- Charissa Mak
- Josh Mehdian

