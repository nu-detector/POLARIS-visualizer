# POLARIS Visualizer

Publication-quality 3D rendering of the POLARIS underwater neutrino detector.

## Installation

```bash
pip install -e .
```

Or install dependencies directly:

```bash
pip install numpy scipy pyvista xarray netCDF4 matplotlib
```

## Usage

### Generate a static image

```bash
python render_polaris.py
```

This produces `polaris_detector.png` (4000x3000 px) suitable for journal publication.

### Options

```
--geo FILE          Detector geometry file (default: three_arm_polaris_fixed.geo)
--bathy FILE        Bathymetry NetCDF file (default: bathymetry_wide.nc)
-o, --output FILE   Output image path (default: polaris_detector.png)
--width N           Image width in pixels (default: 4000)
--height N          Image height in pixels (default: 3000)
-i, --interactive   Open interactive 3D window instead of saving image
```

### Interactive mode

```bash
python render_polaris.py --interactive
```

Opens a 3D window where you can explore the scene.

#### Mouse controls

| Action | Control |
|--------|---------|
| Rotate | Left mouse drag |
| Pan / translate | Middle mouse drag, or Shift + left drag |
| Zoom | Scroll wheel, or right mouse drag |

#### Keyboard shortcuts

| Key | Action |
|-----|--------|
| `c` | Print current camera position to terminal |
| `r` | Reset camera to fit all objects |
| `q` | Quit (prints final camera position for reuse) |

When you close the interactive window, the final camera position is printed to the terminal. You can paste these values into the script to reproduce the exact same view in batch mode.

## Input files

- **`.geo` file** -- Detector geometry with one line per optical module: `x y z string_id dom_id`. Coordinates in meters.
- **`bathymetry_wide.nc`** -- SRTM15+ ocean bathymetry (NetCDF). Fetched from NOAA ERDDAP for the region east of Taiwan.
