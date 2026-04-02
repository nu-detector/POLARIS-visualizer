# POLARIS Visualizer - Development Guide

## Project overview

3D rendering tool for the POLARIS underwater neutrino detector (Taiwan/Japan). Produces publication-quality illustrations using PyVista (VTK-based) with real ocean bathymetry.

## Key files

- `render_polaris.py` — main rendering script, single-file architecture
- `three_arm_polaris_fixed.geo` — detector geometry (62 strings, 10 DOMs each, 3-arm layout)
- `bathymetry_wide.nc` — SRTM15+ seafloor data from NOAA ERDDAP (east of Taiwan, 23.5N 123.5E)

## .geo file format

```
### Metadata ###
Medium: water
### Modules ###
x  y  z  string_id  dom_id
```

- Coordinates in meters, z is depth (negative = below sea level)
- 62 strings (id 0-61), 10 DOMs per string (id 0-9), 50m DOM spacing
- Z range: -2925 to -2475 m (450m instrumented depth)
- 3 arms radiating from center at ~120 degree separation
- Arm 1 (horizontal): y ~ +/-100m, x from -2500 to +2500, 500m between string pairs
- Arms 2 & 3: extend diagonally to y ~ +/-2215m

## Bathymetry data

Fetched from NOAA ERDDAP SRTM15+ service (15 arcsec resolution, ~464m/pixel):
```
https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus.nc?z[(LAT_MIN):1:(LAT_MAX)][(LON_MIN):1:(LON_MAX)]
```
- Returns NetCDF with variables: `z` (depth), `latitude`, `longitude`
- Read with `xarray.open_dataset()`
- Convert lat/lon to meters: `lat_m = (lat - center) * 111320`, `lon_m = (lon - center) * 111320 * cos(center_lat)`

## PyVista rendering notes

### Off-screen rendering
- Set `pv.OFF_SCREEN = True` before creating plotter for headless/batch rendering
- `pv.Plotter(off_screen=True)` also needed
- Use `plotter.screenshot(path)` to save

### Performance with many small objects
- Use `pv.PolyData(points).glyph(geom=pv.Sphere(radius=r))` to render many spheres efficiently from a point cloud, instead of creating individual sphere meshes
- `pv.MultiBlock()` to group line meshes

### Transparency
- `opacity=0.25` with `smooth_shading=True` for ghost-like terrain
- Terrain contour lines via `terrain.contour(isosurfaces=N)` add depth perception to transparent surfaces
- Custom colormaps via `matplotlib.colors.LinearSegmentedColormap.from_list()`, register with `matplotlib.colormaps.register()`

### Labels in 3D space
- `plotter.add_point_labels()` for labels anchored to 3D positions
- Use `shape=None, point_size=0, always_visible=True` for clean floating text
- `plotter.add_text()` for 2D overlay text (title)

### Camera
- Position is computed from spherical coordinates: distance, elevation (degrees from horizontal), azimuth
- Interactive mode: press `c` to print camera position, paste values back into script
- `plotter.camera.position`, `.focal_point`, `.up` control the view

### Terrain from bathymetry
- `scipy.ndimage.zoom(z, factor, order=3)` upsamples the grid for smoother rendering
- `scipy.interpolate.RegularGridInterpolator` to query seafloor height at arbitrary (x,y) — used for placing anchors on terrain
- Spatial compression (`x * scale_factor`) makes terrain features more dramatic relative to detector
- Vertical exaggeration (`centered * 1.5`) enhances relief

## Running

```bash
# Static high-res image
python render_polaris.py

# Interactive 3D window
python render_polaris.py --interactive

# Custom geometry
python render_polaris.py --geo other_geometry.geo --bathy other_bathy.nc -o output.png
```

## Design decisions

- Terrain is semi-transparent (opacity 0.25) with contour lines so detector strings have strong contrast
- Anchor blocks sit on interpolated seafloor surface; anchor-to-DOM lines vary in length per terrain
- DOM positions are never modified from the .geo file — only anchor positions adapt to terrain
- Taipei 101 (508m) used as scale reference — recognizable and regionally relevant
- Distance indicators (500m between string pairs, 200m within pairs) shown in red above detector
