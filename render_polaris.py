#!/usr/bin/env python3
"""Render a publication-quality 3D illustration of the Polaris neutrino detector."""

import argparse
import os
import numpy as np
from collections import defaultdict

import pyvista as pv
import requests
import xarray as xr
from scipy.ndimage import zoom as ndzoom
from scipy.interpolate import RegularGridInterpolator

# Default bathymetry region: east of Taiwan (23.15-23.85N, 123.15-123.85E)
BATHY_LAT_CENTER = 23.5
BATHY_LON_CENTER = 123.5
BATHY_HALF_DEG = 0.35


def download_bathymetry(output_path):
    """Download SRTM15+ bathymetry from NOAA ERDDAP."""
    lat_min = BATHY_LAT_CENTER - BATHY_HALF_DEG
    lat_max = BATHY_LAT_CENTER + BATHY_HALF_DEG
    lon_min = BATHY_LON_CENTER - BATHY_HALF_DEG
    lon_max = BATHY_LON_CENTER + BATHY_HALF_DEG
    url = (
        f"https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus.nc"
        f"?z[({lat_min}):1:({lat_max})][({lon_min}):1:({lon_max})]"
    )
    print(f"Downloading bathymetry from NOAA ERDDAP...")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    with open(output_path, 'wb') as f:
        f.write(resp.content)
    print(f"Saved bathymetry to {output_path} ({len(resp.content)} bytes)")


def parse_geo_file(path):
    """Parse .geo file and return dict of strings with their DOM positions."""
    strings = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('Medium') or not line:
                continue
            parts = line.split()
            if len(parts) == 5:
                x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                sid, did = int(parts[3]), int(parts[4])
                strings[sid].append((x, y, z, did))
    for sid in strings:
        strings[sid].sort(key=lambda t: t[3])
    return strings


def get_string_positions(strings):
    """Get unique (x, y) position for each string."""
    positions = {}
    for sid, doms in strings.items():
        positions[sid] = (doms[0][0], doms[0][1])
    return positions


def load_bathymetry(nc_path, center_lat, center_lon):
    """Load bathymetry NetCDF and convert to meters relative to center."""
    ds = xr.open_dataset(nc_path)
    z = ds['z'].values
    lats = ds.latitude.values
    lons = ds.longitude.values
    y_m = (lats - center_lat) * 111320
    x_m = (lons - center_lon) * 111320 * np.cos(np.radians(center_lat))
    ds.close()
    return x_m, y_m, z


def build_terrain_mesh(x_m, y_m, z, upsample=4):
    """Build a PyVista StructuredGrid from bathymetry, upsampled for smoothness."""
    # Upsample with cubic interpolation for smoother terrain
    z_fine = ndzoom(z, upsample, order=3)
    x_fine = np.linspace(x_m[0], x_m[-1], z_fine.shape[1])
    y_fine = np.linspace(y_m[0], y_m[-1], z_fine.shape[0])
    X, Y = np.meshgrid(x_fine, y_fine)
    grid = pv.StructuredGrid(X, Y, z_fine)
    return grid


def build_taipei101(base_center, height=508):
    """Build a simplified Taipei 101 tower as stacked tapered boxes + spire."""
    x0, y0, z0 = base_center
    meshes = pv.MultiBlock()

    # Base podium
    podium = pv.Box(bounds=(x0 - 40, x0 + 40, y0 - 40, y0 + 40, z0, z0 + 25))
    meshes.append(podium)

    # 8 tapered sections (the iconic design)
    section_height = 50
    z_cur = z0 + 25
    width_bottom = 35
    taper = 2.5  # each section narrows at bottom, widens at top

    for i in range(8):
        w_bot = width_bottom - i * 1.5
        w_top = w_bot + taper
        # Approximate with a box at average width
        w_avg = (w_bot + w_top) / 2
        box = pv.Box(bounds=(
            x0 - w_avg, x0 + w_avg,
            y0 - w_avg, y0 + w_avg,
            z_cur, z_cur + section_height
        ))
        meshes.append(box)
        z_cur += section_height

    # Upper narrowing section
    for i in range(3):
        w = 15 - i * 4
        h = 20
        box = pv.Box(bounds=(
            x0 - w, x0 + w,
            y0 - w, y0 + w,
            z_cur, z_cur + h
        ))
        meshes.append(box)
        z_cur += h

    # Spire
    spire = pv.Line((x0, y0, z_cur), (x0, y0, z0 + height))
    meshes.append(spire)

    return meshes


def build_scale_bar(origin, length, direction, tick_height=30):
    """Build a scale bar as line segments."""
    end = origin + np.array(direction) * length
    bar = pv.Line(origin, end)
    tick1 = pv.Line(origin - np.array([0, 0, tick_height]),
                    origin + np.array([0, 0, tick_height]))
    tick2 = pv.Line(end - np.array([0, 0, tick_height]),
                    end + np.array([0, 0, tick_height]))
    return bar + tick1 + tick2


def render(geo_path, bathy_path, output_path, resolution=(4000, 3000),
           interactive=False, animate=False):
    pv.OFF_SCREEN = not interactive

    if not os.path.exists(bathy_path):
        download_bathymetry(bathy_path)

    strings = parse_geo_file(geo_path)
    positions = get_string_positions(strings)

    all_z = [z for doms in strings.values() for (_, _, z, _) in doms]
    z_min, z_max = min(all_z), max(all_z)
    all_xy = np.array(list(positions.values()))
    cx, cy = all_xy.mean(axis=0)

    # --- Bathymetry ---
    x_m, y_m, bathy_z = load_bathymetry(bathy_path, BATHY_LAT_CENTER, BATHY_LON_CENTER)

    # Shift so the *highest* point of the seafloor sits 150m below detector bottom
    bathy_max = np.nanmax(bathy_z)
    seafloor_target = z_min - 150
    bathy_z_shifted = bathy_z - bathy_max + seafloor_target

    # Compress spatial extent to make terrain more dramatic relative to detector
    spatial_scale = 0.5
    x_scaled = x_m * spatial_scale
    y_scaled = y_m * spatial_scale

    # Exaggerate vertical relief for visual impact
    bathy_centered = bathy_z_shifted - np.nanmean(bathy_z_shifted)
    bathy_exaggerated = np.nanmean(bathy_z_shifted) + bathy_centered * 1.5

    terrain = build_terrain_mesh(x_scaled, y_scaled, bathy_exaggerated, upsample=4)

    # Build interpolator for seafloor height at any (x, y) position
    seafloor_interp = RegularGridInterpolator(
        (y_scaled, x_scaled), bathy_exaggerated,
        method='linear', bounds_error=False,
        fill_value=np.nanmean(bathy_exaggerated)
    )

    # --- Detector ---
    dom_points = []
    for doms in strings.values():
        for x, y, z, _ in doms:
            dom_points.append([x, y, z])
    dom_points = np.array(dom_points)
    dom_cloud = pv.PolyData(dom_points)

    string_lines = pv.MultiBlock()
    anchor_lines = pv.MultiBlock()  # lines from seafloor anchor to first DOM
    anchor_blocks = []
    for doms in strings.values():
        top = np.array([doms[-1][0], doms[-1][1], doms[-1][2]])
        bot = np.array([doms[0][0], doms[0][1], doms[0][2]])
        # Instrumented string line (between DOMs)
        string_lines.append(pv.Line(bot, top))
        # Interpolate seafloor z at this string's (x, y)
        floor_z = float(seafloor_interp((bot[1], bot[0])))
        # Anchor block sitting on seafloor
        bsz = 10
        anchor = pv.Box(bounds=(
            bot[0] - bsz, bot[0] + bsz,
            bot[1] - bsz, bot[1] + bsz,
            floor_z - bsz, floor_z + bsz
        ))
        anchor_blocks.append(anchor)
        # Line from anchor top to bottom DOM (variable length)
        anchor_top = np.array([bot[0], bot[1], floor_z + bsz])
        anchor_lines.append(pv.Line(anchor_top, bot))

    # --- Taipei 101 for scale --- next to central strings, on the ground
    t101_x = 1800
    t101_y = 400  # just beside the center string pair
    t101_z = float(seafloor_interp((t101_y, t101_x))) + 30  # on the ground
    taipei101 = build_taipei101((t101_x, t101_y, t101_z), height=508)

    # --- Scale bar (1 km) --- placed close to detector edge
    bar_origin = np.array([cx - 500, all_xy[:, 1].min() - 300, z_min])
    scale_bar = build_scale_bar(bar_origin, 1000, [1, 0, 0], tick_height=50)

    # --- 500m indicator between two adjacent string pairs along arm 1 ---
    # Arm 1 strings at y≈±100: e.g. string at x=0 and x=500 (500m apart)
    ind500_z = z_max + 60
    ind500_start = np.array([1000, -100, ind500_z])
    ind500_end = np.array([1500, -100, ind500_z])
    ind500_bar = build_scale_bar(ind500_start, 500, [1, 0, 0], tick_height=30)

    # --- 200m indicator between two strings within a pair ---
    # e.g. string pair at x=0: y=-100 and y=100
    ind200_z = z_max + 60
    ind200_start = np.array([1000, -100, ind200_z])
    ind200_end = np.array([1000, 100, ind200_z])
    ind200_bar = build_scale_bar(ind200_start, 200, [0, 1, 0], tick_height=30)

    # --- Vertical indicator over instrumented string length (450m) ---
    # Place next to an outer string on arm 1 for visibility
    indv_x = 2500 + 150  # just outside the rightmost string
    indv_y = -100
    indv_bar = build_scale_bar(
        np.array([indv_x, indv_y, z_min]), int(z_max - z_min),
        [0, 0, 1], tick_height=30
    )

    # --- Render ---
    plotter = pv.Plotter(window_size=resolution,
                         off_screen=not interactive)
    plotter.set_background('white', top='#D6E8F7')

    # Terrain - mostly transparent with contour-like shading
    from matplotlib.colors import LinearSegmentedColormap
    # Go from nearly transparent light tone to darker lines at depth changes
    seafloor_colors = [
        '#1A1A2E',  # dark (valleys) - will show through transparency
        '#3D3D5C',  # medium dark
        '#6B6B8A',  # muted blue-gray
        '#9E9EB5',  # lighter
        '#C8C8D8',  # pale (ridges)
    ]
    seafloor_cmap = LinearSegmentedColormap.from_list('seafloor', seafloor_colors, N=256)
    import matplotlib
    matplotlib.colormaps.register(cmap=seafloor_cmap, name='seafloor')

    plotter.add_mesh(terrain, scalars=terrain.points[:, 2],
                     cmap='seafloor',
                     clim=[bathy_exaggerated.min(), bathy_exaggerated.max()],
                     show_scalar_bar=False, opacity=0.25,
                     smooth_shading=True)

    # Add wireframe contour lines on terrain for depth visualization
    contours = terrain.contour(isosurfaces=15, scalars=terrain.points[:, 2])
    if contours.n_points > 0:
        plotter.add_mesh(contours, color='#2A2A4A', line_width=1.2, opacity=0.35)

    # Diffuse lighting
    light = pv.Light(position=(5000, -5000, 5000), focal_point=(0, 0, -2700),
                     color='white', intensity=0.8)
    light.positional = False  # directional / diffuse
    plotter.add_light(light)
    plotter.enable_lightkit()  # adds soft ambient fill

    # Anchor blocks on seafloor
    for ab in anchor_blocks:
        plotter.add_mesh(ab, color='#444444', opacity=0.9, smooth_shading=True)

    # Anchor lines (seafloor to first DOM, variable length)
    for al in anchor_lines:
        plotter.add_mesh(al, color='#555555', line_width=1.5, opacity=0.6)

    # String lines (instrumented section between DOMs)
    for block in string_lines:
        plotter.add_mesh(block, color='#222222', line_width=2.5, opacity=0.85)

    # DOMs
    dom_glyphs = dom_cloud.glyph(geom=pv.Sphere(radius=12),
                                  scale=False, orient=False)
    plotter.add_mesh(dom_glyphs, color='#D4A017', opacity=0.9,
                     smooth_shading=True)

    # Taipei 101
    for block in taipei101:
        if isinstance(block, pv.PolyData) and block.n_cells > 0:
            if block.n_points == 2:  # spire line
                plotter.add_mesh(block, color='#555555', line_width=2)
            else:
                plotter.add_mesh(block, color='#708090', opacity=0.85,
                                 smooth_shading=True)
    # Taipei 101 label
    t101_label_pt = np.array([t101_x, t101_y, t101_z + 508 + 80])
    plotter.add_point_labels(
        pv.PolyData(t101_label_pt.reshape(1, 3)),
        ['Taipei 101\n(508 m)'], font_size=52, text_color='#333333',
        shape=None, render_points_as_spheres=False,
        point_size=0, always_visible=True,
    )

    # 1 km scale bar
    plotter.add_mesh(scale_bar, color='black', line_width=4)
    bar_mid = bar_origin + np.array([500, -120, 0])
    plotter.add_point_labels(
        pv.PolyData(bar_mid.reshape(1, 3)),
        ['1 km'], font_size=60, text_color='black',
        shape=None, render_points_as_spheres=False,
        point_size=0, always_visible=True,
    )

    # 500m indicator
    plotter.add_mesh(ind500_bar, color='#AA0000', line_width=3)
    ind500_mid = (ind500_start + ind500_end) / 2 + np.array([0, -80, 0])
    plotter.add_point_labels(
        pv.PolyData(ind500_mid.reshape(1, 3)),
        ['500 m'], font_size=56, text_color='#AA0000',
        shape=None, render_points_as_spheres=False,
        point_size=0, always_visible=True,
    )

    # 200m indicator
    plotter.add_mesh(ind200_bar, color='#AA0000', line_width=3)
    ind200_mid = (ind200_start + ind200_end) / 2 + np.array([-120, 0, 0])
    plotter.add_point_labels(
        pv.PolyData(ind200_mid.reshape(1, 3)),
        ['200 m'], font_size=56, text_color='#AA0000',
        shape=None, render_points_as_spheres=False,
        point_size=0, always_visible=True,
    )

    # Vertical string length indicator (450m)
    plotter.add_mesh(indv_bar, color='#AA0000', line_width=3)
    indv_mid = np.array([indv_x + 100, indv_y, (z_min + z_max) / 2])
    plotter.add_point_labels(
        pv.PolyData(indv_mid.reshape(1, 3)),
        [f'{int(z_max - z_min)} m'], font_size=56, text_color='#AA0000',
        shape=None, render_points_as_spheres=False,
        point_size=0, always_visible=True,
    )

    # Title
    plotter.add_text('POLARIS', position=(0.88, 0.05), font_size=48,
                     color='black', font='times', viewport=True)

    # Camera - higher perspective, focal point lowered to include anchors
    cam_dist = 10000
    cam_elev = 45
    cam_azim = -40
    focal_z = z_min - 100  # lower focal to include seafloor anchors
    elev_rad = np.radians(cam_elev)
    azim_rad = np.radians(cam_azim)
    cam_x = cx + cam_dist * np.cos(elev_rad) * np.cos(azim_rad)
    cam_y = cy + cam_dist * np.cos(elev_rad) * np.sin(azim_rad)
    cam_z = focal_z + cam_dist * np.sin(elev_rad)
    focal = np.array([cx, cy, focal_z])

    # Move camera 20% closer to focal point
    cam_pos = np.array([3003.409, 6783.348, 1651.120])
    foc_pos = np.array([266.620, 503.765, -2972.534])
    cam_pos = cam_pos + 0.20 * (foc_pos - cam_pos)

    plotter.camera.position = tuple(cam_pos)
    plotter.camera.focal_point = tuple(foc_pos)
    plotter.camera.up = (-0.19071551529371714, -0.5267804585977117, 0.8283296087101055)

    plotter.enable_anti_aliasing('msaa')

    if animate:
        # Drone-like orbit around the 3-arm intersection (detector center)
        orbit_center = np.array([0.0, 0.0, (z_min + z_max) / 2])
        orbit_radius = 9000  # distance from center in xy plane
        orbit_elev = 45  # degrees from horizontal
        orbit_z_offset = orbit_radius * np.sin(np.radians(orbit_elev))
        orbit_r_xy = orbit_radius * np.cos(np.radians(orbit_elev))
        start_angle = np.arctan2(cam_pos[1] - orbit_center[1],
                                  cam_pos[0] - orbit_center[0])

        n_frames = 720
        anim_output = output_path.rsplit('.', 1)[0] + '.mp4'

        plotter.open_movie(anim_output, framerate=30)
        print(f"Recording {n_frames} frames to {anim_output}...")

        for i in range(n_frames):
            angle = start_angle + 2 * np.pi * i / n_frames
            # Gentle vertical bobbing like a drone
            z_bob = orbit_z_offset + 150 * np.sin(2 * np.pi * i / n_frames)
            cx_i = orbit_center[0] + orbit_r_xy * np.cos(angle)
            cy_i = orbit_center[1] + orbit_r_xy * np.sin(angle)
            cz_i = orbit_center[2] + z_bob
            plotter.camera.position = (cx_i, cy_i, cz_i)
            # Focal point slightly below orbit center so camera looks more downward
            plotter.camera.focal_point = (orbit_center[0], orbit_center[1],
                                          orbit_center[2] - 500)
            plotter.camera.up = (0, 0, 1)
            plotter.write_frame()

        plotter.close()
        print(f'Saved animation to {anim_output}')
    elif interactive:
        print("Interactive mode — press 'c' to print camera position, 'q' to quit")
        plotter.show()
        # Print final camera position for copying back into the script
        pos = plotter.camera.position
        foc = plotter.camera.focal_point
        up = plotter.camera.up
        print(f"\n# Camera position to paste into script:")
        print(f"plotter.camera.position = {pos}")
        print(f"plotter.camera.focal_point = {foc}")
        print(f"plotter.camera.up = {up}")
    else:
        plotter.screenshot(output_path)
        plotter.close()
        print(f'Saved rendering to {output_path}')


def main():
    parser = argparse.ArgumentParser(description='Render Polaris detector')
    parser.add_argument('--geo', default='three_arm_polaris_fixed.geo',
                        help='Path to detector geometry file')
    parser.add_argument('--bathy', default='bathymetry_wide.nc',
                        help='Path to bathymetry NetCDF file')
    parser.add_argument('--output', '-o', default='polaris_detector.png',
                        help='Output image path (default: polaris_detector.png)')
    parser.add_argument('--width', type=int, default=4000,
                        help='Image width in pixels (default: 4000)')
    parser.add_argument('--height', type=int, default=2400,
                        help='Image height in pixels (default: 2400)')
    parser.add_argument('--interactive', '-i', action='store_true',
                        help='Open interactive 3D window')
    parser.add_argument('--animate', '-a', action='store_true',
                        help='Record drone orbit animation (MP4)')
    args = parser.parse_args()

    render(args.geo, args.bathy, args.output,
           resolution=(args.width, args.height),
           interactive=args.interactive, animate=args.animate)


if __name__ == '__main__':
    main()
