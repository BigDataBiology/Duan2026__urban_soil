"""Figure 1a: map of Jiangsu + Shanghai (Nanjing highlighted), drawn with geopandas.

Data sources (downloaded on first run into ./geo_data):
  * Jiangsu prefectures and Shanghai districts: Aliyun DataV GeoAtlas
  * Neighbouring provinces (grey background): Natural Earth 1:50m admin-1

All geometry is projected to Web Mercator (EPSG:3857) and then mapped onto the
final canvas, so every coordinate, tolerance and line width below is in
millimetres (line widths: points) of the W_MM x H_MM output. The view is fitted
to Jiangsu + Shanghai, leaving MARGIN_MM on the tightest side.
"""
from pathlib import Path
import urllib.request

import geopandas as gpd
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['svg.fonttype'] = 'none'  # keep text editable in the SVG
from shapely.affinity import affine_transform
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import linemerge, nearest_points, unary_union

HERE = Path(__file__).resolve().parent
DATA = HERE / 'geo_data'

URLS = {
    'datav_320000_full.json': 'https://geo.datav.aliyun.com/areas_v3/bound/320000_full.json',
    'datav_310000_full.json': 'https://geo.datav.aliyun.com/areas_v3/bound/310000_full.json',
    'ne_50m_admin_1.zip': 'https://naciscdn.org/naturalearth/50m/cultural/ne_50m_admin_1_states_provinces.zip',
}

# Output size
W_MM, H_MM = 66, 62
MM_PER_INCH = 25.4
DPI = 300  # for the PNG only
MARGIN_MM = 1.5  # space between Jiangsu/Shanghai and the edge (on the tightest side)

# Colours
BG_LAND = '#f2f2f2'
BG_EDGE = '#c8c8c8'
JS_FILL = '#f5f5f5'
NANJING = '#1a9c76'
SH_FILL = '#eacb78'
PUDONG = '#e4aa01'
OUTLINE = '#262626'
DOTTED = '#3c3c3c'
GRID = '#b4b4b4'
GRID_LABEL = '#555555'

NANJING_CODE = '320100'
PUDONG_CODE = '310115'


def fetch():
    DATA.mkdir(exist_ok=True)
    for fname, url in URLS.items():
        out = DATA / fname
        if not out.exists():
            print(f'Downloading {url}')
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as r:
                out.write_bytes(r.read())


def canvas_transform(bounds):
    """Affine matrix (shapely order) from EPSG:3857 metres to canvas mm.

    `bounds` (EPSG:3857) is scaled to fit inside the canvas minus MARGIN_MM and
    centred. The y axis points up, as in matplotlib.
    """
    x0, y0, x1, y1 = bounds
    k = min((W_MM - 2 * MARGIN_MM) / (x1 - x0),
            (H_MM - 2 * MARGIN_MM) / (y1 - y0))  # mm per metre
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    return [k, 0, 0, k, W_MM / 2 - k * cx, H_MM / 2 - k * cy]


def apply(gdf, mat):
    gdf = gdf.copy()
    gdf['geometry'] = gdf.geometry.apply(lambda g: affine_transform(g, mat))
    return gdf


def drop_specks(gdf, min_area=1e6):
    """Remove tiny islands (< min_area m^2 in EPSG:3857, i.e. ~1 km^2)."""
    gdf = gdf.explode(index_parts=False)
    gdf = gdf[gdf.area >= min_area]
    return gdf.dissolve(by='adcode', as_index=False, sort=False)


def merged(gdf, max_hole=0.6):
    """Union of all units, with tiny holes (gaps between adjacent polygons, mm^2) filled."""
    u = gdf.union_all()
    def fill(poly):
        return Polygon(poly.exterior, [r for r in poly.interiors if Polygon(r).area > max_hole])
    return unary_union([fill(p) for p in getattr(u, 'geoms', [u])])


def internal_borders(gdf, tol=0.08):
    """Each shared border exactly once (so the dot pattern is not doubled)."""
    outer = merged(gdf).boundary
    lines = unary_union(list(gdf.boundary))  # nodes and de-duplicates shared edges
    inner = lines.difference(outer.buffer(tol))
    return gpd.GeoSeries([linemerge(inner) if inner.geom_type == 'MultiLineString' else inner])


def load():
    fetch()
    jiangsu = gpd.read_file(DATA / 'datav_320000_full.json')
    shanghai = gpd.read_file(DATA / 'datav_310000_full.json')
    for g in (jiangsu, shanghai):
        g['adcode'] = g['adcode'].astype(str)
    ne = gpd.read_file(DATA / 'ne_50m_admin_1.zip', bbox=(105, 20, 130, 45))
    neighbours = ne[(ne.admin == 'China') & ~ne.name.isin(['Jiangsu', 'Shanghai'])]

    jiangsu = drop_specks(jiangsu.to_crs(3857))
    shanghai = drop_specks(shanghai.to_crs(3857))
    mat = canvas_transform(pd.concat([jiangsu, shanghai]).total_bounds)
    jiangsu = apply(jiangsu, mat)
    shanghai = apply(shanghai, mat)
    neighbours = apply(neighbours.to_crs(3857), mat)

    bg_fill, bg_lines = stitch_background(neighbours, merged(pd.concat([jiangsu, shanghai])))
    # only keep the background that is visible (slightly beyond, so lines run off the edge)
    clip = box(-1.5, -1.5, W_MM + 1.5, H_MM + 1.5)
    bg_fill = bg_fill.clip(clip)
    bg_lines = bg_lines.clip(clip)
    return jiangsu, shanghai, bg_fill, bg_lines, mat


def graticule(mat, step):
    """Meridians and parallels every `step` degrees that fall on the canvas.

    In Web Mercator these are straight vertical/horizontal lines, so each is
    returned as a single canvas coordinate: ([(lon, x_mm)], [(lat, y_mm)]).
    """
    from math import ceil, floor
    corners = gpd.GeoSeries(
        [Point(affine_inverse(mat, 0, 0)), Point(affine_inverse(mat, W_MM, H_MM))],
        crs=3857).to_crs(4326)
    (lon0, lat0), (lon1, lat1) = [(p.x, p.y) for p in corners]
    lons = [step * i for i in range(ceil(lon0 / step), floor(lon1 / step) + 1)]
    lats = [step * i for i in range(ceil(lat0 / step), floor(lat1 / step) + 1)]
    pts = gpd.GeoSeries([Point(lon, lats[0] if lats else lat0) for lon in lons]
                        + [Point(lons[0] if lons else lon0, lat) for lat in lats],
                        crs=4326).to_crs(3857)
    pts = [affine_transform(p, mat) for p in pts]
    return ([(lon, p.x) for lon, p in zip(lons, pts)],
            [(lat, p.y) for lat, p in zip(lats, pts[len(lons):])])


def affine_inverse(mat, x, y):
    """Canvas mm -> EPSG:3857 metres (inverse of the diagonal `canvas_transform`)."""
    a, _, _, e, xoff, yoff = mat
    return (x - xoff) / a, (y - yoff) / e


def stitch_background(neighbours, fg, tol=1.0):
    """Make the (Natural Earth) grey layer meet the (DataV) Jiangsu/Shanghai outline.

    The two datasets disagree slightly along the shared border, so grey border
    lines would peek out next to the black outline (or leave white gaps).
    Within `tol` mm of the outline: grey lines are removed, lines that end
    there are extended onto the outline, and the gap is filled with grey land.
    """
    near = fg.buffer(tol)
    land = neighbours.union_all()
    # morphological closing: fills slivers between the two layers without spilling into the sea
    closed = land.union(fg).buffer(tol).buffer(-tol)
    gap = near.intersection(closed).difference(fg)
    fill = gpd.GeoSeries([land.union(gap).difference(fg)])

    outline = fg.boundary
    lines = unary_union(list(neighbours.boundary)).difference(near)
    stitched = []
    for ls in getattr(lines, 'geoms', [lines]):
        coords = list(ls.coords)
        for end in (0, -1):
            p = Point(coords[end])
            if p.distance(outline) < tol * 1.1:
                q = nearest_points(outline, p)[0]
                if end == 0:
                    coords.insert(0, (q.x, q.y))
                else:
                    coords.append((q.x, q.y))
        stitched.append(LineString(coords))
    return fill, gpd.GeoSeries(stitched)


STYLE = dict(
    bg_lw=0.7,          # line widths in points
    js_lw=1.0,
    sh_lw=0.9,
    dot_lw=0.5,
    sh_dot_lw=0.4,
    dot_on=1.5,         # dash pattern in units of the dotted line width
    dot_off=3.0,
    grid_step=1,        # graticule spacing in degrees (None: no graticule)
    grid_label_step=2,  # only lines at multiples of this (degrees) are labelled
    grid_lw=0.15,
    grid_fontsize=5,
    grid_pad=0.6,       # label distance from the canvas edge and from its line (mm)
)


def draw(data, style=STYLE):
    jiangsu, shanghai, bg_fill, bg_lines, mat = data
    st = dict(STYLE, **style)
    fig = plt.figure(figsize=(W_MM / MM_PER_INCH, H_MM / MM_PER_INCH), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W_MM)
    ax.set_ylim(0, H_MM)
    ax.set_axis_off()

    # 1. grey background provinces
    bg_fill.plot(ax=ax, color=BG_LAND, linewidth=0)
    bg_lines.plot(ax=ax, color=BG_EDGE, linewidth=st['bg_lw'])

    # 1b. graticule: drawn over the background and hidden by Jiangsu/Shanghai
    # (same zorder as their fills, drawn earlier)
    if st['grid_step']:
        lons, lats = graticule(mat, st['grid_step'])
        for lon, x in lons:
            ax.axvline(x, color=GRID, linewidth=st['grid_lw'], zorder=1)
        for lat, y in lats:
            ax.axhline(y, color=GRID, linewidth=st['grid_lw'], zorder=1)
        fg = merged(pd.concat([jiangsu, shanghai]))
        label_graticule(ax, lons, lats, fg, st)

    dots = (0, (st['dot_on'], st['dot_off']))
    # 2. Jiangsu: prefectures with dotted internal borders, solid outline
    jiangsu.plot(ax=ax, color=JS_FILL, linewidth=0)
    internal_borders(jiangsu).plot(ax=ax, color=DOTTED, linewidth=st['dot_lw'], linestyle=dots)
    gpd.GeoSeries([merged(jiangsu)]).boundary.plot(
        ax=ax, color=OUTLINE, linewidth=st['js_lw'], joinstyle='round', capstyle='round')
    # Nanjing is painted last, on top of (and half covering) the outline
    jiangsu[jiangsu.adcode == NANJING_CODE].plot(ax=ax, color=NANJING, linewidth=0)

    # 3. Shanghai: districts, Pudong highlighted
    sh_colors = [PUDONG if c == PUDONG_CODE else SH_FILL for c in shanghai.adcode]
    shanghai.plot(ax=ax, color=sh_colors, linewidth=0)
    internal_borders(shanghai).plot(ax=ax, color=DOTTED, linewidth=st['sh_dot_lw'], linestyle=dots)
    gpd.GeoSeries([merged(shanghai)]).boundary.plot(
        ax=ax, color=OUTLINE, linewidth=st['sh_lw'], joinstyle='round', capstyle='round')
    return fig


def label_graticule(ax, lons, lats, fg, st):
    """Label meridians along the bottom edge and parallels along the right edge
    (only those at multiples of `grid_label_step` degrees).

    Each label is placed beside its line (longitudes: right of it, else left;
    latitudes: above it, else below). A label that would run off the canvas,
    overlap Jiangsu/Shanghai (`fg`, canvas mm) or overlap an earlier label is
    tried at the next position and dropped if none fits.
    """
    pad, fs = st['grid_pad'], st['grid_fontsize']
    renderer = ax.figure.canvas.get_renderer()
    inv = ax.transData.inverted()
    canvas = box(0, 0, W_MM, H_MM)
    blocked = fg.buffer(0.3)

    def place(text, candidates):
        nonlocal blocked
        for x, y, ha, va in candidates:
            t = ax.text(x, y, text, fontsize=fs, color=GRID_LABEL, ha=ha, va=va, zorder=10)
            (x0, y0), (x1, y1) = inv.transform(t.get_window_extent(renderer))
            extent = box(x0, y0, x1, y1)
            if canvas.contains(extent) and not extent.intersects(blocked):
                blocked = blocked.union(extent.buffer(pad))
                return
            t.remove()

    def labelled(values):
        step = st['grid_label_step']
        return [(v, c) for v, c in values if abs(v / step - round(v / step)) < 1e-9]

    for lon, x in labelled(lons):
        place(f'{lon:g}°E', [(x + pad, pad, 'left', 'bottom'),
                             (x - pad, pad, 'right', 'bottom')])
    for lat, y in labelled(lats):
        place(f'{lat:g}°N', [(W_MM - pad, y + pad, 'right', 'bottom'),
                             (W_MM - pad, y - pad, 'right', 'top')])


def main():
    fig = draw(load())
    for ext in ('png', 'svg', 'pdf'):
        fig.savefig(HERE / f'fig1a_map.{ext}', dpi=DPI, facecolor='white')
    print('Wrote', HERE / 'fig1a_map.{png,svg,pdf}')


if __name__ == '__main__':
    main()
