#!/usr/bin/env python3
import argparse, os, sys, textwrap
from osgeo import gdal, osr

def detect_pair(path_like):
    """Return (hh_path, hv_path) from a directory or a file list string."""
    if os.path.isdir(path_like):
        files = [os.path.join(path_like, f) for f in os.listdir(path_like)]
    else:
        # path_like may be a single file or wildcard expanded by shell
        files = [p for p in path_like.split(",") if p.strip()]
        if len(files) == 1 and os.path.isdir(files[0]):
            files = [os.path.join(files[0], f) for f in os.listdir(files[0])]
    hh = [f for f in files if f.endswith("HHHH.tif")]
    hv = [f for f in files if f.endswith("HVHV.tif")]
    if len(hh) != 1 or len(hv) != 1:
        raise SystemExit("Could not uniquely find HHHH.tif and HVHV.tif in: " + path_like)
    return hh[0], hv[0]

def read_geo(meta_path):
    ds = gdal.Open(meta_path, gdal.GA_ReadOnly)
    if ds is None:
        raise SystemExit(f"Failed to open {meta_path}")
    gt = ds.GetGeoTransform()
    w, h = ds.RasterXSize, ds.RasterYSize
    proj_wkt = ds.GetProjectionRef()
    return ds, (w, h), gt, proj_wkt

def same_grid(info_a, info_b, tol=1e-9):
    (wA,hA), gtA, projA = info_a[1], info_a[2], info_a[3]
    (wB,hB), gtB, projB = info_b[1], info_b[2], info_b[3]
    if (wA,hA)!=(wB,hB): return False
    if projA.strip()!=projB.strip(): return False
    return all(abs(a-b) <= tol for a,b in zip(gtA, gtB))

VRT_TMPL = """<VRTDataset rasterXSize="{w}" rasterYSize="{h}">
  <SRS>{srs}</SRS>
  <GeoTransform>{gt0}, {gt1}, {gt2}, {gt3}, {gt4}, {gt5}</GeoTransform>

  <VRTRasterBand dataType="Float32" band="1" subClass="VRTSourcedRasterBand">
    <Description>HHHH_power</Description>
    <SimpleSource>
      <SourceFilename relativeToVRT="0">{hh}</SourceFilename>
      <SourceBand>1</SourceBand>
    </SimpleSource>
  </VRTRasterBand>

  <VRTRasterBand dataType="Float32" band="2" subClass="VRTSourcedRasterBand">
    <Description>HVHV_power</Description>
    <SimpleSource>
      <SourceFilename relativeToVRT="0">{hv}</SourceFilename>
      <SourceBand>1</SourceBand>
    </SimpleSource>
  </VRTRasterBand>

  <VRTRasterBand dataType="Float32" band="3" subClass="VRTDerivedRasterBand">
    <Description>HH_dB</Description>
    <PixelFunctionLanguage>Python</PixelFunctionLanguage>
    <PixelFunctionType>HH_to_dB</PixelFunctionType>
    <PixelFunctionCode><![CDATA[
import numpy as np
_eps = 1e-32
def HH_to_dB(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize, raster_ysize, buf_radius, gt, **kwargs):
    hh = in_ar[0].astype(np.float32)
    out_ar[:] = 10.0 * np.log10(np.maximum(hh, _eps))
]]></PixelFunctionCode>
    <ComplexSource><SourceBand>1</SourceBand></ComplexSource>
  </VRTRasterBand>

  <VRTRasterBand dataType="Float32" band="4" subClass="VRTDerivedRasterBand">
    <Description>HV_dB</Description>
    <PixelFunctionLanguage>Python</PixelFunctionLanguage>
    <PixelFunctionType>HV_to_dB</PixelFunctionType>
    <PixelFunctionCode><![CDATA[
import numpy as np
_eps = 1e-32
def HV_to_dB(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize, raster_ysize, buf_radius, gt, **kwargs):
    hv = in_ar[0].astype(np.float32)
    out_ar[:] = 10.0 * np.log10(np.maximum(hv, _eps))
]]></PixelFunctionCode>
    <ComplexSource><SourceBand>2</SourceBand></ComplexSource>
  </VRTRasterBand>

  <VRTRasterBand dataType="Float32" band="5" subClass="VRTDerivedRasterBand">
    <Description>HH_minus_HV_dB</Description>
    <PixelFunctionLanguage>Python</PixelFunctionLanguage>
    <PixelFunctionType>HH_minus_HV</PixelFunctionType>
    <PixelFunctionCode><![CDATA[
import numpy as np
_eps = 1e-32
def HH_minus_HV(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize, raster_ysize, buf_radius, gt, **kwargs):
    hh = np.maximum(in_ar[0].astype(np.float32), _eps)
    hv = np.maximum(in_ar[1].astype(np.float32), _eps)
    out_ar[:] = 10.0*np.log10(hh) - 10.0*np.log10(hv)
]]></PixelFunctionCode>
    <ComplexSource><SourceBand>1</SourceBand></ComplexSource>
    <ComplexSource><SourceBand>2</SourceBand></ComplexSource>
  </VRTRasterBand>

  <VRTRasterBand dataType="Float32" band="6" subClass="VRTDerivedRasterBand">
    <Description>HH_over_HV_lin</Description>
    <PixelFunctionLanguage>Python</PixelFunctionLanguage>
    <PixelFunctionType>HH_over_HV</PixelFunctionType>
    <PixelFunctionCode><![CDATA[
import numpy as np
_eps = 1e-32
def HH_over_HV(in_ar, out_ar, xoff, yoff, xsize, ysize, raster_xsize, raster_ysize, buf_radius, gt, **kwargs):
    hh = in_ar[0].astype(np.float32)
    hv = in_ar[1].astype(np.float32)
    out_ar[:] = np.clip(hh / np.maximum(hv, _eps), 0.0, 1e6)
]]></PixelFunctionCode>
    <ComplexSource><SourceBand>1</SourceBand></ComplexSource>
    <ComplexSource><SourceBand>2</SourceBand></ComplexSource>
  </VRTRasterBand>

</VRTDataset>
"""

def main():
    ap = argparse.ArgumentParser(
        description="Create a polarimetric combo VRT (HH, HV, HH_dB, HV_dB, HH–HV dB, HH/HV).",
        epilog=textwrap.dedent("""\
            Examples:
              python make_pol_combo_vrt.py --in "/path/to/folder" --out nisar_pol_combo.vrt
              python make_pol_combo_vrt.py --hh /path/..HHHH.tif --hv /path/..HVHV.tif --out nisar_pol_combo.vrt
        """)
    )
    ap.add_argument("--in", dest="in_path", help="Directory or comma-separated file list containing exactly one HHHH.tif and one HVHV.tif")
    ap.add_argument("--hh", help="Explicit path to HHHH.tif")
    ap.add_argument("--hv", help="Explicit path to HVHV.tif")
    ap.add_argument("--out", required=True, help="Output VRT path")
    args = ap.parse_args()

    if (args.hh and args.hv) is None and args.in_path is None:
        ap.error("Provide either --in DIR_OR_LIST or both --hh and --hv.")

    if args.in_path:
        hh_path, hv_path = detect_pair(args.in_path)
    else:
        if not (args.hh and args.hv):
            ap.error("Provide both --hh and --hv when not using --in.")
        hh_path, hv_path = args.hh, args.hv

    # Read metadata
    ds_hh, size_hh, gt_hh, wkt_hh = read_geo(hh_path)
    ds_hv, size_hv, gt_hv, wkt_hv = read_geo(hv_path)

    if not same_grid((ds_hh, size_hh, gt_hh, wkt_hh), (ds_hv, size_hv, gt_hv, wkt_hv)):
        raise SystemExit("ERROR: HHHH and HVHV rasters differ in size, CRS, or geotransform.")

    w, h = size_hh
    gt = gt_hh

    # Clean WKT (optional)
    srs = osr.SpatialReference()
    srs.ImportFromWkt(wkt_hh)
    srs_wkt_pretty = srs.ExportToWkt()

    xml = VRT_TMPL.format(
        w=w, h=h,
        srs=srs_wkt_pretty,
        gt0=f"{gt[0]:.15f}", gt1=f"{gt[1]:.15f}", gt2=f"{gt[2]:.15f}",
        gt3=f"{gt[3]:.15f}", gt4=f"{gt[4]:.15f}", gt5=f"{gt[5]:.15f}",
        hh=os.path.abspath(hh_path),
        hv=os.path.abspath(hv_path),
    )

    with open(args.out, "w", encoding="utf-8") as f:
        f.write(xml)

    print(f"Wrote VRT: {args.out}")
    print("Bands:")
    print("  1: HHHH_power")
    print("  2: HVHV_power")
    print("  3: HH_dB")
    print("  4: HV_dB")
    print("  5: HH_minus_HV_dB")
    print("  6: HH_over_HV_lin")

if __name__ == "__main__":
    main()
