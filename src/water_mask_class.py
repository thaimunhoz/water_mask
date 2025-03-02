import os
import ee
import geemap
import logging
import rasterio
import numpy as np
import geopandas as gpd
import rioxarray as rxr
from typing import Optional
from rasterio.mask import mask
from rasterio.merge import merge
from shapely.geometry import box, shape
from rasterio.features import shapes
from rasterio.transform import Affine
from rasterio.crs import CRS

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class WaterMaskClass:

    """
    A class to extract water masks from satellite imagery using either an image-based approach
    or the JRC Global Surface Water dataset.
    """

    def __init__(self):
        pass

    def get_jrc_gcer(
        self,
        bbox_geom: gpd.GeoDataFrame,
        percentage_occurrence: int,
        output_path: str,
        database_path: str,
    ) -> str:
        """
        Generate a water mask using the JRC Global Surface Water dataset from the GCER database.

        Args:
            bbox_geom (gpd.GeoDataFrame): Input bounding box geometry.
            percentage_occurrence (int): Percentage of water occurrence threshold.
            output_path (str): Path to save the output shapefile.
            database_path (str): Path to the JRC dataset directory.

        Returns:
            str: Path to the generated water mask shapefile.
        """
        if not os.path.exists(database_path):
            raise FileNotFoundError(f"Database path {database_path} does not exist.")

        files = [
            f for f in os.listdir(database_path)
            if os.path.isfile(os.path.join(database_path, f)) and f.endswith(".tif")
        ]

        if not files:
            raise ValueError(f"No TIFF files found in {database_path}.")

        images_list = []
        for count, file in enumerate(files):
            file_path = os.path.join(database_path, file)
            try:
                with rasterio.open(file_path) as src:
                    raster_bounds = src.bounds
                    bounds_aux = bbox_geom.bounds

                    if (
                        bounds_aux["minx"].values[0] > raster_bounds[2]
                        or bounds_aux["maxx"].values[0] < raster_bounds[0]
                        or bounds_aux["miny"].values[0] > raster_bounds[3]
                        or bounds_aux["maxy"].values[0] < raster_bounds[1]
                    ):
                        continue

                    subset, out_transform = mask(src, bbox_geom.geometry, crop=True, indexes=1)
                    subset_clip = subset > percentage_occurrence

                    out_meta = src.meta.copy()
                    out_meta.update(
                        {
                            "driver": "GTiff",
                            "height": subset_clip.shape[0],
                            "width": subset_clip.shape[1],
                            "transform": out_transform,
                            "nodata": 0,
                        }
                    )

                    output_file = os.path.join(output_path, f"water_mask_jrc_{count}.tif")
                    with rasterio.open(output_file, "w", **out_meta) as dest:
                        dest.write(subset_clip.astype(int), 1)

                    images_list.append(output_file)
            except Exception as e:
                logger.warning(f"Failed to process {file_path}: {e}")
                continue

        if not images_list:
            raise RuntimeError("No valid images were processed.")

        if len(images_list) > 1:
            mosaic, out_trans = merge([rasterio.open(fp) for fp in images_list])
            out_meta = rasterio.open(images_list[0]).meta.copy()
            out_meta.update(
                {
                    "driver": "GTiff",
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": out_trans,
                }
            )

            mosaic_path = os.path.join(output_path, "water_mask_jrc.tif")
            with rasterio.open(mosaic_path, "w", **out_meta) as dest:
                dest.write(mosaic)

            return mosaic_path
        else:
            return images_list[0]

    def raster2shp(
        self, raster_numpy: np.ndarray, raster_transform: Affine, raster_crs: CRS
    ) -> gpd.GeoDataFrame:
        """
        Convert a raster array to a shapefile.

        Args:
            raster_numpy (np.ndarray): Raster array.
            raster_transform (Affine): Raster transform.
            raster_crs (CRS): Raster coordinate reference system.

        Returns:
            gpd.GeoDataFrame: GeoDataFrame containing the vectorized raster.
        """
        results = list(
            {"properties": {"raster_val": v}, "geometry": s}
            for s, v in shapes(np.asarray(raster_numpy, dtype=np.int16), transform=raster_transform)
            if v  # Only take shapes with raster_val = True (i.e., v=1)
        )
        geometries = [shape(feature["geometry"]) for feature in results]

        return gpd.GeoDataFrame(geometry=geometries, crs=raster_crs)

    def get_wm_jrc(
        self, image_path: str, percentage_occurrence: int
    ) -> gpd.GeoDataFrame:
        """
        Generate a water mask using the JRC Global Surface Water dataset.

        Args:
            image_path (str): Path to the input image.
            percentage_occurrence (int): Percentage of water occurrence threshold.

        Returns:
            gpd.GeoDataFrame: GeoDataFrame containing the water mask.
        """
        if not os.path.exists(image_path):
            raise FileNotFoundError(f"Image path {image_path} does not exist.")

        xda_nir = rxr.open_rasterio(image_path)
        bounds = xda_nir.rio.bounds()
        bounding_box = box(*bounds)
        gdf_input = gpd.GeoDataFrame({"geometry": [bounding_box]}, crs=xda_nir.rio.crs).to_crs(epsg=4326)

        feature_collection = geemap.geopandas_to_ee(gdf_input)
        gsw = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
        occurrence = gsw.select("occurrence")
        water_mask = occurrence.gt(percentage_occurrence)

        output_crs = xda_nir.rio.crs.to_epsg()
        water_mask_bbox = water_mask.clip(feature_collection).reproject(f"EPSG:{output_crs}", None, 30)

        output_tif = "water_mask_jrc.tif"

        geemap.ee_export_image(
            water_mask_bbox,
            filename=output_tif,
            region=feature_collection.geometry().bounds(),
        )

        try:
            with rasterio.open(output_tif) as src:
                wm_jrc = src.read(1)
                gdf = self.raster2shp(wm_jrc, src.transform, src.crs)
        except Exception as e:
            logger.warning(f"Failed to process JRC mask: {e}")

        # Remove the temporary file
        os.remove(output_tif)

        return gdf

    def get_mask_water(
        self, mask_method: str, img_path: str
    ) -> None:
        """
        Generate a water mask using either an image-based approach or the JRC dataset.

        Args:
            mask_method (str): Method to use ("image_based" or "jrc_based").
            img_path (Optional[str]): Path to the input image folder (required for "image_based" method).

        Raises:
            ValueError: If `img_path` is not provided for the "image_based" method.
        """

        if mask_method == "image_based":
            if not img_path:
                raise ValueError("img_path is required for the image_based method.")

            band_path = [i for i in os.listdir(img_path) if i.endswith(".tif") or i.endswith(".TIF")]
            green_band = next((band for band in band_path if "B03" in band or "B3" in band), None)
            nir_band = next((band for band in band_path if "B05" in band or "B5" in band), None)
            swir_band = next((band for band in band_path if "B06" in band or "B6" in band or "B11" in band or "B12" in band), None)

            if not all([green_band, nir_band, swir_band]):
                raise ValueError("Required bands (B03, B05, B06/B11/B12) not found in the image folder.")

            xda_green = rxr.open_rasterio(os.path.join(img_path, green_band)) #/ 10000
            xda_swir = rxr.open_rasterio(os.path.join(img_path, swir_band)) #/ 10000

            green_band_glint = np.where((xda_green - xda_swir) < 0, xda_green, (xda_green - xda_swir))
            swir_mask = xda_swir < 0.03
            mndwi = (green_band_glint - xda_swir) / (green_band_glint + xda_swir)
            mndwi_mask = mndwi > 0.3
            water_mask = swir_mask & mndwi_mask

            with rasterio.open(os.path.join(img_path, green_band)) as src:
                _transform = src.transform
                _crs = src.crs

            water_shp = self.raster2shp(water_mask, _transform, _crs)
        else:
            water_shp = self.get_wm_jrc(img_path, 75)

        water_shp.to_file("water_mask.shp")

        logger.info(f"Water mask saved!")