import os
import ee
import geemap
import warnings
import rasterio
import numpy as np
import geopandas as gpd
import rioxarray as rxr
from rasterio.mask import mask
from rasterio.merge import merge
from shapely.geometry import box
from shapely.geometry import shape
from rasterio.features import shapes
warnings.filterwarnings("ignore", category=RuntimeWarning)

class WaterMaskClass:

    def __init__(self):
        pass

    def get_jrc_gcer(self, bbox_geom, percentage_occurrence, output_path, database_path):

        '''
        Generates the water mask using the JRC Global Surface Water dataset from GCER database
        Input:
            bbox_geom (GeoDataFrame): input shapefile
            percentage_occurrence (int): percentage of water occurrence
            output_path (str): path to save the output shapefile
        Return:
            str: path to the water mask shapefile
        '''

        #image_path = r'Z:\guser\tml\mypapers\hls_synthetic\JRC_DATA'
        image_path =  database_path

        files = [f for f in os.listdir(image_path) if os.path.isfile(os.path.join(image_path, f)) and f.endswith('.tif')]

        count = 0  # in case of more than one image intersect the lake
        images_list = []

        for file in files:

            with rasterio.open(image_path + '/' + file) as src:

                raster_bounds = src.bounds
                geom = bbox_geom.geometry
                bounds_aux = bbox_geom.bounds

                if (bounds_aux[['minx']].values > raster_bounds[2] or bounds_aux[['maxx']].values < raster_bounds[0] or
                        bounds_aux[['miny']].values > raster_bounds[3] or bounds_aux[['maxx']].values < raster_bounds[
                            1]):
                    continue

                else:
                    try:
                        subset, out_transform = mask(src, geom, nodata=src.nodata, crop=True, indexes=1)
                    except:
                        continue

                    subset_clip = subset > percentage_occurrence

                    out_meta = src.meta.copy()
                    out_meta.update({"driver": "GTiff",
                                     "height": subset_clip.shape[0],
                                     "width": subset_clip.shape[1],
                                     "transform": out_transform,
                                     "nodata": 0})

                    with rasterio.open(output_path + '/' + 'water_mask_jrc_' + str(count) + '.tif', "w",
                                       **out_meta) as dest:
                        dest.write(subset_clip.astype(int), 1)
                        dest.nodata = 0

                    images_list.append(output_path + '/' + 'water_mask_jrc_' + str(count) + '.tif')
                    count += 1

        if count > 1:
            src_files_to_mosaic = []

            for fp in images_list:
                src = rasterio.open(fp)
                src_files_to_mosaic.append(src)

            # Merge the tiffs into one mosaic
            mosaic, out_trans = merge(src_files_to_mosaic)

            out_meta = src_files_to_mosaic[0].meta.copy()

            out_meta.update({
                "driver": "GTiff",
                "height": mosaic.shape[1],
                "width": mosaic.shape[2],
                "transform": out_trans,
            })

            with rasterio.open(output_path + '/' + 'water_mask_jrc.tif', "w", **out_meta) as dest:
                dest.write(mosaic)

            return output_path + '/' + 'water_mask_jrc.tif'

        else:
            return output_path + '/' + 'water_mask_jrc_0.tif'

    def raster2shp(self, raster_numpy, raster_transform, raster_crs):

        '''
        Converts a raster to a shapefile
        Input:
            raster_numpy (numpy array): raster array
            raster_transform (affine.Affine): raster transform
            raster_crs (str): raster crs
        Return:
            GeoDataFrame: shapefile
        '''

        results = (
            {'properties': {'raster_val': v}, 'geometry': s}
            for s, v in shapes(raster_numpy.astype(np.int16), transform=raster_transform)
            if v  # Only take shapes with raster_val = True (i.e., v=1)
        )

        geometries = [shape(feature['geometry']) for feature in results]

        gdf = gpd.GeoDataFrame(geometry=geometries, crs=raster_crs)

        return gdf

    def get_wm_jrc(self, image_path, percentage_occurrence, output_path):

        '''
        Generates the water mask using the JRC Global Surface Water dataset
        Input:
            gdf_input (GeoDataFrame): input shapefile
            percentage_occurrence (int): percentage of water occurrence
            output_path (str): path to save the output shapefile
            output_crs (int): output crs
        Return:
            str: path to the water mask shapefile
        '''

        ee.Initialize()

        # Return the bouding box of the input image
        xda_nir = rxr.open_rasterio(image_path)
        bounds = xda_nir.rio.bounds()
        xmin, ymin, xmax, ymax = bounds
        bounding_box = box(xmin, ymin, xmax, ymax)
        gdf_input = gpd.GeoDataFrame({"geometry": [bounding_box]}, crs=xda_nir.rio.crs).to_crs(epsg=4326)

        output_crs = xda_nir.rio.crs.to_epsg()

        feature_collection = geemap.geopandas_to_ee(gdf_input)

        gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
        occurrence = gsw.select('occurrence')
        water_mask = occurrence.gt(percentage_occurrence)

        water_mask_bbox = water_mask.clip(feature_collection).reproject('EPSG:' + str(output_crs), None, 30)

        # Save the JRC water mask as .tiff file
        geemap.ee_export_image(water_mask_bbox, filename=output_path + '/' + 'water_mask_jrc.tif',
                               region=feature_collection.geometry().bounds())

        try:
            with rasterio.open(output_path + '/' + 'water_mask_jrc.tif') as src:
                wm_jrc = src.read(1)
        except:
            jrc_img_path = self.get_jrc_gcer(gdf_input, percentage_occurrence, output_path)

            with rasterio.open(jrc_img_path) as src:
                wm_jrc = src.read(1)

        gdf = self.raster2shp(wm_jrc, src.transform, src.crs)

        # Clean the water_mask_jrc files

        return gdf

    def get_mask_algal(self, mask_method, output_path, img_path = None):

        '''
        Classify algal blooms based on the NDCI index
        Input:
            img_path (str): path to the image folder
            ndci_threshold (float): NDCI threshold
        Return:
            xarray: algal bloom mask
        '''

        if mask_method == "image_based":

            band_path = [i for i in os.listdir(img_path) if i.endswith('.tif')]

            green_band = [band for band in band_path if 'B03' in band]
            nir_band = [band for band in band_path if 'B05' in band]
            swir_band = [band for band in band_path if 'B06' in band or 'B11' in band or 'B12' in band]

            xda_green = rxr.open_rasterio(img_path + '/' + green_band[0])/10000
            xda_nir = rxr.open_rasterio(img_path + '/' + nir_band[0])/10000
            xda_swir = rxr.open_rasterio(img_path + '/' + swir_band[0])/10000

            # Apply glint correction
            green_band_glint = np.where((xda_green - xda_swir) < 0, xda_green, (xda_green - xda_swir))

            # 2. Define mask for SWIR band to exclude land pixels
            swir_mask = xda_swir < 0.03

            # Calculate MNDWI index to mask land
            mndwi = (green_band_glint - xda_swir) / (green_band_glint + xda_swir)
            mndwi_mask = mndwi > 0.3

            water_mask = swir_mask & mndwi_mask

            water_shp = self.raster2shp(water_mask, xda_nir.rio.transform, xda_nir.rio.crs)

        else:
            # 3. Apply JRC water mask
            water_shp = self.get_wm_jrc(img_path, 75, img_path)

        water_shp.to_file(output_path)