import os
import pyproj
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from osgeo import gdal, osr
from WBT.whitebox_tools import WhiteboxTools
from sys import getsizeof

path_to_file = (os.getcwd()).replace("\\", "/")
wbt = WhiteboxTools()
name_reservoir = input("Введите название водохранилища: ")
number_of_isolines = int(input("Введите число изолиний: "))
step = int(input("Введите шаг глубины: "))

if os.path.exists(f"{path_to_file}/{name_reservoir}"):
    pass
else:
    os.mkdir(f"{path_to_file}/{name_reservoir}")

if os.path.exists(f"{path_to_file}/DRM"):
    pass
else:
    os.mkdir(f"{path_to_file}/DRM")

if os.path.exists(f"{path_to_file}/DEMs"):
    pass
else:
    os.mkdir(f"{path_to_file}/DEMs")

if os.path.exists(f"{path_to_file}/DEMs"):
    pass
else:
    os.mkdir(f"{path_to_file}/DEMs")

if os.path.exists(f"{path_to_file}/Исходная БД"):
    pass
else:
    os.mkdir(f"{path_to_file}/Исходная БД")


tension = -2
points = 2

def filter_ndarray(three_d_array, points):
    shape = three_d_array.shape
    
    # Flatten the 3D array into a 2D plane
    two_d_array = three_d_array.reshape(shape[0], -1)

    if points == 0:
        # Create an array with z indices
        z_indices_array = np.empty_like(three_d_array)

        for z_index, _ in np.ndenumerate(three_d_array):
            z_indices_array[z_index] = z_index[0]

        return z_indices_array
    
    elif points == 1:
        # Z-axis indices with minimum value for each 2D plane cell
        min_indices = np.argmin(two_d_array, axis=0)
        
        # Expand indices back to 2D space
        min_indices_2d = min_indices.reshape((shape[1], shape[2]))
        
        # Array of minimum values
        min_values_array = np.min(three_d_array, axis=0)
        
        return min_indices_2d, min_values_array
    
    else:
        # N indices of z-axis with minimum value for each cell of 2D plane
        sorted_indices = np.argsort(two_d_array, axis=0)[:points]
        sorted_indices = sorted_indices.astype(np.int8)
        
        # Expand indices back to 2D space
        sorted_indices_2d = sorted_indices.reshape((points, shape[1], shape[2]))
        
        # Arrays with minimum values
        sorted_values_array = np.sort(three_d_array, axis=0)[:points]
        
        return sorted_indices_2d, sorted_values_array

# Function of ICW interpolation
def icw_interpolation(value_array, distance_array, power):
    collapsed_array = np.zeros_like(value_array)
    numerator = 0
    denominator = 0
    for z in range(value_array.shape[0]):
        temp_z = value_array[z]
        dist_z = distance_array[z]

        numerator += np.power(dist_z, power) * temp_z
        denominator += np.power(dist_z, power)

    
    collapsed_array = numerator / denominator
    
    numerator, denominator = None, None

    return collapsed_array

def my_callback(value):
    if not "%" in value:
        print(value)

image_ds = gdal.Open(f"{path_to_file}/DEMs/topo_70.tif")
my_table = pd.read_excel(f"{path_to_file}/Исходная БД/7_IFA_rusgidro_main_DB_3D_ZEYSKOE_7feb.xlsx", skiprows = 2)
my_table = my_table.assign(ht_real = my_table['wat_level'] - my_table['horiz_sampl'])

# Get the geotransform and projection from the original raster
projection = image_ds.GetProjection()

# Get the EPSG code from the projection
srs = osr.SpatialReference(wkt=projection)
epsg = srs.GetAttrValue('AUTHORITY', 1)

# Define the WGS-84 (EPSG:4326) and the target CRS
wgs84 = pyproj.CRS("EPSG:4326")
target_crs = pyproj.CRS(f"EPSG:{epsg}")

# Create a transformer to convert from WGS-84 to the target CRS
transformer = pyproj.Transformer.from_crs(wgs84, target_crs, always_xy=True)

# Reproject the coordinates in the DataFrame
my_table[['lon', 'lat']] = my_table.apply(lambda row: pd.Series(transformer.transform(row['lon'], row['lat'])), axis=1)

band = image_ds.GetRasterBand(1)
dem = band.ReadAsArray()

no_data = band.GetNoDataValue()
geotransform = image_ds.GetGeoTransform()
min_x = geotransform[0]
min_y = geotransform[3] + geotransform[5] * image_ds.RasterYSize
max_x = geotransform[0] + geotransform[1] * image_ds.RasterXSize
max_y = geotransform[3]

rows, cols = dem.shape
driver = gdal.GetDriverByName("GTiff")
print(rows, cols)

# Get the pixel width and height from the geotransform
pixel_width = abs(geotransform[1])  # Absolute value to handle negative scaling
pixel_height = abs(geotransform[5])

# Calculate the area of one pixel
pixel_area = pixel_width * pixel_height

level_max = np.nanmax(dem)

# calculating Area in square metres
dem1 = dem[dem != no_data]
dem_new = dem1[dem1 < level_max]
count = len(dem_new)
S = count*pixel_area

niz_biefs = ["Нижний бьеф", "НБ", "нижний бьеф", "нижний бьеф плотины"]
all_res = pd.DataFrame()
dict_max_res = {}


unique_expeditions = list(my_table["expedition"].unique())
for one_exp in unique_expeditions:
    filtered_df = my_table[my_table["expedition"] == one_exp]
    filtered_df = filtered_df[~filtered_df["region"].isin(niz_biefs)]
    unique_stations = list(filtered_df["station_num"].unique())

    w = 0.5*np.sqrt(S/len(unique_stations))



    unique_stations_df = filtered_df.drop_duplicates(subset=["station_num"]).copy()
    min_dist_arr = np.array([])
    for station in unique_stations:
        distance_arr = np.array([])
        # Get the coordinates of the station dynamically
        station_coords = unique_stations_df.loc[filtered_df["station_num"] == station, ['lon', 'lat']].values
        station_coords = station_coords[0]
        point_x, point_y = station_coords[0], station_coords[1]
        for station in unique_stations:
            station_coords_2 = filtered_df.loc[filtered_df["station_num"] == station, ['lon', 'lat']].values
            station_coords_2 = station_coords_2[0]
            point_x_2, point_y_2 = station_coords_2[0], station_coords_2[1]
            
            # Calculate the distances between the station point and the current row using Euclidean distance
            distance = np.sqrt((point_x - point_x_2) ** 2 + (point_y - point_y_2) ** 2)
            distance_arr = np.append(distance_arr, distance)
        distance_arr = distance_arr[distance_arr != 0]
        min_dist = np.min(distance_arr)
        min_dist_arr = np.append(min_dist_arr, min_dist)
    min_dist_arr = np.unique(min_dist_arr)
    min_mean_kot = np.mean(min_dist_arr)/2
    threshold_distance = np.percentile(min_dist_arr, 5)
    
    print("w:", w)
    print("w_percentile:", threshold_distance)
    print("kotel:", min_mean_kot)
    key = one_exp
    dict_max_res[key]= pixel_width/(pixel_width*10/threshold_distance)

    the_dict_res = {"reservoir": name_reservoir, "expedition": one_exp, "res_kotel": min_mean_kot, "res_unif_distr": w, "res_kotel_finest": threshold_distance, "area": S/1000000}
    res_dat = pd.DataFrame(the_dict_res, index=[0])
    all_res = pd.concat([all_res, res_dat], ignore_index=True)
    
all_res.to_excel(f"{path_to_file}/table_res.xlsx")

new_res = pixel_width / (pixel_width * 10 / threshold_distance)
int_new_res = int(new_res)
# Specify the resampling algorithm and resample the dataset
resample_alg = gdal.GRA_Bilinear
resampled_ds = gdal.Warp(f"{path_to_file}/topo_{int_new_res}.tif", f"{path_to_file}/DEMs/topo_70.tif", xRes=new_res, yRes=new_res, resampleAlg=resample_alg)

image_ds_res = gdal.Open(f"{path_to_file}/topo_{int_new_res}.tif")
band_res = image_ds_res.GetRasterBand(1)
dem_res = band_res.ReadAsArray()
rows2, cols2 = dem_res.shape
print(rows2, cols2)
projection_res = image_ds_res.GetProjection()
no_data = band_res.GetNoDataValue()
geotransform_res = image_ds_res.GetGeoTransform()
print(geotransform_res)

all_means = pd.DataFrame()
# removing points out of reservoir

for one_exp in unique_expeditions[0:1]:

    filtered_df = my_table[my_table["expedition"] == one_exp]
    filtered_df = filtered_df[~filtered_df["region"].isin(niz_biefs)]
    unique_stations = list(filtered_df["station_num"].unique())

    max_elev = np.nanmax(filtered_df["wat_level"])

    all_stations_df = pd.DataFrame()
    print(unique_stations)
    for stat in unique_stations:
        print(stat)
        one_station = filtered_df[filtered_df["station_num"] == stat]
        df = one_station[["expedition", "station_num", "lat", "lon", "wat_level", "wat_depth", "horiz_sampl", "ht_real", "tmp", "oxy_mg", "ch4_wat", "bottom_flag", "surface_flag"]]
        new_df = df.sort_values(by = "ht_real", ascending=False)
        new_df = new_df.drop_duplicates(subset=['ht_real'], keep='last')
        
        # Set 'Height' column as the index
        new_df.set_index('ht_real', inplace=True)

        # Linear interpolation nodata values based on 'Height' index
        new_df['tmp'] = new_df['tmp'].interpolate(method='index', limit_direction='both')
        new_df['ch4_wat'] = new_df['ch4_wat'].interpolate(method='index', limit_direction='both')
        new_df['oxy_mg'] = new_df['oxy_mg'].interpolate(method='index', limit_direction='both')

        new_df.reset_index(inplace=True)
        
        # Linear interpolation of hydrological vertical

        if new_df["bottom_flag"].isin([1]).any():
            min_elev = new_df["ht_real"].min()
        else: 
            #min_elev = max_elev - new_df["wat_depth"].mean()
            min_elev = new_df["ht_real"].min()
        try:
            n = int((max_elev-min_elev)/step)+1
        except ValueError:
            print("Err!", min_elev, max_elev)
            min_elev = max_elev
            n = 1
        else:
            print(min_elev)
            print(max_elev)
        
        expedition = [one_exp]*n
        stations = [stat]*n
        lats = [new_df["lat"].mean()]*n
        lons = [new_df["lon"].mean()]*n
        lt = []
        lt_2 = []
        for i in range(n):
            elev = max_elev - i*step
            lt.append(elev)
            lt_2.append(i*step)
        sampls_heights = new_df["ht_real"].tolist()
        print(lt)
        sampls_tmps = new_df["tmp"].tolist()
        sampls_oxys = new_df["oxy_mg"].tolist()
        sampls_met = new_df["ch4_wat"].tolist()
        if len(sampls_heights) == 1:
            ht = sampls_heights[0]
            tmp = sampls_tmps[0]
            ox = sampls_oxys[0]
            mt = sampls_met[0]
            list_of_tmps = []
            list_of_oxys = []
            list_of_met = []
            y_tmp = None
            y_ox = None
            y_met = None
            for k in lt:
                if abs(k - ht) < step/2:
                    y_tmp = tmp
                    y_ox = ox
                    list_of_tmps.append(y_tmp)
                    list_of_oxys.append(y_ox)
                    list_of_met.append(y_met)
                    break
            # if y_tmp == None and y_ox == None:
            #     continue
            # elif y_tmp == None:
            #     the_dict = {"expedition": one_exp, "stations_num": stat, "lat": new_df["lat"].mean(), "lon": new_df["lon"].mean(), "ht_real": ht, "horiz_depth": new_df["horiz_sampl"].mean(), "oxy_mg": list_of_oxys}
            # elif y_ox == None:
            #     the_dict = {"expedition": one_exp, "stations_num": stat, "lat": new_df["lat"].mean(), "lon": new_df["lon"].mean(), "ht_real": ht, "horiz_depth": new_df["horiz_sampl"].mean(), "tmp": list_of_tmps}
            # else:
            #     the_dict = {"expedition": one_exp, "stations_num": stat, "lat": new_df["lat"].mean(), "lon": new_df["lon"].mean(), "ht_real": ht, "horiz_depth": new_df["horiz_sampl"].mean(), "tmp": list_of_tmps, "oxy_mg": list_of_oxys}
            the_dict = {"expedition": one_exp, "stations_num": stat, "lat": new_df["lat"].mean(), "lon": new_df["lon"].mean(), "ht_real": ht, "horiz_depth": new_df["horiz_sampl"].mean(), "tmp": list_of_tmps, "oxy_mg": list_of_oxys, "ch4_wat": list_of_met}
        else:
            error_oc = False
            a_ht = sampls_heights[0]
            b_ht = sampls_heights[1]
            a_tmp = sampls_tmps[0]
            b_tmp = sampls_tmps[1]
            a_ox = sampls_oxys[0]
            b_ox = sampls_oxys[1]
            a_mt = sampls_met[0]
            b_mt = sampls_met[1]
            y_tmp = None
            y_ox = None
            y_met = None
            list_of_tmps = []
            list_of_oxys = []
            list_of_met = []
            iterator = 1
            for k in lt:
                print("a", a_ht, "k", k, "b", b_ht)
                if k >= a_ht:
                    y_tmp = a_tmp
                    y_ox = a_ox
                    y_met = a_mt
                elif b_ht < k < a_ht:
                    y_tmp = a_tmp + ((k - a_ht)/(b_ht - a_ht)) * (b_tmp - a_tmp)
                    y_ox = a_ox + ((k - a_ht)/(b_ht - a_ht)) * (b_ox - a_ox)
                    y_met = a_mt + ((k - a_ht)/(b_ht - a_ht)) * (b_mt - a_mt)
                elif k == b_ht:
                    y_tmp = b_tmp
                    y_ox = b_ox
                    y_met = b_mt
                else:
                    print("Warning!")
                    while b_ht > k:
                        print("Iteration")
                        a_ht = b_ht
                        iterator += 1
                        try:
                            b_ht = sampls_heights[iterator]
                        except IndexError:
                            list_of_tmps.append(sampls_tmps[iterator-1])
                            list_of_oxys.append(sampls_oxys[iterator-1])
                            list_of_met.append(sampls_met[iterator-1])
                            error_oc = True
                            break
                        else:
                            a_tmp = b_tmp
                            b_tmp = sampls_tmps[iterator]
                            a_ox = b_ox
                            b_ox = sampls_oxys[iterator]
                            a_mt = b_mt
                            b_mt = sampls_met[iterator]
                            print("a", a_ht, "k", k, "b", b_ht)
                    if error_oc:
                        break
                    if b_ht == k:
                        y_tmp = b_tmp
                        y_ox = b_ox
                    else:
                        y_tmp = a_tmp + ((k - a_ht)/(b_ht - a_ht)) * (b_tmp - a_tmp)
                        y_ox = a_ox + ((k - a_ht)/(b_ht - a_ht)) * (b_ox - a_ox)
                        y_met = a_mt + ((k - a_ht)/(b_ht - a_ht)) * (b_mt - a_mt)
                    # except ZeroDivisionError:
                    #     while b_ht == a_ht:
                    #         iterator += 1
                    #         a_ht = b_ht
                    #         b_ht = sampls_heights[iterator]
                    #         a_tmp = b_tmp
                    #         b_tmp = sampls_tmps[iterator]
                    #         a_ox = b_ox
                    #         b_ox = sampls_oxys[iterator]
                    #         print("ZeroDivisionError: (одинаковые высоты)")
                
                list_of_tmps.append(y_tmp)
                list_of_oxys.append(y_ox)
                list_of_met.append(y_met)
            the_dict = {"expedition": one_exp, "stations_num": stat, "lat": lats[:len(list_of_tmps)], "lon": lons[:len(list_of_tmps)], "ht_real": lt[:len(list_of_tmps)], "horiz_depth": lt_2[:len(list_of_tmps)], "tmp": list_of_tmps, "oxy_mg": list_of_oxys, "ch4_wat": list_of_met}
        print("Exp: ", one_exp, "St: ", stat, "list_of_oxys ", list_of_oxys, "\nlist_of_tmps ", list_of_tmps)
        interpolated_station_df = pd.DataFrame(the_dict)
        all_stations_df = pd.concat([all_stations_df, interpolated_station_df], ignore_index=True)
    interpolated_station_df = None
    #all_stations_df = all_stations_df.dropna(subset=['tmp', 'oxy_mg'])
    try:
        unique_levels = list(all_stations_df["ht_real"].unique())
        print(unique_levels)
    except KeyError:
        print("ERRRR!")
        continue
    else:
        pass
    
    # Interpolation for each level, ICW-interpolation
    model_array_temp = None
    layer = 0
    for level in unique_levels:

        # calculating Volume in cubic kilometres
        dem_res1 = dem_res[dem_res != no_data]
        dem_res_new = dem1[dem1 < level]
        dem_sum = (np.sum(level - dem_res_new) * pixel_area)/1000000000

        # calculating Area in square kilometres
        count = len(dem_res_new)
        S = count*pixel_area/1000000

        nan_mask = np.logical_or(np.isclose(dem_res, no_data), ~np.less(dem_res, level))
        copy_array = np.where(nan_mask == True, 0, 0)
        copy_array_1 = np.where(nan_mask == True, 0, 1)
        
        one_level = all_stations_df[all_stations_df["ht_real"] == level]
        print(len(one_level))

        # Delete duplicate rows in the 'lat' and 'lon' columns
        one_level_cleaned = one_level.drop_duplicates(subset=['lat', 'lon'])
        one_level_copy = one_level_cleaned.copy()

        if len(one_level_cleaned) != 1:
            outRaster = driver.Create(f"{path_to_file}/mask1.tif", cols2, rows2, 1, gdal.GDT_Int16)
            outRaster.SetGeoTransform(geotransform_res)
            outRaster.SetProjection(projection_res)

            # Write the modified array to the new GeoTIFF file
            outband = outRaster.GetRasterBand(1)
            outband.WriteArray(copy_array_1)

            # Set NoData value to nan
            outband.SetNoDataValue(0)

            # Close the raster files
            outband.FlushCache()
            outRaster = None
            
            # icw_array_tmp = 0
            # icw_array_oxy = 0
            # array_chisl_tmp = 0
            # array_znam_tmp = 0
            # array_chisl_oxy = 0
            # array_znam_oxy = 0
            print(len(one_level_cleaned))
            print(one_level_cleaned)
            my_array = np.zeros((len(one_level_cleaned), rows2, cols2), dtype=np.float16)

            index_num = 0
            dict_of_temp = {}
            dict_of_oxy = {}
            dict_of_met = {}

            for index, row in one_level_cleaned.iterrows():
                
                point_tmp = row["tmp"]
                dict_of_temp.update({index_num: point_tmp})
                
                point_oxy = row["oxy_mg"]
                dict_of_oxy.update({index_num: point_oxy})
                print(point_oxy)

                point_met = row["ch4_wat"]
                dict_of_met.update({index_num: point_met})

                point_x = row['lon']
                point_y = row['lat']

                # Calculate the column and row indices of the point in the raster array
                inv_geo_transform = gdal.InvGeoTransform(geotransform_res)
                pixel_x, pixel_y = map(int, gdal.ApplyGeoTransform(inv_geo_transform, point_x, point_y))
                
                if rows <= pixel_y or cols <= pixel_x:
                    continue

                print(f"Row {index}: Column index - {pixel_x}, Row index - {pixel_y}")

                outRaster = driver.Create(f"{path_to_file}/mask2.tif", cols2, rows2, 1, gdal.GDT_Int16)
                outRaster.SetGeoTransform(geotransform_res)
                outRaster.SetProjection(projection_res)

                # Write the modified array to the new GeoTIFF file
                outband = outRaster.GetRasterBand(1)
                copy_array[pixel_y, pixel_x] = 10
                outband.WriteArray(copy_array)

                # Close the raster files
                outband.FlushCache()
                outRaster = None
                copy_array[pixel_y, pixel_x] = 0

                source = f'{path_to_file}/mask2.tif'
                cost = f'{path_to_file}/mask1.tif'
                wbt.cost_distance(source, cost, f"{path_to_file}/output2.tif", f"{path_to_file}/backlink.tif", callback = my_callback)
                
                
                target_image = gdal.Open(f"{path_to_file}/output2.tif")
                target_band = target_image.GetRasterBand(1)
                target_array = target_band.ReadAsArray()/1000

                my_array[index_num] = target_array

                print("t_arr", np.dtype(target_array.dtype))
                print("m_arr", np.dtype(my_array.dtype))
                

                target_image = None
                del target_array
                os.remove(f"{path_to_file}/mask2.tif")
                os.remove(f"{path_to_file}/output2.tif")
                

                index_num += 1
                
        
            x_grid, y_grid = np.meshgrid(np.linspace(min_x, max_x, image_ds_res.RasterXSize),
                                    np.linspace(max_y, min_y, image_ds_res.RasterYSize))

            my_array = np.where(nan_mask == True, np.inf, my_array)

            if points == 0:
                z_i = filter_ndarray(my_array, points)
                print("Z индексы:")
                print(z_i)
                temp_values = np.vectorize(dict_of_temp.get)(z_i)
                oxy_values = np.vectorize(dict_of_oxy.get)(z_i)
                print("Массив температуры:")
                print(temp_values)
                temp_interpol = icw_interpolation(temp_values, my_array, tension)
                oxy_interpol = icw_interpolation(oxy_values, my_array/1000, tension)
                print("Массив интерполяции:")
                print(temp_interpol)
                temper_mean = np.nanmean(temp_interpol)
                print(temper_mean)
            else:
                if points <= len(one_level_cleaned):
                    min_indices, min_distances = filter_ndarray(my_array, points)
                else:
                    min_indices, min_distances = filter_ndarray(my_array, len(one_level_cleaned))
                print("Минимальные индексы:")
                print(min_indices)
                print("Минимальные дистанции:")
                print(min_distances)
                temp_values = np.vectorize(dict_of_temp.get, otypes=[np.float16])(min_indices)
                oxy_values = np.vectorize(dict_of_oxy.get, otypes=[np.float16])(min_indices)
                met_values = np.vectorize(dict_of_met.get, otypes=[np.float16])(min_indices)
                print("Массив температуры:")
                print(temp_values)
                print("Массив кислорода:")
                print(oxy_values)
                print("Массив метана:")
                print(met_values)
                if points != 1 and len(one_level_cleaned) != 1:
                    temp_interpol = icw_interpolation(temp_values, min_distances, tension)
                    oxy_interpol = icw_interpolation(oxy_values, min_distances, tension)
                    met_interpol = icw_interpolation(met_values, min_distances, tension)
                    print("Массив интерполяции:")
                    print(temp_interpol)
                    print("Массив интерполяции:")
                    print(oxy_interpol)
                    print("Массив интерполяции:")
                    print(met_interpol)
                    inter_new_temp = temp_interpol[np.isfinite(temp_interpol)]  # Filter out nan and inf values
                    inter_new_oxy = oxy_interpol[np.isfinite(oxy_interpol)] # Filter out nan and inf values
                    inter_new_met = met_interpol[np.isfinite(met_interpol)] # Filter out nan and inf values
                    temper_mean = np.mean(inter_new_temp)
                    oxy_mean = np.mean(inter_new_oxy)
                    met_mean = np.mean(inter_new_met)
                else:
                    temp_interpol = np.where(nan_mask == True, np.nan, temp_values)
                    oxy_interpol = np.where(nan_mask == True, np.nan, oxy_values)
                    met_interpol = np.where(nan_mask == True, np.nan, met_values)
                    inter_new_temp = temp_interpol[np.isfinite(temp_values)]
                    inter_new_oxy = oxy_interpol[np.isfinite(oxy_values)]
                    inter_new_met = met_interpol[np.isfinite(met_interpol)]
                    temper_mean = np.mean(inter_new_temp)
                    oxy_mean = np.mean(inter_new_oxy)
                    met_mean = np.mean(inter_new_met)
        else:
            for index, row in one_level_cleaned.iterrows():
                temper_mean = row["tmp"]
                oxy_mean = row["oxy_mg"]
                met_mean = row["ch4_wat"]
                temp_interpol = np.where(nan_mask == True, np.nan, temper_mean)
                oxy_interpol = np.where(nan_mask == True, np.nan, oxy_mean)
                met_interpol = np.where(nan_mask == True, np.nan, met_mean)
                

        the_dict_1 = {"reservoir": name_reservoir, "expedition": one_exp,  "ht_real": level, "tmp_mean": temper_mean, "oxy_mg": oxy_mean, "ch4_wat": met_mean, "volume": dem_sum, "area": S}
        means = pd.DataFrame(the_dict_1, index=[0])
        all_means = pd.concat([all_means, means], ignore_index=True)

        one_level_copy2 = one_level_copy.copy()
        one_level_final_temp = one_level_cleaned.dropna(subset=['tmp'])
        one_level_final_oxy = one_level_copy.dropna(subset=['oxy_mg'])
        one_level_final_met = one_level_copy2.dropna(subset=['ch4_wat'])

        xs_t = one_level_final_temp['lon'].values
        ys_t = one_level_final_temp['lat'].values

        xs_o = one_level_final_oxy['lon'].values
        ys_o = one_level_final_oxy['lat'].values

        xs_mt = one_level_final_met['lon'].values
        ys_mt = one_level_final_met['lat'].values

        z_values_t = temp_interpol
        z_values_o = oxy_interpol
        z_values_mt = met_interpol
        
        # Визуализация результатов интерполяции

        

        if len(one_level_cleaned) != 1:
            fig = plt.figure()
            fig.set_size_inches(8,8)
            plt.contourf(x_grid, y_grid, z_values_t, levels=number_of_isolines, cmap='viridis')
            plt.colorbar()
            plt.scatter(xs_t, ys_t, color = "b", edgecolors='w')
            plt.xlabel('lon')
            plt.ylabel('lat')
            plt.title(f'reservoir: {name_reservoir}\n exp: {one_exp}\n parameter: T, °C\n method: ICW, {points} points\n level: {one_level_cleaned["ht_real"].mean():.2f} m BS')
            plt.savefig(f"{path_to_file}/{name_reservoir}/{name_reservoir}_{one_exp}_T_ICW_{points}_points_{int(one_level_cleaned['ht_real'].mean())}.png")


            fig = plt.figure()
            fig.set_size_inches(8,8)
            plt.contourf(x_grid, y_grid, z_values_o, levels=number_of_isolines, cmap='viridis')
            plt.colorbar()
            plt.scatter(xs_o, ys_o, color = "r", edgecolors='w')
            plt.xlabel('lon')
            plt.ylabel('lat')
            plt.title(f'reservoir: {name_reservoir}\n exp: {one_exp}\n parameter: Oxygen, mg/l\n method: ICW, {points} points\n level: {one_level_cleaned["ht_real"].mean():.2f} m BS')
            plt.savefig(f"{path_to_file}/{name_reservoir}/{name_reservoir}_{one_exp}_Oxy_ICW_{points}_points_{int(one_level_cleaned['ht_real'].mean())}.png")
        
            fig = plt.figure()
            fig.set_size_inches(8,8)
            plt.contourf(x_grid, y_grid, z_values_mt, levels=number_of_isolines, cmap='Reds')
            plt.colorbar()
            plt.scatter(xs_mt, ys_mt, color = "r", edgecolors='w')
            plt.xlabel('lon')
            plt.ylabel('lat')
            plt.title(f'reservoir: {name_reservoir}\n exp: {one_exp}\n parameter: Methane, mkl/l\n method: ICW, {points} points\n level: {one_level_cleaned["ht_real"].mean():.2f} m BS')
            plt.savefig(f"{path_to_file}/{name_reservoir}/{name_reservoir}_{one_exp}_methane_ICW_{points}_points_{int(one_level_cleaned['ht_real'].mean())}.png")

        try:
            os.remove(f"{path_to_file}/backlink.tif")
        except FileNotFoundError:
            pass
        try:
            os.remove(f"{path_to_file}/mask1.tif")
        except FileNotFoundError:
            pass
        
        
        if layer == 0:
            model_array_temp = temp_interpol.reshape(1, temp_interpol.shape[0], temp_interpol.shape[1])
            model_array_oxy = oxy_interpol.reshape(1, oxy_interpol.shape[0], oxy_interpol.shape[1])
            model_array_meth = met_interpol.reshape(1, met_interpol.shape[0], met_interpol.shape[1])
        else:
            appendant_temp = temp_interpol.reshape(1, temp_interpol.shape[0], temp_interpol.shape[1])
            appendant_oxy = oxy_interpol.reshape(1, oxy_interpol.shape[0], oxy_interpol.shape[1])
            appendant_meth = met_interpol.reshape(1, met_interpol.shape[0], met_interpol.shape[1])
            model_array_temp = np.append(model_array_temp, appendant_temp, axis=0)
            model_array_oxy = np.append(model_array_oxy, appendant_oxy, axis=0)
            model_array_meth = np.append(model_array_meth, appendant_meth, axis=0)

        layer += 1
        # Clear variables
        my_array = None 
        min_indices, min_distances = None, None
        inter_new_temp, temp_interpol, temp_values = None, None, None
        inter_new_oxy, oxy_interpol, oxy_values = None, None, None
        temper_mean, oxy_mean = None, None
        z_values_o, z_values_t = None, None
        copy_array, copy_array_1 = None, None


        # Close Matplotlib figures
        plt.close('all')

model_array_temp.shape
model_array_oxy.shape
model_array_meth.shape


print(getsizeof(model_array_temp)/1024/1024)
print(getsizeof(model_array_oxy)/1024/1024)
print(getsizeof(model_array_meth)/1024/1024)
# model_array_temp.tofile(f"{path_to_file}/model.dat")
np.save(f'{path_to_file}/DRM/test_temp.npy', model_array_temp)
np.save(f'{path_to_file}/DRM/test_oxy.npy', model_array_oxy)
np.save(f'{path_to_file}/DRM/test_meth.npy', model_array_meth)

all_means.to_excel(f"{path_to_file}/Итоговые вертикали/{name_reservoir}_all.xlsx", index = False)