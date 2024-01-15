# To be typed cmd: python H5_2_Tif_coversion_sample.py filepath
# Load necessary packages into Python
import h5py, glob, sys, getopt, argparse, re, os
import numpy as np
# Import gdal_array to match numpy data type names to gdal data type names
import osgeo
from osgeo import gdal, gdal_array, osr
from fnmatch import fnmatch
#------------------------------------------------------------------------------
# Define Script and handle errors
def main(argv):
    parser = argparse.ArgumentParser()
    try:
        opts, args = getopt.getopt(argv,"hi:",["input_directory"])   
        if len(sys.argv[1:])==0:
            class MyParser(argparse.ArgumentParser):
                def error(self, message):
                    sys.stderr.write('error: %s\n' % message)
                    self.print_help()
                    sys.exit(2)
            parser=MyParser()
            # Add command line argument for input directory            
            parser.add_argument('input_directory', nargs='+')
            args=parser.parse_args()
            # below are a series of potential common errors and responses
        elif "'" in sys.argv[1] or '"' in sys.argv[1]:
            parser.error('error: Do not include quotes in input dir argument')
        elif len(sys.argv) > 2:
            parser.error('error: Only 1 Argument is allowed (input_directory)')
        elif sys.argv[1][-1] != '/' and sys.argv[1][-1] != '\\':
            parser.error('error: Please end directory location with / or \\')
    except getopt.GetoptError:
        print('error: Invalid option passed as argument')      
        print('HDF5toGeoTIFF.py <input_directory>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('HDF5toGeoTIFF.py <input_directory>')
            sys.exit()
    try:
        os.chdir(sys.argv[1])
    except FileNotFoundError:
        print('error: input_directory does not exist or cannot be found')
        sys.exit(2)
#------------------------------------------------------------------------------
    # Set input/current working directory from user defined argument
    in_dir = sys.argv[1]
    
    # Create and set output directory
    out_dir = os.path.normpath((os.path.split(in_dir)[0] + os.sep + 
        'output_py/'))+ '\\'
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
                
    # Create a list of All VIIRS Surface Reflectance HDF-EOS5 files in the dir
    file_list = glob.glob(in_dir + '3DIMG**.h5')
    if len(file_list) == 0:
        print('Error: no valid VNP09 h5 files were found in this directory')
        sys.exit(2)
    # The projection information can not be obtained as a WKT of Proj4 string 
    # from the VIIRS file. Proj info is hard coded to match proj of MODIS tile.
    # projInfo[0] = Sinusoidal info, projInfo[1] = CMG (geo) info
    projInfo = 'PROJCS["unnamed",GEOGCS["Unknown datum based upon the custom spheroid", DATUM["Not specified (based on custom spheroid)", SPHEROID["Custom spheroid",6371007.181,0]],PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]], PROJECTION["Sinusoidal"],PARAMETER["longitude_of_center",0],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]',\
               'GEOGCS["Unknown datum based upon the Clarke 1866 ellipsoid", DATUM["Not specified (based on Clarke 1866 spheroid)", SPHEROID["Clarke 1866",6378206.4,294.9786982139006]], PRIMEM["Greenwich",0], UNIT["degree",0.0174532925199433]]'
    format = "GTiff"
    xmin=44.5
    xmax=105.5
    ymin=-10
    ymax=45.5
    gdal.UseExceptions()
    sourcepath=r'filepath\filename.tif'
    SourceDS = gdal.Open(sourcepath, gdal.GA_ReadOnly)
    #--------------------------------------------------------------------------
    # Function to get geoinformation from the StructMetadata object
    def GetGeographicInfo(input_file):
        # Get info from the StructMetadata Object
        f_Metadata = np.array(input_file['Projection_Information'])
        f_Metadata_byte2str = [s.decode('utf-8') for s in f_Metadata if isinstance('upper_left_lat_lon',s.decode('utf-8'))]
        ulcLon = float(f_Metadata_byte2str.split(',',2)[1])
        ulcLat = float(f_Metadata_byte2str.split(',',2)[0])
        return((ulcLon,  0, ulcLat, 0))
     #-------------------------------------------------------------------------
    # Function to read all datasets in the VIIRS HDF-EOS5 file
    def GetDatasetList(input_file):
        all_h5_objs = []
        input_file.visit(all_h5_objs.append)
        all_datasets = [str(obj) for obj in all_h5_objs if \
                        isinstance(f[obj],h5py.Dataset)]
        return(all_datasets)

    #--------------------------------------------------------------------------
    # Batch process all files in input directory
    for vnp in file_list:
        # Maintain original filename convention    
        vnp_name = vnp.split('\\')[-1][:-3]  
        # Read in the VIIRS HDF-EOS5 file
        f = h5py.File(vnp, "r")
        # Retrieve Geolocation information for the file    
        #geoInfo = GetGeographicInfo(f)
        # Retrieve list of VIIRS datasets    
        dsList = GetDatasetList(f)
        print('Processing: {}.h5'.format(vnp_name))
    #--------------------------------------------------------------------------
        # Loop through each dataset in the file and output as GeoTIFF
        for ds in dsList:
            dsName = ds.split('/')[-1]
            # Create array and read dimensions           
            dsArray=np.array(f['INS'])
            imr_1,nRow,nCol=dsArray.shape
            dsArray.shape=(nRow,nCol)
            lat=np.array(f['Y'])
            lon=np.array(f['X'])
            # # Cell size not specified in the metadata of VIIRS version 001
            xmin,ymin,xmax,ymax=[np.min(lon),np.min(lat),np.max(lon),np.max(lat)]
            yres = (ymax-ymin)/float(nRow)   
            xres = (xmax-xmin)/float(nCol)
            # xres = 4000.458635
            # yres = 3997.98769
            dataType = gdal_array.NumericTypeCodeToGDALTypeCode(dsArray.dtype)
            geotransform=(xmin,xres,0.0,ymax,0.0,-yres)       
                    
            # Set output coordinate referense system information 
            Projection = osr.SpatialReference()
            Projection.ImportFromWkt(SourceDS.GetProjectionRef())           
            driver = gdal.GetDriverByName(format)
            out_ds = driver.Create('{}/{}_{}.tif'.format(out_dir, vnp_name, \
                        dsName), nCol, nRow, 1, dataType)
            out_ds.SetGeoTransform(geotransform)
            out_ds.SetProjection(Projection.ExportToWkt())
            out_ds.GetRasterBand(1).WriteArray(dsArray)
            out_ds = None
            del dataType, ds, dsArray, nCol, nRow, xres, yres
        print('Output location: {}'.format(out_dir))
#------------------------------------------------------------------------------
if __name__ == "__main__":
   main(sys.argv[1:])
#------------------------------------------------------------------------------   
