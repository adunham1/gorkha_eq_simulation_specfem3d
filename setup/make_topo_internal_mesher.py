
#code that will generate interface_topo500m_fullsource.dat needed for mesh generation but is too large for github

#imports 
import pandas as pd
import numpy as np
from scipy.interpolate import griddata, RectBivariateSpline
import os
import matplotlib.pyplot as plt
from pyproj import Proj

myProj = Proj("+proj=utm +zone=45, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

#define desired bounds of mesh, must be within the DEM files

xmin = 211900.00
xmax = 442100.00
ymin = 2999900.00
ymax = 3185000.00

#convert to lat lon just to better visualize bounds (I think about the world in lat/lon not UTM)
minlon,minlat = myProj(xmin,ymin,inverse = True)
maxlon, maxlat = myProj(xmax,ymax,inverse = True)

print(minlon,minlat,maxlon, maxlat)


# In[67]:
#use gmt to convert 
get_ipython().run_cell_magic('bash', '', "BOUNDS='-R84.3/86.5/27/28.8'\ngrd='/core2/amd95/SPECFEM3D/himalaya_landslides_specfem3d/500mtopo_hpc_parameters/_data/him_dem_500m_fullsource.grd'\ngrd2xyz $grd > /core2/amd95/SPECFEM3D/himalaya_landslides_specfem3d/500mtopo_hpc_parameters/_data/him_500m_utm_gridforinternalmesh_full.xyz")


# In[68]:

res = 500 #define grid size

#name of grid file we just made in xyz format
grd = '/core2/amd95/SPECFEM3D/himalaya_landslides_specfem3d/500mtopo_hpc_parameters/_data/him_'+ str(res) + 'm_utm_gridforinternalmesh_full.xyz'

#load in xyz
df = pd.read_csv(grd, delim_whitespace=True, names = ['x','y','z'])
df = df[df.z > 0] #gets rid of nan values with value of -99999

#make x and y ranges with sampling equal to res and dimensions of mesh
xrange = np.arange(xmin,xmax+res, res)
yrange = np.arange(ymin,ymax+res, res)

#make mesh grid and grid the values with griddata function
X,Y = np.meshgrid(xrange,yrange)
topo_grd = griddata((df.x.values,df.y.values),df.z.values, (X,Y))

#convert the 2d array into dataframe with values in correct order for topo file
df_grd = pd.DataFrame()
for i in range(len(topo_grd)):
    df_tmp = pd.DataFrame({'x':xrange,'y':[yrange[i]]*len(xrange),'z':topo_grd[i]})
    df_grd = df_grd.append(df_tmp)
#plt.imshow(topo_grd)

#save file!
np.savetxt('/core2/amd95/SPECFEM3D/himalaya_landslides_specfem3d/500mtopo_hpc_parameters/DATA/meshfem3d_files/interface_topo' + str(res) + 'm_fullsource.dat',np.c_[df_grd.z])


# In[60]:

#print useful information for meshing
print('NX for topo: ' + str(df_grd.x.nunique()))
print('NY for topo: ' + str(df_grd.y.nunique()))
vsmin = 3560
depth = 21000
#target frequency:
f0 = 2
#target element size
#elem = vsmin/f0/1.5
#or use element size of choosing#
elem =500
nx = (max(df_grd.x)-min(df_grd.x))/elem
ny = (max(df_grd.y)-min(df_grd.y))/elem
nz = depth/elem
nproc = 8
print('Target frequency: ' + str(f0) + 'hz')
print('Target element size: '+ str(elem) + 'm')
print('NX: ' + str(nx))
print('NY: ' + str(ny))
print('NZ: ' + str(nz))
print('Number processors: ' + str(nproc))

print('Use value for NEX_XI: '+ str(round(nx/nproc) *nproc))
print('Use value for NEX_ETA: '+ str(round(ny/nproc) *nproc))


# In[ ]:



