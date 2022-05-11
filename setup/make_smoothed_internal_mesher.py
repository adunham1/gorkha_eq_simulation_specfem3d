
# coding: utf-8

# ## Make topo file for internal mesher

# In[36]:

#imports 
import pandas as pd
import numpy as np
from scipy.interpolate import griddata, RectBivariateSpline
import os
import matplotlib.pyplot as plt
from pyproj import Proj

myProj = Proj("+proj=utm +zone=45, +north +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

#define bounds of model
xmin = 211900.00
xmax = 442100.00
ymin = 2999900.00
ymax = 3185000.00
minlon,minlat = myProj(xmin,ymin,inverse = True)
maxlon, maxlat = myProj(xmax,ymax,inverse = True)

print(minlon,minlat,maxlon, maxlat)


# In[38]:

get_ipython().run_cell_magic('bash', '', "grd='/core2/amd95/SPECFEM3D/himalaya_landslides_specfem3d/_data/20km_smooth_fullsource.grd'\ngrd2xyz $grd > /core2/amd95/SPECFEM3D/himalaya_landslides_specfem3d/_data/20km_smooth_fullsource.xyz")


# In[39]:

res = 1000
grd = '/core2/amd95/SPECFEM3D/himalaya_landslides_specfem3d/_data/20km_smooth_fullsource.xyz'

df = pd.read_csv(grd, delim_whitespace=True, names = ['x','y','z'])
df = df[df.z > 0]


# In[40]:

xrange = np.arange(xmin,xmax+res, res)
yrange = np.arange(ymin,ymax+res, res)

X,Y = np.meshgrid(xrange,yrange)

topo_grd = griddata((df.x.values,df.y.values),df.z.values, (X,Y))

df_grd = pd.DataFrame()
for i in range(len(topo_grd)):
    df_tmp = pd.DataFrame({'x':xrange,'y':[yrange[i]]*len(xrange),'z':topo_grd[i]})
 # df_tmp = pd.DataFrame({'x':xrange,'y':[yrange[i]]*len(xrange),'z':topo_grd[i]-2000}) #use this when you need to make the 2000m smoothed surface for the 500m mesh
    df_grd = df_grd.append(df_tmp)
#plt.imshow(topo_grd)
#plt.colorbar()

print('NX for topo: ' + str(df_grd.x.nunique()))
print('NY for topo: ' + str(df_grd.y.nunique()))
vsmin = 3560
depth = 21000
#target frequency:
f0 = 2
#target element size
#elem = vsmin/f0/1.5
#or use element size of choosing#
elem = 1000
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

#np.savetxt('interface_topo' + str(res) + 'm_larger.dat',np.c_[df_grd.z])
np.savetxt('interface_20kmsmooth_' + str(res) + 'm_fullsource.dat',np.c_[df_grd.z])

