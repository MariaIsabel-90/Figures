import numpy as np                              #array manipulation and math
from matplotlib import pyplot as plt            #plotting the data
from netCDF4 import Dataset,num2date            #read netcdf
import matplotlib.patches                        
import matplotlib.dates as mpD                  #dates in plots
import cartopy.crs as ccrs                      #mapping module
import cartopy.feature as cfeature              #to put features in a map
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER  #latitude, longitude formats
from datetime import datetime                   #transform dates to python dates
from datetime import timedelta   

nc = Dataset('/net/centaure/local/tmp/1/mbarros/EXP_2005to2008/TRA02_SEL01_DET05_eddies_r3km_2005M1_2008M12.nc','r')

#getting the variables from nc 
lat_all = nc.variables['lat'][:]
lon_all = nc.variables['lon'][:] #tem q colocar [:]
vort_all = nc.variables['Vorticity'][:]
Energy_all = nc.variables['Energy'][:]
Ampl_all = nc.variables['Amplitude'][:]
U_all = nc.variables['U'][:]
V_all = nc.variables['V'][:]
time_all = nc.variables['time'][:]
radius_all = nc.variables['Radius'][:]
ID_all = nc.variables['ID'][:]
nc.close()

path = '/data0/project/benguela/BUIC/OUTPUT/R3KM/BENGSAFE_R3KM_V2c_EXTRACTED_eddydetect_fulllong/roms_avg_Y2005M6.nc'
nc2 = Dataset(path,'r')

#lon, lat and f only need to be read once
lon = nc2.variables['lon_rho'][:] #tem q colocar [:]
lat = nc2.variables['lat_rho'][:]
zeta = nc2.variables['zeta'][:]

nc2.close()

#contruct the distance matrix
#to calculate the distances in a sphere

def dist_spheric(lat1,lon1,lat2,lon2): # function that takes angles in deg
    #return mastrix dist of size NI*ND
    #NI: size of lat2
    #ND size of lat1
    # arguments must be vectors (NOT 2D arrays)
    R= 6371.008*10**3
    NI=np.size(lat2)
    ND=np.size(lat1)

    l=np.abs(np.dot(np.transpose(lon2),np.ones((1,ND)))- np.dot(np.ones((NI,1)),lon1))
    l[l>=180]=360-l[l>=180]

    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)
    l=np.radians(l)

    dist=R*np.arctan2(np.sqrt((np.sin(l)*np.dot(np.transpose(np.cos(lat2)),np.ones((1,ND))))**2 \
                        +(np.dot(np.transpose(np.sin(lat2)),np.cos(lat1)) \
                         - np.dot(np.transpose(np.cos(lat2)),np.sin(lat1))*np.cos(l))**2 \
                        ),
                          np.dot(np.transpose(np.sin(lat2)),np.sin(lat1))+\
                          np.dot(np.transpose(np.cos(lat2)),np.cos(lat1))*np.cos(l) \
                          )

    return dist

ID_inridge=[];

#to look for each ID if they at some point got close to the ridge
#first we list the list of different eddies IDs from eddies bigger than 100km#
IDs=np.unique(ID_all)
#now we investigate each eddy, choosing only the ones that got close to the ridge. The criteria was the maximum radius, in order to not exclude the biggest eddies
for id in IDs:
    loneddy=np.array([lon_all[ID_all==id]])                 #separating the positions of the eddy of an specific ID
    lateddy=np.array([lat_all[ID_all==id]])   
    dist_ridge=dist_spheric(lateddy,loneddy,lat_line,lon_line)    #calculating the distance between each position
    inridge=np.array(np.where(dist_ridge<=np.max(radius_all)))
    if inridge.size==0:
        pass
    else:
        ID_inridge.append(id)
    
id_50=[]

for id in ID_inridge:
    raeddy=radius_all[ID_all==id]   #calculating the distance between each position
    is50=np.array(np.where(raeddy>50000))
    if is50.size==0:
        pass
    else:
        id_50.append(id)
#plot the trajectories of these eddies
#features 
proj=ccrs.PlateCarree()  
land_10m=cfeature.NaturalEarthFeature(category='physical',
                                      name='land',
                                      scale='10m',
                                      edgecolor='k',
                                      facecolor='lightgrey')
#the plot in loop
#for id in range(len(id_50)):
 #   fig=plt.figure(figsize=(12,10))
  #  ax=plt.axes(projection=proj)
   # #ax.add_feature(cfeature.COASTLINE)     #Add Coastlines
    #gl=ax.gridlines(linestyle=':',draw_labels=True)   #Add Gridlines
    #gl.xlabels_top = False
    #gl.ylabels_right = False
    #gl.xformatter = LONGITUDE_FORMATTER
    #gl.yformatter = LATITUDE_FORMATTER
    #gl.xlabel_style = {'size': 14, 'color': 'black'}
    #gl.ylabel_style = {'size': 14, 'color': 'black'}
    #ax.add_feature(land_10m)
    #Let's plot Bathymetry
    #ctf = ax.contour(lon,lat,h,cmap='gray')
    #loneddy=np.array(lon_all[ID_all==id_50[id]])                                  #separating the positions of the eddy of an specific ID
    #lateddy=np.array(lat_all[ID_all==id_50[id]]) 
    #ax.plot(loneddy,lateddy)
    #id_traj=id_50[id]
    #plt.savefig("Map_traj_{id_traj}.png".format(id_traj=np.int(id_traj)))
                                      #separating the positions of the eddy of an specific ID


