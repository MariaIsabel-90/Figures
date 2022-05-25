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
time_all = nc.variables['time'][:]
ID_all = nc.variables['ID'][:]

nc.close()
#separate the eddy trajectories
id_plot_winter=[1109,1182,1288,3506,3818,3841,5962,8128,8560,7774]
#to know the dates and eddy location for generating the maps
#the main vortex
teddy1=time_all[ID_all==id_plot_winter[-1]]
lateddy1=lat_all[ID_all==id_plot_winter[-1]]
loneddy1=lon_all[ID_all==id_plot_winter[-1]]
#the vortex after splitting
teddy2=time_all[ID_all==id_plot_winter[7]]
lateddy2=lat_all[ID_all==id_plot_winter[7]]
loneddy2=lon_all[ID_all==id_plot_winter[7]]

#
# Convert julian date to gregorian date and vice versa
dateref=datetime(1900,1,1)
def jd_to_date(jd, date_ref):
    date = [date_ref + timedelta(days=x) for x in jd]
    return date

def greg_to_jd(dategreg,dateref):
    tmp=dategreg-dateref
    return np.float(tmp.days) + np.float(tmp.seconds)/86400
#the dates in gregorian callendar
gregsd = jd_to_date(time_all,dateref)
gregsdn = np.array(gregsd)

#plot the trajectories of these eddies
#features 
proj=ccrs.PlateCarree()  
land_10m=cfeature.NaturalEarthFeature(category='physical',
                                      name='land',
                                      scale='10m',
                                      edgecolor='k',
                                      facecolor='lightgrey')
#the plot 

path = '/data0/project/benguela/BUIC/OUTPUT/R3KM/BENGSAFE_R3KM_V2c_EXTRACTED_eddydetect_fulllong/'
files =['roms_avg_Y2008M3.nc','roms_avg_Y2008M4.nc','roms_avg_Y2008M5.nc','roms_avg_Y2008M6.nc','roms_avg_Y2008M7.nc','roms_avg_Y2008M8.nc','roms_avg_Y2005M9.nc','roms_avg_Y2005M10.nc']
nc2 = Dataset(path+files[0],'r')
lon = nc2.variables['lon_rho'][:] #tem q colocar [:]
lat = nc2.variables['lat_rho'][:]
bathy = nc2.variables['h'][:]
nc2.close()

#dates for plotting specific maps
juldbeg=39520.0
juldend=39773.0
juld=np.arange(juldbeg,juldend+1)
for ii in range(len(files)):
    nc3 = Dataset(path+files[ii],'r')
    time= nc3['scrum_time'][:]
    time=time/(24*3600);
    zeta = nc3.variables['zeta'][:]
    greg= jd_to_date(time,dateref)
    for jj in range(len(time)):
        intt=np.floor(time[jj])==np.floor(juld)
        if intt.any()==True:
            fig=plt.figure(figsize=(12,10))
            ax=plt.axes(projection=proj)
            ax.add_feature(cfeature.COASTLINE)     #Add Coastlines
            gl=ax.gridlines(linestyle=':',draw_labels=True)   #Add Gridlines
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER
            gl.xlabel_style = {'size': 14, 'color': 'black'}
            gl.ylabel_style = {'size': 14, 'color': 'black'}
            ax.add_feature(land_10m)
            ctf2 = ax.contour(lon,lat,zeta[jj,:,:],np.linspace(-1.25,1.25,25),cmap='jet')
            ax.contour(lon,lat,bathy,cmap='gray')
            ax.scatter(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],color='k')
            ax.text(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],str(7774),fontsize=14)
            intt2= np.floor(time[jj])==np.floor(teddy2)
            if intt2.any()==False:
                pass
            else:
                ax.scatter(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],color='k')
                ax.text(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],str(8128),fontsize=14)                      
            ax.set_title(str(greg[jj]),fontsize=16)
            plt.savefig("/home/mbarros/Documents/traj/Map{dia}.png".format(dia=np.int(time[jj]))) 
            plt.close()
    nc3.close()

