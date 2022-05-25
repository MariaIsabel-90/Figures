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

#the days that will be plotted
days_plot=[39520,39558,39573,39576]

#to load the data that dont change on time
#Plotting the eddy number 7774 from 
path = '/data0/project/benguela/BUIC/OUTPUT/R3KM/BENGSAFE_R3KM_V2c_EXTRACTED_eddydetect_fulllong/'
files = ['roms_avg_Y2008M3.nc','roms_avg_Y2008M4.nc','roms_avg_Y2008M5.nc','roms_avg_Y2005M10.nc']
nc2 = Dataset(path+files[0],'r')
lon = nc2.variables['lon_rho'][:] #tem q colocar [:]
lat = nc2.variables['lat_rho'][:]
bathy = nc2.variables['h'][:]
nc2.close()

#the plot
nc3 = Dataset(path+files[0],'r')
time = nc3['scrum_time'][:]
time = time/(24*3600);
zeta = nc3.variables['zeta'][:]
greg = jd_to_date(time,dateref)


#features 
proj=ccrs.PlateCarree()  
land_10m=cfeature.NaturalEarthFeature(category='physical',
                                      name='land',
                                      scale='10m',
                                      edgecolor='k',
                                      facecolor='lightgrey')
for jj in range(len(time)):
        intt=np.floor(time[jj])==np.floor(days_plot)
        if intt.any()==True:
            fig=plt.figure(figsize=(12,10))
            ax1 = plt.plot(projection=proj)
            ax1.add_feature(cfeature.COASTLINE)     #Add Coastlines
            gl1=ax1.gridlines(linestyle=':',draw_labels=True)   #Add Gridlines
            gl1.xlabels_top = False
            gl1.ylabels_right = False
            gl1.xformatter = LONGITUDE_FORMATTER
            gl1.yformatter = LATITUDE_FORMATTER
            gl1.xlabel_style = {'size': 14, 'color': 'black'}
            gl1.ylabel_style = {'size': 14, 'color': 'black'}
            ax1.add_feature(land_10m)
            ctf1 = ax1.contour(lon,lat,zeta[jj,:,:],np.linspace(-1.25,1.25,25),cmap='jet')
            ax1.contour(lon,lat,bathy,cmap='gray')
            ax1.scatter(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],color='k')
            ax1.text(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],str(7774),fontsize=14)
            intt2= np.floor(time[jj])==np.floor(teddy2)
            if intt2.any()==False:
                pass
            else: 
                ax1.scatter(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],color='k')
                ax1.text(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],str(8128),fontsize=14) 
nc3.close()
nc3 = Dataset(path+files[1],'r')
time = nc3['scrum_time'][:]
time = time/(24*3600);
zeta = nc3.variables['zeta'][:]
greg = jd_to_date(time,dateref)
for jj in range(len(time)):
        intt=np.floor(time[jj])==np.floor(days_plot)
        if intt.any()==True:
            ax2 = plt.subplot(2, 2, 2, projection=proj)
            ax2.add_feature(cfeature.COASTLINE)     #Add Coastlines
            gl2=ax2.gridlines(linestyle=':',draw_labels=True)   #Add Gridlines
            gl2.xlabels_top = False
            gl2.ylabels_right = False
            gl2.xformatter = LONGITUDE_FORMATTER
            gl2.yformatter = LATITUDE_FORMATTER
            gl2.xlabel_style = {'size': 14, 'color': 'black'}
            gl2.ylabel_style = {'size': 14, 'color': 'black'}
            ax2.add_feature(land_10m)
            ctf2 = ax2.contour(lon,lat,zeta[jj,:,:],np.linspace(-1.25,1.25,25),cmap='jet')
            ax2.contour(lon,lat,bathy,cmap='gray')
            ax2.scatter(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],color='k')
            ax2.text(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],str(7774),fontsize=14)
            intt2= np.floor(time[jj])==np.floor(teddy2)
            if intt2.any()==False:
                pass
            else: 
                ax1.scatter(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],color='k')
                ax1.text(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],str(8128),fontsize=14) 
nc3.close()

nc3 = Dataset(path+files[2],'r')
time = nc3['scrum_time'][:]
time = time/(24*3600);
zeta = nc3.variables['zeta'][:]
greg = jd_to_date(time,dateref)
for jj in range(len(time)):
        intt=np.floor(time[jj])==np.floor(days_plot)
        if intt.any()==True:
            ax3 = plt.subplot(2, 2, 3, projection=proj)
            ax3.add_feature(cfeature.COASTLINE)     #Add Coastlines
            gl3=ax3.gridlines(linestyle=':',draw_labels=True)   #Add Gridlines
            gl3.xlabels_top = False
            gl3.ylabels_right = False
            gl3.xformatter = LONGITUDE_FORMATTER
            gl3.yformatter = LATITUDE_FORMATTER
            gl3.xlabel_style = {'size': 14, 'color': 'black'}
            gl3.ylabel_style = {'size': 14, 'color': 'black'}
            ax3.add_feature(land_10m)
            ctf3 = ax3.contour(lon,lat,zeta[jj,:,:],np.linspace(-1.25,1.25,25),cmap='jet')
            ax3.contour(lon,lat,bathy,cmap='gray')
            ax3.scatter(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],color='k')
            ax3.text(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],str(7774),fontsize=14)
            intt2= np.floor(time[jj])==np.floor(teddy2)
            if intt2.any()==False:
                pass
            else: 
                ax3.scatter(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],color='k')
                ax3.text(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],str(8128),fontsize=14) 
nc3.close()

nc3 = Dataset(path+files[3],'r')
time = nc3['scrum_time'][:]
time = time/(24*3600);
zeta = nc3.variables['zeta'][:]
greg = jd_to_date(time,dateref)
for jj in range(len(time)):
        intt=np.floor(time[jj])==np.floor(days_plot)
        if intt.any()==True:
            ax4 = plt.subplot(2, 2, 4, projection=proj)
            ax4.add_feature(cfeature.COASTLINE)     #Add Coastlines
            gl4=ax4.gridlines(linestyle=':',draw_labels=True)   #Add Gridlines
            gl4.xlabels_top = False
            gl4.ylabels_right = False
            gl4.xformatter = LONGITUDE_FORMATTER
            gl4.yformatter = LATITUDE_FORMATTER
            gl4.xlabel_style = {'size': 14, 'color': 'black'}
            gl4.ylabel_style = {'size': 14, 'color': 'black'}
            ax4.add_feature(land_10m)
            ctf4 = ax4.contour(lon,lat,zeta[jj,:,:],np.linspace(-1.25,1.25,25),cmap='jet')
            ax4.contour(lon,lat,bathy,cmap='gray')
            ax4.scatter(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],color='k')
            ax4.text(loneddy1[np.floor(teddy1)==np.floor(time[jj])],lateddy1[np.floor(teddy1)==np.floor(time[jj])],str(7774),fontsize=14)
            intt2= np.floor(time[jj])==np.floor(teddy2)
            if intt2.any()==False:
                pass
            else: 
                ax4.scatter(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],color='k')
                ax4.text(loneddy2[np.floor(teddy2)==np.floor(time[jj])],lateddy2[np.floor(teddy2)==np.floor(time[jj])],str(8128),fontsize=14) 
nc3.close()
plt.show()
