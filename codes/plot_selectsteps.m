clear all 
close all
%This script is based on okubo_aviso, get_one_eddy,eddy_detect_model,
%eddy_tracking_param
roms_dir=['/data0/project/benguela/BUIC/OUTPUT/R3KM/' ...
          'BENGSAFE_R3KM_V2c_EXTRACTED_eddydetect_fulllong/'];
      
model='roms';
filetype='avg';  % 'his' or 'avg'
suffix='.nc';
% Domain limits 
lonmin=8.25;
lonmax=28.8;
latmin=-44.5;
latmax=-24.5;
%
% Duration
%
Ymin=2005;
Ymax=2005;

Mmin=12;
Mmax=12;

dzeta=0.02;   % Interval [m] between the contours 
              % (should be around the precision of altimetry (~2cm ?))
%
Rmax=300;     % Maximum radius [km] of a close curved detected 
              % (to prevent taking an ocean gyre as a giant eddy)
	      % (should be larger than the largest mesoscale eddies: 300-400 km?)
%
Omin=-2e-12;  % Threshold for Okubo-Weiss parameter detection [s-2]
%Omin=0;       % (Chelton (2007) used -2e12 , but here it work also with 0 !...)
%
Nhanning=3;   % Number of Hanning filter pass on the Okubo-Weiss parameter
              % (1 or 2 passes might help a bit still...)

grd_file=[roms_dir,'roms','_',filetype,'_Y',num2str(Ymin),...
                                       'M',num2str(Mmin),suffix];



nc=netcdf(grd_file);
lon=nc{'lon_rho'}(:);
lat=nc{'lat_rho'}(:);
%
imin=find(lon(1,:)>=lonmin, 1 );
imax=max(find(lon(1,:)<=lonmax));
jmin=min(find(lat(:,1)>=latmin));
jmax=max(find(lat(:,1)<=latmax));
%
lon=lon(jmin:jmax,imin:imax);
lat=lat(jmin:jmax,imin:imax);
%

mask=nc{'mask_rho'}(jmin:jmax,imin:imax);
pm=nc{'pm'}(jmin:jmax,imin:imax);
pn=nc{'pn'}(jmin:jmax,imin:imax);
f=nc{'f'}(jmin:jmax,imin:imax);
%read time
time=nc{'scrum_time'}(31);
time=time/(24*3600);
%read the SSH
zeta=squeeze(nc{'zeta'}(31,jmin:jmax,imin:imax));
close(nc)
clear nc roms_dir model file_type suffix grd_file
%%Step 1- Local SSH maximun and minimun

%to get the local maximuns and local minimuns
[Imax,Jmax,LONmax,LATmax,Zmax]=get_locmax(lon,lat,zeta,mask); 
[Imin,Jmin,LONmin,LATmin,Zmin]=get_locmin(lon,lat,zeta,mask);

zeta(mask==0)=mean(zeta(mask==1));
%
% Disregard data on the shelves
%
zeta(mask==0)=NaN;
%
% get the minimum eddy radius detectable
%
[dx,dy]=get_dx(lon,lat);
ds=dx.*dy;
dxmax=max(dx(:)/1000);
dymax=max(dy(:)/1000);
Rmin=max([dxmax dymax]);
DX=mean([mean(dx(:)) mean(dy(:))]);
npts=round(1e3*Rmax/DX);
%
% Get the Okubo-Weiss parameter, Vorticity and Eddy Kinetic Energy
%xi is the relative vorticity
[oku,vort]=okubo_aviso(pm,pn,f,zeta);
vort(isnan(vort))=0;
oku(isnan(oku))=0;
eke=ke_aviso(pm,pn,f,zeta);
U=sqrt(2*eke);
%
% Apply Nhanning times a Hanning filter on the Okubo-Weiss parameter
%
%oku_h1=oku;
oku_h3=oku;
%oku_h5=oku;
%Nhanning=3;
for n=1:Nhanning
  oku_h3=hanning(oku);
%  oku(mask0==0)=0;
end
%
%oku_h1(mask==0)=NaN;
oku_h3(mask==0)=NaN;
%oku_h3(mask==0)=NaN;
oku(mask==0)=NaN;

%To get the local maximum and minimum for Okubo-Weiss
[ImaxO,JmaxO,LONmaxO,LATmaxO,OMAX]=get_locmax(lon,lat,oku_h3,mask); 
[IminO,JminO,LONminO,LATminO,OMIN]=get_locmax(lon,lat,oku_h3,mask); 
%Ploting the local maxima and minimun in SSH
% 
                   
%close all
ismax=1;
i=0;
for nmax=1:length(Imax)
    
    %if isfinite(Imax(nmax))
     %  [Imin,Imax,~,~,~, ~,~,~, ~,~,~,~,~,~]= get_one_eddy(lon,lat,mask,zeta,oku_h3,vort,eke,U,ds,...
      %     Omin,npts,...
       %    Imin,Jmin,LONmin,LATmin,Zmin,...
        %   Imax,Jmax,LONmax,LATmax,Zmax,...
         %  nmax,dzeta,ismax);
        
     %plot_SSH(lon,lat,mask,zeta,oku_h3,vort,eke,U,ds,...
      %                npts,...
       %               Imin,Jmin,Zmin,...
        %              Imax,Jmax,Zmax,...
         %             nmax,dzeta,ismax)
     [xc,yc,Rc]=get_xy2(lon,lat,mask,zeta,eke,U,ds,...
                      npts,...
                      Imin,Jmin,LONmin,LATmin,Zmin,...
                      Imax,Jmax,LONmax,LATmax,Zmax,...
                      nmax,dzeta,ismax);
     
      if isfinite(xc)
          i=i+1;
          xx(1:length(xc),i)=xc;
          yy(1:length(yc),i)=yc;
          Reddy(:,i)=Rc;
      end
    
end
xx(xx==0)=nan;
yy(yy==0)=nan;
%Zmin and Zmax for plot and dz for plot
zmin=-1.00;
dz2=0.1;
zmax=2.00;

figure(1)
 m_proj('mercator',...
          'lon',[lonmin lonmax],...
          'lat',[latmin latmax]); 
 h2=m_contour(lon,lat,zeta,[zmin:dz2:zmax],'k');
 coastfile='coastline_l.mat';
 
 m_usercoast(coastfile,'patch',[.9 .9 .9]);
 hold on;
 title(['SSH Y ',num2str(Ymax),' - M ',num2str(12),...
        ' - D ',num2str(31)],'Fontsize',14)
 hold on;
 Ix=(xx(1,:));
 for i=1:length(Ix)
     m_plot(xx(:,i),yy(:,i),'r')
     hold on;
 end
 hold on;
 










