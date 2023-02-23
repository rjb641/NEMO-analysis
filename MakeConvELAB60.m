clc
clear
%/mnt/storage6/clark/ANHA4/Anallysis/MakeConvLAB60.m
CaseConf='2_ANHA4-ECP017'
YearStart=2010
YearEnd=2018
daystart=1
monstart=1
countstart=0
SKIP=1
CountSkip=0
MakeFig=1


CaseConf='2_ANHA4-ECP017'
Mesh_hgr='/mnt/storage6/clark/NEMO_meshes/AGRIF/LAB60/LabSea-domain/2_LAB60_mesh_hgr.nc'
Mesh_zgr='/mnt/storage6/clark/NEMO_meshes/AGRIF/LAB60/LabSea-domain/2_LAB60_mesh_zgr.nc'
Maskfile='/mnt/storage6/clark/NEMO_meshes/AGRIF/LAB60/LabSea-domain/2_LabSeaMask.nc'
Height=2000
   path='/mnt/storage6/clark/ANHA4/ANHA4-ECP017-LAB60/'
LatCenter=55
LonCenter=-45
Radius=20

e3t=GetNcVar(Mesh_zgr,'e3t_1d');
e3t_ps=GetNcVar(Mesh_zgr,'e3t_ps');
g=9.81; % gravity

Tmask=GetNcVar(Maskfile,'tmask');
%find min/max X and Y of this file for smaller opening of files later
[Y X]=size(Tmask)
if strcmp(CaseConf,'2_ANHA4-ECP017')
   Tmask=GetNcVar('/mnt/storage6/clark/NEMO_meshes/AGRIF/LAB60/LabSea-domain/2_LAB60_mask.nc','tmask');
   [Z Y X]=size(Tmask)
   Tmask=squeeze(Tmask(1,:,:));
   size(Tmask)
end
[ignore1 Z]=size(e3t)

Minx=0;
Maxx=0;
Miny=0;
Maxy=0;
for i=1:X
   if(Minx==0)
      if(sum(Tmask(:,i)>=1))
         Minx=i
      end
   end
end
for i=X:-1:1
   if(Maxx==0)
      if(sum(Tmask(:,i)>=1))
         Maxx=i
      end
   end
end

for i=1:Y
   if(Miny==0)
      if(sum(Tmask(i,:)>=1))
         Miny=i
      end
   end
end
for i=Y:-1:1
   if(Maxy==0)
      if(sum(Tmask(i,:)>=1))
         Maxy=i
      end
   end
end
for z=1:Z
   if(sum(e3t(1:z))<Height)
      Zcount=z
   end
end
clear('Tmask');
Minx=Minx
Miny=Miny
Xcount=Maxx-Minx
Ycount=Maxy-Miny
Levels=Zcount



ovVar='hdept'; 
ovfile='/mnt/storage4/clark/ANHA4-ECP007/2_mesh_zgr.nc'; 
conLev=[0 500 1000 1500 2000 2500 3000 3500 4000 4500 5000];
conVar=GetNcVar(ovfile,ovVar,[Minx Miny 0],[Xcount Ycount 1]);
conVar=permute(conVar,[2,1]);

xx=GetNcVar(Mesh_hgr,'nav_lon',[Minx Miny],[Xcount Ycount]);
yy=GetNcVar(Mesh_hgr,'nav_lat',[Minx Miny],[Xcount Ycount]);
MinLat=min(min(yy))
MaxLat=max(max(yy))
MinLon=min(min(xx))
MaxLon=max(max(xx))

LatCenter=(MinLat+MaxLat)/2
LonCenter=(MinLon+MaxLon)/2
Radius=sqrt((MaxLat-MinLat)^2+(MaxLon-MinLon)^2)



xx=permute(xx,[2,1]);
yy=permute(yy,[2,1]);
Length=GetNcVar(Mesh_hgr,'e1t',[Minx Miny 0],[Xcount Ycount 1]);
Width=GetNcVar(Mesh_hgr,'e2t',[Minx Miny 0],[Xcount Ycount 1]);
e3t=GetNcVar(Mesh_zgr,'e3t_1d',[0 0],[Zcount 1]);
e3t_ps=GetNcVar(Mesh_zgr,'e3t_ps',[Minx Miny 0],[Xcount Ycount 1]);

TimeTag=['y',num2str(YearStart),'m11d05'];
File=[path,CaseConf,'_',TimeTag,'_gridT.nc'];

Salinity=GetNcVar(File,'vosaline',[Minx Miny 0 0],[Xcount Ycount Zcount 1]);
Tmask=Salinity;
Tmask(Tmask>0)=1;
Hmask=zeros(size(Tmask));
Depth=zeros(size(Tmask));
%make e3t
[ZZ YY XX]=size(Depth)

for z=1:Levels-1
z
   for i=1:XX
      for j=1:YY
         if(Tmask(z,j,i)==1)
            Depth(z,j,i)=e3t(z); 
         end
         if(and(Tmask(z,j,i)==1, Tmask(z+1,j,i)==0));
            Depth(z,j,i)=e3t_ps(j,i);
         end
      end
   end
end
Tmask(Tmask==0)=nan;

disp('done making depth and tmask')
for z=1:Levels
z
   for i=1:XX
      for j=1:YY
         if(nansum(Depth(1:z,j,i),1)<=Height)
            Hmask(z,j,i)=1;
         end
      end
   end
end

%   for i=1:XX
%      for j=1:YY
%         if(nansum(Depth(:,j,i),1)<=Height)
%            Hmask(:,j,i)=0;
%         end
%      end
%   end
%%
tmask=Tmask(1,:,:);
tmask=permute(tmask,[2,3,1]);
disp('done makeing hmask')

disp('Can loop over YearStart and Year End now')

count=countstart;
count2=0
for year=YearStart:YearEnd
   for mon=monstart:12
      for day=daystart:31
         TimeTag=['y',num2str(year),'m',num2str(mon,'%02d'),'d',num2str(day,'%02d')];
         File=[path,CaseConf,'_',TimeTag,'_gridT.nc'];
         if exist(File,'file')
            File
            count=count+1
            Salinity=GetNcVar(File,'vosaline',[Minx Miny 0 0],[Xcount Ycount Zcount 1]);
            Temperature=GetNcVar(File,'votemper',[Minx Miny 0 0],[Xcount Ycount Zcount 1]);
            PDens=(sw_dens0(Salinity,Temperature)); % get density
            PDens(PDens<1000)=nan;
            MaxDens=PDens.*Hmask.*Tmask;
            MaxDens=max(MaxDens,[],1);
            PotDens=PDens.*Hmask.*Tmask;
            PotDens(PotDens==0)=nan;
            Inner=nansum(PotDens.*Depth);
            %Outer=Height*MaxDens;
            Outer=nansum(Depth,1).*MaxDens;
            Diff=abs(squeeze(Outer-Inner));
            Conv=(g./(Length.*Width.*tmask)).*(Diff.*Length.*Width.*tmask);
            Conv=permute(Conv,[2,1]);
            %Conv(Conv>10000)=nan;
            Conv(Conv==0)=nan;
            min(min(Conv));
            max(max(Conv))
 
            if (MakeFig==1)
            disp('making figure now')
            if(count==1)
                            f=figure('visible','on')
            else
                           f=figure('visible','off') 
            end
            m_proj('lambert','lat',[MinLat MaxLat],'long',[MinLon MaxLon],'rect','off');
            hold on
            %%PH
            v_2=[0 500 1000 1500 2000 2500 3000 3500 4000];
            vmax=4000;
            vmin=0;
            centerPoint=0;
            scalingIntensity=1.5;
            
%             %             m_pcolor(xx,yy,Conv)
             m_pcolor(xx,yy,Conv)
             shading flat
             load nclcolormap
             cmap=nclcolormap.BlAqGrYeOrReVi200;
             colormap(cmap);
             xbar=1:length(cmap);
             xbar=xbar-(centerPoint-vmin)*length(xbar)/(vmax-vmin);
             xbar = scalingIntensity * xbar/max(abs(xbar));
             xbar = sign(xbar).* exp(abs(xbar));
             xbar = xbar - min(xbar); xbar = xbar*3999/max(xbar)+1; 
             newMap = interp1(xbar, cmap, 0:4000);
             colormap(newMap); 
            
            title(['ECP017 Convective Energy: ',TimeTag])
            caxis([0 4000]) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hcb=colorbar
            set(hcb, 'YTick', v_2)
            m_gshhs_i('patch',[0.8 0.8 0.8]); set(findobj('tag','m_gshhs_i'),'linestyle','none')
            m_grid('box','fancy','tickdir','in')
            m_contour(xx,yy,conVar,conLev,'linecolor','k','ShowText','off','linestyle','-','LineWidth',.5);
            set(gcf,'color','w');
            set(gcf, 'InvertHardCopy', 'off');
            savename=['LAB60-fig/LS60_ECP017_ConvE_',num2str(count,'%04d'),'.png']
            eval(['print -dpng -r300 ',savename])
close all
          end %count skip
       end % file exists
     end %day
      daystart=1;
   end %month
   monstart=1;
end %year

clc
clear
%end
