clear all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% radar_pdf.m
%   create a bivariate histogram (a pdf) of radar variables such as
%   Zh and Zdr.
%
%       Bonnie Brown 23 December 2014, University of Hawaii
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Only goes to 5.5 km (hardcoded) to stay below freezing level (we only do
% ZDR for liquid water in WRF at the moment.

%%%%%%%%%%%%%%%%% User Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cl = 'r';
%wrfpth = '/Users/brbrown/WRF/ana2014/';  % path to WRF output%
%wrfpth = '/santaana/WRF/arthur2014/';
wrfpth = '/santaana/WRF/Ana2014/';
%dppth = '/Users/mmbell/Development/dualpol_sim/'; %path to simulated dual pol products
%dppth = '/Users/brbrown/WRF/ana2014/specbin/';
%dppth = '/santaana/WRF/Ana2014/specbin/';
dppth = [wrfpth 'wsm6/'];


dpfl = 'dualpol_radar.nc';
%dpfl = 'arthur_morrison.nc';
%
% REFL_10CM, ZDR, ZV, ZH all for d04
%

mcphys = 'wsm6/';                          % microphysics subdirectory
%flnm = 'wrfout_d04_2014-07-03_12:00:00';    % wrfout filename
flnm = 'wrfout_d04_2014-10-19_00:00:00';
varr1 = 'ZDR';   % variable name in ncf radar files
varr2 = 'REF';
varw1 = 'ZDR';   % variable name in WRF
varw2 = 'ZH';
dualpol = 1;    %read wrf variable from simulated dual pol file?

domask = 0;
mask = 1; % 1 for convective, 0 for stratiform
load /Users/brbrown/Documents/MATLAB/cs_partition_wsm6.mat


vnam1 = 'Differential Reflectivity Z_{DR} (dB)';     % variable string
vnam2 = 'Horizontal Reflectivity Z_{h} (dBZ)';

%st = 36;    %starting time index
%et = 48;    %ending time index
st = 60; et = 72;

% for gridded radar observations %
useradar = 0;
rdr = 'KLTX';
dt = '20140703';
%obspth = ['/Users/brbrown/data/arthur2014/radar/' rdr '/output/' dt '/' dt '/'];
obspth = '/zephyr/andrew41/Thesis/Volumes/Ana/PHKI-3hr/ParamFile3/Composite/';
radht = [.25:.25:14]*1000; % height grid in meters
nl = length(radht);
%obsbase = ['/Users/brbrown/data/ana2014/radar/PHKI/20141019/'];
obsbase = '/Users/brbrown/matlab_scripts/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

km = 1000;
g = 9.81;

if ~useradar
    pth = strcat(wrfpth,mcphys);
    titlestr = strcat('ARTHUR initialized 2014-07-03 12:00:00 PDF for ',mcphys,' hours ',num2str(st*15/60),' to ',num2str(et*15/60));
    tmp1 = read_wrf_full(pth,flnm,varw1,st,et-st,[dppth dpfl]);
    tmp2 =  read_wrf_full(pth,flnm,varw2,st,et-st,[dppth dpfl]);
    data_zdr = tmp1(:,:,1:11,:);
    data_zh = tmp2(:,:,1:11,:);
    %{
    tt=0;
    for ti = st:et
        tt=tt+1;
        data_zdr(:,:,1,tt) = read_wrf_height(pth,flnm,'ZDR',ti,2000,[dppth dpfl]);
        data_zh(:,:,1,tt) = read_wrf_height(pth,flnm,'REFL_10CM',ti,2000,[]);
    end
    %}
    if domask
    for t = 1:et-st
            csmask = csmask_write(:,:,st+(t-1));
            for l = 1:13
                tz = data_zh(:,:,l,t);
                tzd = data_zdr(:,:,l,t);
                tz(csmask~=mask) = NaN;
                tzd(csmask~=mask) = NaN;
                data_zh(:,:,l,t) = tz;
                data_zdr(:,:,l,t) = tzd;
            end
    end
    end
elseif useradar
   %  titlestr=strcat('ARTHUR initialized 2014-07-03 12:00:00 PDF for ',rdr,' radar, hours ',num2str(st*15/60),' to ',num2str(et*15/60));
    titlestr='ANA 19 Oct 2014 1430 to 1730 UTC - PHKI, PHMO radar composite';
     % first get filenames of all gridded radar files in directory (should
    % be only times you want)
    %{
    [names,nf] = get_radar_filenames(rdr,dt,obsbase,'grid');
    for f=1:nf 
        tmp1 = ncread([obspth names(f,:)],varr1,[1 1 1],[inf inf inf]);
        tmp2 = ncread([obspth names(f,:)],varr2,[1 1 1],[inf inf inf]);
        data_zdr(:,:,:,f) = tmp1(:,:,1:16);
        data_zh(:,:,:,f) = tmp2(:,:,1:16);
        
        if domask
            csmask = ncread([obspth names(f,:)],'CSMASK',[1 1],[inf inf]);
            for l = 1:22
                tz = data_zh(:,:,l,f);
                tzd = data_zdr(:,:,l,f);
                tz(csmask~=mask) =  NaN;
                tzd(csmask~=mask) = NaN;
                data_zh(:,:,l,f) = tz;
                data_zdr(:,:,l,f) = tzd;
            end
        end
        
    end
    %}
    tmp1 = ncread('radar_mosaic_3.nc',varr1);
    tmp2 = ncread('radar_mosaic_3.nc',varr2);
    data_zdr = tmp1(:,:,1:16,:);
    data_zh = tmp2(:,:,1:16,:);
    
else
    disp('Select to use radar observations or WRF data please using the boolean variable useradar');
end

[nx,ny,nz,nt] = size(data_zdr);
zdr = reshape(data_zdr,nx*ny*nz*nt,1);
zh = reshape(data_zh,nx*ny*nz*nt,1);
[counts,bins] = hist3([zh zdr],'Edges',{[-10:5:60] [-1:0.5:3.5]});
%[counts,bins] = hist3([zh zdr],'Edges',{[-7.5:5:62.5] [-.75:0.5:3.75]});


tot = max(max(counts));
%tot = sum(sum(counts));
freq = 100*counts./tot;

figure(2), hold on

subplot(3,2,1);

%[c,h] = contour(bins{1},bins{2},freq',[0.5],cl,'linewidth',3);
%clabel(c,h,'labelspacing',288,'color',cl,'fontsize',7)

contourf(bins{1},bins{2},freq',[5:5:100],'linestyle','none')
hold on
contour(bins{1},bins{2},freq',[50 50],'k','linewidth',2);
caxis([5 100]);
set(gca,'box','on','linewidth',2,'layer','top','fontweight','bold','fontsize',12)
[c,h]=contour(bins{1},bins{2},freq',[.01 .1 1 2.5],'color',[0.6 0.6 0.6],'linewidth',2);
clabel(c,h,'color',[0.6 0.6 0.6],'labelspacing',288)

xlabel(vnam2); ylabel(vnam1);
%{
cf = [0.1 0.2 0.5 1 2 5 10:10:100];

cb = colorbar;
caxis([-2.3026 4.6052])
%caxis([-2.3026 3.4012])
contourf(bins{1},bins{2},log(freq'),log(cf));
get(gca); set(gca,'linewidth',2,'fontweight','bold','fontsize',12,'box','on');

set(colorbar,'linewidth',2,'YTick',log([0.1 1 5 20 50 80]),'YTickLabel',{'0.1','1','5','20','50','80'});
%}

