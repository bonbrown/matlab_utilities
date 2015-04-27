clear all; clc;
%close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  cfad.m
%   create a cfad (contour frequency by altitude diagram, see
%       Yuter & Houze (1995, MWR), Houze et al (2007, QJRMS).
%
%   Bonnie Brown 19 Feb 2013 Univ. of Washington
%   Heavily modified BRB October 2014 Univ. of Hawaii
%       (uses new MATLAB NetCDF lib, interpolates to height, took out
%       annulus sorting for now)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% User Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%wrfpth = '/Users/brbrown/WRF/arthur2014/';  % path to WRF output
wrfpth = '/santaana/WRF/Ana2014/';
%dppth = '/Users/mmbell/Development/dualpol_sim/'; %path to simulated dual pol products
dppth = '/santaana/WRF/Ana2014/morrison/';
dpfl = 'arthur_wsm6.nc';
%dpfl = 'dualpol_radar.nc';
%
% REFL_10CM, ZDR, ZV, ZH all for d04
%
domask = 0;
mask = 0; % 1 for convective, 0 for stratiform
%load /Users/brbrown/Documents/MATLAB/cs_partition_specbin.mat
load cs_partition_kltx_kmhx_composite.mat;


mcphys = 'morrison/';                          % microphysics subdirectory
flnm = 'wrfout_d04_2014-10-19_00:00:00';    % wrfout filename
%varr = 'ZDR';   % variable name in ncf radar files
%varw = 'REFL_10CM';   % variable name in WRF
varw = 'W';
%varw = 'ZDR';
dualpol = 0;    %read wrf variable from simulated dual pol file?

%vnam = 'vertical velocity (m s^{-1})';     % variable string
%vnam = 'reflectivity (dBz)';
vnam = 'differential reflectivity (dB)';

% variable bin size (not the contour size that you'll see on the CFAD)%
%bins = [-15:5:65];     % standard for reflectivity
%bins = [-5:0.2:5];      % for differential reflectivity
bins = [-20:0.5:20];    % for vertical velocity

%cf = [2.5:2.5:100];    % frequency contours for plotting%
%cf1 = [2.5:2.5:100]; % for cfad method 1
cf1 = [0.1 0.5 1 2.5 5:5:100];
cf0 = [0.1 0.2 0.5 1 2 5 10 15 20 25 30];    % for cfad method 2
st = 36;    %starting time index
et = 48;    %ending time index
toplev = 22; %max height level for WRF (5.2km ~ 14)

cfad_method = 0; % 0 for Yuter & Houze (frequency normalized by # pts in level), 1 for Houze et al (frequency normalized by max frequency at any level)

% for gridded radar observations %
useradar = 0;
rdr = 'KMHX';
dt = '20140703';
obspth = ['/Users/brbrown/data/arthur2014/radar/' rdr '/output/' dt '/' dt '/'];
radht = [.25:.25:14]*1000; % height grid in meters
radht = radht(1:22);
nl = length(radht);
obsbase = ['/Users/brbrown/data/arthur2014/radar/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%counts = zeros(nl,length(
km = 1000;
g = 9.81;

if ~useradar
    pth = strcat(wrfpth,mcphys);
    titlestr = strcat('ARTHUR initialized 2014-07-03 12:00:00 CFAD of ',vnam,' for ',mcphys,' hours ',num2str(st*15/60),' to ',num2str(et*15/60));
    % create a regular height grid using the base geopotential (though it is
    % stored in a 4-D array, it only varies in height) - (unstagger while computing height)
    % --> consider later going to even spaced 1km height bins
    phb = ncread([pth flnm],'PHB',[1 1 1 st],[inf inf inf et-st]);
    rhght = ( phb(:,:,2:end,:) + phb(:,:,1:end-1,:) ) ./ (2*g) ;
    rhght = rhght(:,:,1:toplev,:);
    rht = squeeze(rhght(50,50,:,1));

    [nx,ny,nz,nt] = size(rhght);
    ht = rht;
    % Get variable of interest at regular height levels, calculate frequency of
    % values at each level using specified bins

    for l = 1:length(rht)
        disp(['regular height level ' num2str(l) ' located at ' num2str(rht(l)) ' meters']);
        t0 = 0;
        for t = st:et
            t0 = t0+1;
            csmask = csmask_write;
            if dualpol %read dual pol variables from different file, they are not included in wrfout
                tmp = read_wrf_height(pth,flnm,varw,t,rht(l),[dppth dpfl]);
                if domask
                    tmp(csmask(:,:,t)~=mask) = NaN;
                end
                data(:,:,l,t0) = tmp;
            else
                tmp = read_wrf_height(pth,flnm,varw,t,rht(l),[]);
                if domask 
                    tmp(csmask(:,:,t)~=mask) = NaN;
                end
                data(:,:,l,t0) = tmp;
            end
        end
        % if you want a time series, move this to inside the t loop and modify
        % indices
        tmp = reshape(squeeze(data(:,:,l,:)),nx*ny*(nt+1),1);
        counts(l,:) = histc(tmp,bins);
    end
elseif useradar
    titlestr=strcat('ARTHUR initialized 2014-07-03 12:00:00 CFAD of ',vnam,' for ',rdr,' radar, hours ',num2str(st*15/60),' to ',num2str(et*15/60));
    % first get filenames of all gridded radar files in directory (should
    % be only times you want)
    [names,nf] = get_radar_filenames(rdr,dt,obsbase,'grid');
    %{
    for f=1:nf %analog of time loop but no height interpolation needed (for now, will be on different height grid)
        reflt = ncread([obspth names(f,:)],varr,[1 1 1 1],[inf inf inf inf]);
        csmask = ncread([obspth names(f,:)],'CSMASK',[1 1],[inf inf]);
        for l=1:nl
            tmp = reflt(:,:,l);
            if domask
                tmp(csmask~=mask) = NaN;
            end
            data(:,:,l,f) = tmp;
        end
    end
    clear data
    %}
    %csmask= csmask_write;
    data = ncread('radar_mosaic_3.nc',varr);
    data = data(:,:,1:22,:);
    [ny,nx,nl,nf] = size(data);
    
    for l =1:nl
        
        if domask     
            for f = 1:nf
            tmp = data(:,:,l,f);
            tmp(csmask(:,:,f)~=mask) = NaN;
            data(:,:,l,f) = tmp;
            end
        end

    tmp1 = reshape(data(:,:,l,:),nx*ny*nf,1);
    counts(l,:) = histc(tmp1,bins);
    end
    
    ht = radht;
else
    disp('Select to use radar observations or WRF data please using the boolean variable useradar');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if cfad_method
    % This method normalizes by max frequency at ALL/ANY level (Houze et al 2007)
    tot = max(max(counts));
    freq = 100*counts./tot;
    cf = cf1;
else
    % This method normalizes by total frequency (number of data points) at each level (Yuter & Houze 1995)
    tot = sum(counts,2); % sum over each level
    TOT = tot*ones(1,length(bins));
    freq = 100*counts./TOT; % get percentage of points at each level in each variable value bin
    cf = cf0;
end

% take log of data and change scale to highlight extreme (<1%)
% probabilities
%{
figure(1), 
subplot(2,2,4); hold on

contourf(bins,ht./km,freq,[2.5:2.5:30],'linestyle','none');
set(gca,'box','on','linewidth',2,'layer','top','fontweight','bold','fontsize',12)
[c,h]=contour(bins,ht./km,freq,[.01 .1 1],'color',[0.6 0.6 0.6],'linewidth',2);
clabel(c,h,'color',[0.6 0.6 0.6],'labelspacing',288)
xlabel(vnam); ylabel('height (km)');
title(titlestr);
%}
%{
cb = colorbar;
%caxis([-2.3026 4.6052])
caxis([-2.3026 3.4012])
[h,c] = contourf(bins,ht./km,log(freq),log(cf));
get(gca); set(gca,'linewidth',2,'fontweight','bold','fontsize',12,'box','on');

set(colorbar,'linewidth',2,'YTick',log([0.1 0.2 0.5 1 2 5 10 15 20 25 30]),'YTickLabel',{'0.1','0.2','0.5','1','2','5','10','15','20','25','30'});

xlabel(vnam);
ylabel('height (km)');
title(titlestr);
%set(gca,'XLim',[-10 10],'XTick',[-10:1:10]);
%}
