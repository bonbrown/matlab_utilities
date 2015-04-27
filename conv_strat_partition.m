clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% conv_strat_partition.m
% to partition pixels in radar observations of hurricanes into convective, 
% stratiform and weak echo. Based on Appendix A of Didlake and Houze (MWR 2009),
% referred to here as DH; initially leaving tuned parameters to their 
% values, though they used airborne ELDORA observations, may not be appropriate for ground based NEXRAD.
% 
%   Bonnie Brown, November 2014
%    University of Hawaii, Manoa
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------- INPUT OPTIONS ------------------------------------%
disp(['Level of calculation (in km)'])
level = 2 % low-level altitude to conduct partitioning at (km)
disp(['Arbitrary parameters for the convective center calculation'])
a = 9     % arbitrary parameter
b = 45    % arbitrary parameter
disp(['Reflectivity threshold for convective points'])
Zti = 42  % convective threshold intensity (dBZ)
disp(['Reflectivity threshold for weak echoes'])
Zwe = 20  % weak echo threshold (dBZ)
km = 1000;

% path and name of file
%{
nc = '/Users/brbrown/wrf/arthur2014/specbin/wrfout_d04_2014-07-03_12:00:00';
rootp = '/Users/brbrown/wrf/arthur2014/specbin/';
%rootp = '/santaana/WRF/arthur2014/wdm6/';
flnm = 'wrfout_d04_2014-07-03_12:00:00';
%nc = strcat(rootp,flnm);
%}

nc = '/Users/brbrown/matlab_scripts/radar_mosaic_3.nc';
rootp = '/Users/brbrown/matlab_scripts/';
flnm = 'radar_mosaic_3.nc';

%----------------------- WRF OUTPUT FIELDS ---------------------------%
%{
REF = ncread(nc,'REFL_10CM');
lat  = ncread(nc,'XLAT');
lon = ncread(nc,'XLONG');

s = size(REF);
sl = size(lat);
ti = s(4);
reflec = zeros(sl);

disp('Interpolating WRF reflectivity to height level')
for t = 1:ti
    reflec(:,:,t) = read_wrf_height(rootp,flnm,'REFL_10CM',t,level*km,[]);
end
%}
REF = ncread(nc,'REF');
reflec = squeeze(REF(:,:,level*4,:));
lat = ncread(nc,'lat0');
lon = ncread(nc,'lon0');
sl = size(reflec);
ti = sl(3);

% preallocate final product
csmask_write = zeros(sl);

for t = 1:ti
%---------------------------- TIME LOOP ----------------------------------%
    %% Background Reflectivity Calculation

	% Zbg (background reflectivity; dBZ) is average of nonnegative and nonzero reflectivity 
	% within a radius of 11km around the grid point
	disp('Beginning background reflectivity calculation')
	refl = squeeze(reflec(:,:,t));
    sr = size(refl);
	Zbg = NaN*ones(sr);
	for n = 1:sr(1)
	    for m = 1:sr(2)
    	    %dist = haversine(lat(n,m,t),lon(n,m,t),lat(:,:,t),lon(:,:,t));  % find great circle distance from each point
            dist = haversine(lat(n,m),lon(n,m),lat(:,:),lon(:,:));
        	tmp = refl(dist<=11.0);                % find only points within 11km
        	Zbg(n,m) = mean(tmp(tmp>0));           % take mean of points within radius only if reflectivity is nonnegative and nonzero
    	end
	end

	disp('Background reflectivity calculation complete')
    %% Convective Center Calculation
    
    % Now define the convective center criterion delta Zcc first introduced by Steiner et al 1995. The cosine function used
	% by DH was introduced in Yuter and Houze (1997). If Z exceeds Zbg by delta Zcc, it is a convective center.

	dZcc = a*cos( (1/b) * (pi.*Zbg/2) );
	delZ = refl - Zbg;
	cc = find(delZ>=dZcc);

	disp('Convective center calculation complete')
    %% Convective Radius Calculation
   
	% define the convective radius R - this is the radius of points around a convective center which are also classified as convective
	R = zeros(sr);
	R(Zbg<20) = 0.5;
	R(Zbg>=20 & Zbg<35) = 0.5 + 3.5 * (Zbg(Zbg>=20 & Zbg<35) - 20)./15;
	R(Zbg>=35) = 4.0;      
	
	disp('Convective radius calculation complete')
    %% Classification 
    
    % Classify all convective centers, points within a convective radius, and points exceeding the convective intensity threshold 
	% as convective points (1). Classify points beneath the Zwe threshold as weak echoes (2), and everything else as stratiform (0). 
	% Missing data should remain missing (-9999)
	%lat1 = squeeze(lat(:,:,t));
	%lon1 = squeeze(lon(:,:,t));
    lat1 = lat;
    lon1 = lon;
	
	csmask = zeros(sr);
	csmask(cc) = 1;
	ri = R(cc);
	for r = 1:length(ri)
	    dist = haversine(lat1(cc(r)),lon1(cc(r)),lat1,lon1);  % find great circle distance from each point
	    pts = find(dist<=ri(r));
	    csmask(pts) = 1;
	end
	csmask(refl>=Zti) = 1;
	csmask(refl<Zwe) = 2;
	csmask(refl<-900) = refl(refl<-900); %keep missing missing
	
	csmask_write(:,:,t) = csmask;
	% the remaining points stay at a value of 0 (zero) indicating stratiform
    c = clock;
	disp('classification into convective, weak echo, stratiform, and missing complete')
    disp(['FINISHED TIME ' num2str(t) ' at Nov ' num2str(c(3)) ' ' num2str(c(4)) ':' num2str(c(5))])
%------------------------- END TIME LOOP ---------------------------------%	
end
save cs_partition_kltx_kmhx_composite.mat csmask_write 