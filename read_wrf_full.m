function [field] = read_wrf_full(pth,fl,vnam,ti,int,altvarfl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% read_wrf_full.m
%
% read a full, 3-D field from WRF wrfout file for a time interval and unstagger to mass points
% if necessary.
%
%       Bonnie Brown, University of Hawaii
%           October 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncfile = strcat(pth,fl);

% get variable at desired times
% get variable at desired time index
if length(altvarfl) > 2
    dat = ncread(altvarfl,vnam,[1 1 1 ti],[inf inf inf int]);
else
    dat = ncread(ncfile,vnam,[1 1 1 ti],[inf inf inf int]);
end
[we,sn,neta] = size(dat);

% get size of unstaggered grid
weu = ncreadatt(ncfile,'/','WEST-EAST_PATCH_END_UNSTAG');
snu = ncreadatt(ncfile,'/','SOUTH-NORTH_PATCH_END_UNSTAG');
nulev = ncreadatt(ncfile,'/','BOTTOM-TOP_PATCH_END_UNSTAG');

% unstagger if necessary 
if we>weu
    field = 0.5*(dat(2:end,:,:) + dat(1:end-1,:,:));
elseif sn>snu
    field = 0.5*(dat(:,2:end,:) + dat(:,1:end-1,:));
elseif neta>nulev
    field = 0.5*(dat(:,:,2:end) + dat(:,:,1:end-1));
else 
    field = dat;
end
