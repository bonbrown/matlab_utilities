function [varh] = read_wrf_height(pth,fl,vnam,ti,hgt,altvarfl) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_wrf_height.m                                                       
%                                                  
% Read a variable from the wrfout file and interpolate to height from
% model native eta coordinates.
%
% inputs:   pth  - path (string)
%           fl   - file name (string)
%           vnam - variable name (string)
%           ti   - time index (integer)
%           hgt  - desired height level in meters (float)
%           altvarfl - alternative filename to source variable of interest
%               from. Still takes height, dimension info from pth+fl, must have same
%               dimensions (for dual pol simulations). Input [], '', or other
%               dummy of length < 2 to bypass (string)
%
% outputs:  varh - the variable field (unstaggered to mass pts) 
%                   interpolated to desired height level (X x Y matrix)
%
%       Bonnie Brown, University of Hawaii
%           October 2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ncfile = strcat(pth,fl);
g = 9.81;
level = hgt;

% WRF variables are output in time by eta by south-north by west-east
% format. MATLAB NetCDF libraries spit them out in opposite format:
% west-east by south-north by eta by time.

% get variable at desired time index
if length(altvarfl) > 2
    dat = ncread(altvarfl,vnam,[1 1 1 ti],[inf inf inf 1]);
else
    dat = ncread(ncfile,vnam,[1 1 1 ti],[inf inf inf 1]);
end
[we,sn,neta,t] = size(dat);

% get size of unstaggered grid
weu = ncreadatt(ncfile,'/','WEST-EAST_PATCH_END_UNSTAG');
snu = ncreadatt(ncfile,'/','SOUTH-NORTH_PATCH_END_UNSTAG');
nulev = ncreadatt(ncfile,'/','BOTTOM-TOP_PATCH_END_UNSTAG');

% unstagger if necessary 
if we>weu
    data = 0.5*(dat(2:end,:,:) + dat(1:end-1,:,:));
elseif sn>snu
    data = 0.5*(dat(:,2:end,:) + dat(:,1:end-1,:));
elseif neta>nulev
    data = 0.5*(dat(:,:,2:end) + dat(:,:,1:end-1));
else 
    data = dat;
end

% get height at given time (read on w pts, compute to mass pts)
phb = ncread(ncfile,'PHB',[1 1 1 ti],[inf inf inf 1]);
ph = ncread(ncfile,'PH',[1 1 1 ti],[inf inf inf 1]) + phb;

hght = ( ph(:,:,2:end) + ph(:,:,1:end-1) ) ./ (2*g) ;

% throw out heights you don't want and interpolate the rest (this snippet
% is borrowed from old Ryan Torn (UWash/UAlbany) script interp_var_hght.m)
varh = zeros(weu,snu);

for ii = 1:weu; for jj = 1:snu

  if ((hght(ii,jj,1) < level) || (hght(ii,jj,nulev) > level))

    % find the nearest height levels
    for kk = 2:nulev
      if hght(ii,jj,kk) > level
        klev = kk-1;
        break;
      end;
    end;

    % linearly interpolate
    m = (data(ii,jj,klev+1) - data(ii,jj,klev)) ./ ...
        (hght(ii,jj,klev+1) - hght(ii,jj,klev)); 
    varh(jj,ii) = m .* (level - hght(ii,jj,klev)) + data(ii,jj,klev);

  else % if height is above or below model, set to NaN

    varh(ii,jj) = NaN;

  end
  
end; end
