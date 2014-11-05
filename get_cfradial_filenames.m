function [cfradfl,nf] = get_cfradial_filenames(radar,date,basepath)

pth = [basepath radar '/output/' date '/'];
cfrfl = fopen([pth 'filenames_cfrad.txt'],'r');

c = 1;
cln = fgetl(cfrfl);

while ischar(cln)
        cfradfl(c,:) = cln;
        %disp(cln)
        c = c + 1;
        cln = fgetl(cfrfl);
end

fclose(cfrfl);
nf = c-1;