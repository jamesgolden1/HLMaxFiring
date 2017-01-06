
clear srf irf sta sta3

params.eyeRadius = 4;
innerRetina=ir(bp,params);
cellType = {'on parasol'};
innerRetina.mosaicCreate('type',cellType{1});

srf = RGB2XWFormat(innerRetina.mosaic{1}.sRFcenter{1,1}-innerRetina.mosaic{1}.sRFsurround{1,1});
irf= innerRetina.mosaic{1}.tCenter;

sta = srf*irf{1}';

movieBig = zeros(256+64,256+64,3*size(sta,2));

% rs = [4 8 12 16 20 24 ]*2; 
rs = [4:2:24]*2;
eccind = 0;
for ecc = .1:.5:5
    eccind = eccind+1
    clear sta3
    for ii = 1:size(sta,2)
        staTemp = XW2RGBFormat(sta(:,ii),186,186);
        sta3(:,:,ii) = imresize(staTemp,[rs(eccind),rs(eccind)]);
    end
%     
%     for xc = [-1 0 1]
%         for yc = [-1 0 1]
    for xc = [-1:.25:1]
        for yc = [-1:.25:1]
            if ~(xc == 0 && yc == 0)
                xcr = xc + 1*(1/8)*rand(1,1);
                ycr = yc + 1*(1/8)*rand(1,1);
                nv = norm([xcr ycr]); xcn = xcr/nv; ycn = ycr/nv;
                xcvecst = xcn*(8+24*ecc) + 128+32 - round(size(sta3,1)/2) + 1;
                xcvecend = xcn*(8+24*ecc) + 128+32 + round(size(sta3,1)/2);
                ycvecst = ycn*(8+24*ecc) + 128+32 - round(size(sta3,1)/2) + 1;
                ycvecend = ycn*(8+24*ecc) + 128+32 + round(size(sta3,1)/2);
                
                rstart = round((400-1)*rand(3,1))+1;
                for tind = 1%:2
                    movieBig(round(xcvecst:xcvecend),round(ycvecst:ycvecend),rstart:rstart+199) = ...
                        movieBig(round(xcvecst:xcvecend),round(ycvecst:ycvecend),rstart:rstart+199) + sta3;
                end
            end
        end
    end
    
end
p.save = true;
p.vname = '/Users/james/Documents/MATLAB/HLMaxFiring/maxFireStim0.avi';
figure; ieMovie(movieBig,p);