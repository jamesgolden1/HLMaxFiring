
clear srf irf sta sta3

params.eyeRadius = 4;
innerRetina=ir(bp,params);
cellType = {'on parasol'};
innerRetina.mosaicCreate('type',cellType{1});

srf = RGB2XWFormat(innerRetina.mosaic{1}.sRFcenter{1,1}-innerRetina.mosaic{1}.sRFsurround{1,1});
irf= innerRetina.mosaic{1}.tCenter;

sta = srf*irf{1}';

movieBig = zeros(256+64,256+64,size(sta,2));

rs = [4 8 12 16 20 24 ]; 

for ecc = 1:5
    ecc
    clear sta3
    for ii = 1:size(sta,2)
        staTemp = XW2RGBFormat(sta(:,ii),93,93);
        sta3(:,:,ii) = imresize(staTemp,[rs(ecc),rs(ecc)]);
    end
    
    for xc = [-1 0 1]
        for yc = [-1 0 1]
%     for xc = [-1 -.5 0 .5 1]
%         for yc = [-1 -.5 0 .5 1]
            if ~(xc == 0 && yc == 0)
                nv = norm([xc yc]); xc = xc/nv; yc = yc/nv;
                xcvecst = xc*(8+24*ecc) + 128+32 - round(size(sta3,1)/2) + 1;
                xcvecend = xc*(8+24*ecc) + 128+32 + round(size(sta3,1)/2);
                ycvecst = yc*(8+24*ecc) + 128+32 - round(size(sta3,1)/2) + 1;
                ycvecend = yc*(8+24*ecc) + 128+32 + round(size(sta3,1)/2);
                movieBig(round(xcvecst:xcvecend),round(ycvecst:ycvecend),:) = sta3;
            end
        end
    end
    
end

figure; ieMovie(movieBig);