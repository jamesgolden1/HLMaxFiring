% t_viveOptimal
% 
% Build the optimal RGC stimulus for one of the HTC Vive

% Vive specs:
% 1080x1200 per eye, 2160x1200 total
% 2160 pixels = 110 degrees
% 19.6364 pixels per degree
% 0.0509 degrees per pixel

% For midgets, 1 pixel per degree under 10 degrees ecc


%%%%%%%%%% TO DO
% This is dendritic field - add factor for RF?
% Check that this is ON midget type spacing

% Typically we are concerned with the spacing within one class, in which
% case density is halved and the spacings should be multiplied by sqrt2.

% scale factor = ['on parasol', 'off parasol', 'on midget', 'off midget', 'sbc]
% scale factor = [18.92 0.85*18.9211 10.76 0.85*10.7629 20]
% scale midget = mean([10.76 0.85*10.7629]) = 9.95
% scaleFactor = [18.92 0.85*18.9211 10.76 0.85*10.7629 20]./9.95
%% Choose cell type
% This sets the size of the STAs for individual cells

% cellType = 'on parasol';
cellType = 'off parasol';
% cellType = 'on midget';
% cellType = 'off midget';
% cellType = 'on sbc';

% Spatial scale factor by type
% scaleFactorArr = ['on parasol', 'off parasol', 'on midget', 'off midget', 'sbc]
scaleFactorArr = [18.92 0.85*18.9211 10.76 0.85*10.7629 20]./9.95;

%% Add paths - isetbio, HLMaxFiring, RemoteDataToolbox
% addpath(genpath(isetbioRootPath))

%% Initialize movie

% Set size of movie, based on Vive specs:
% 1080x1200 per eye, 2160x1200 total

zeroPad = 0;
szCols = 1080+zeroPad; szRows = 1200+zeroPad;
frames = 250;

% Binocular field of view (FOV)
fovCols = 110; % horizontal fov in degrees for vive: http://doc-ok.org/?p=1414
fovRows = 122;

% % movieBig = 128*uint8(ones(szRows,szCols,frames));
% movieBig = single(zeros(szRows,szCols,frames));

% 2160 pixels = 110 degrees
% 19.6364 pixels per degree
pixelsPerDegree = szCols/fovCols;
% 0.0509 degrees per pixel
degreesPerPixel = fovCols/szCols;

%% Get temporal response from physiology data

% Use the RemoteDataToolbox
% rdt = RdtClient('isetbio');

% Get data for chosen cell type
switch (cellType)
    case{'on parasol'}
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_ONParasol', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
        load('dat/mosaicGLM_apricot_ONParasol.mat');
        scaleFactor = scaleFactorArr(1);
    case{'off parasol'}
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_OFFParasol', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
        load('dat/mosaicGLM_apricot_OFFParasol.mat');
        scaleFactor = scaleFactorArr(2);
    case{'on midget'}
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_ONMidget', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
        load('dat/mosaicGLM_apricot_ONMidget.mat');
        scaleFactor = scaleFactorArr(3);
    case{'off midget'}
%         rdt.crp('resources/data/rgc/apricot');      
%         data = rdt.readArtifact('mosaicGLM_apricot_OFFMidget', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
        load('dat/mosaicGLM_apricot_OFFMidget.mat');
        scaleFactor = scaleFactorArr(4);
    case{'on sbc','sbc'}
%         rdt.crp('resources/data/rgc/apricot');
%         data = rdt.readArtifact('mosaicGLM_apricot_sbc', 'type', 'mat');
%         mosaicGLM = data.mosaicGLM;
        load('dat/mosaicGLM_apricot_sbc.mat');
        scaleFactor = scaleFactorArr(5);
end

% Take mean impulse response function (IRF) over all the cells in the mosaic
for i = 1:length(mosaicGLM)
    irfTemp = mosaicGLM{i}.linearfilters.Stimulus.time_rk1;
    [mv,mi] = max(abs(irfTemp)); irfSign(i) = sign(irfTemp(mi));
    
    irf(i,:) = irfSign(i)*mosaicGLM{i}.linearfilters.Stimulus.time_rk1;
    
end
figure; plot([0:-1+size(irf,2)]/60.5,mean(irf));
irfMean = mean(irf);

% Make the peak IRF negative for off cells
switch (cellType)
    case{'offparasol','offmidget'}
        irfMean = -irfMean;
end

%% Test spatial response
rfRadius = 5;
[so, spatialRFonedim, magnitude1STD] = buildSpatialRF(rfRadius);
% magnitude1STD
% figure; imagesc(so);

%% Get lookup table for how large RF is 
% There is an asymmetry in the size of RGC RFs over the retina

% CHECK SUPERIOR/INFERIOR
rgcDiameterLUT = scaleFactor*sqrt(2)*watsonRGCSpacing(szCols,szCols,fovRows)';

% degStart = -27.5; degEnd = 27.5;
% degarr = [-27.5 : (degEnd-degStart)/1080 : 27.5];
% contourf(degarr,degarr,rgcDiameterLUT',[0:max(rgcDiameterLUT(:))/20:max(rgcDiameterLUT(:))] ); axis square
% title(sprintf('Human Midget RGC RF Size (degrees)')); colorbar; 

%% Add to big movie

% Initialize movie
movieBig = single(zeros(szRows,szCols,frames));

eccind = 0;

% Set center - Vive display is 55 degrees per eye
degCenterX = fovRows/2; degCenterY = fovCols/2;

% Set array of eccentricity values at which to put STA stimulus
eccArr = [1:1:fovCols/2 - 2];
% eccArr = [.1:.25:5 5:.5:10 10:1:27.5];

% Scale the sRF appropriately by cell type
switch (cellType)
    case{'onparasol'}
        sRFmultFactor = 2;
    case{'offparasol'}
        sRFmultFactor = 1.8;
    case{'onsbc','sbc'}
        sRFmultFactor = 2.2;
    otherwise 
        sRFmultFactor = 1;
end

% Loop over eccentricities, put STA stim at individual positions

nAngles = 128/2;
angleNoise = 0;%0.5;
eccNoise   = 0;%2.5;


for ecc = eccArr;
    eccind = eccind+1; 

    % Display how much has been done, how much to go
    [eccind length(eccArr)]    
    
    % xc = ecc*sin(theta), yc = ecc*cos(theta)
%     for xc = [-1:.25/1:1]
%         for yc = [-1:.25/1:1]

    for theta = [0 : 2*pi / nAngles : 2*pi - 2*pi/nAngles]

        % Convert theta to x,y coordinates with some angular noise
        xc = cos(theta) + angleNoise*(1/2)*(rand(1,1)-.5); 
        yc = sin(theta) + angleNoise*(1/2)*(rand(1,1)-.5);
            if ~(xc == 0 && yc == 0)

                % Add some noise to the position
                xcr = xc;% + 0*ecc*(1/2)*(rand(1,1)-.5);
                ycr = yc;% + 0*ecc*(1/2)*(rand(1,1)-.5);
                nv = norm([xcr ycr]); 
                xcn = xcr/nv; ycn = ycr/nv;

                % Get position in pixels instead of degrees
                degX = degCenterX + ecc*xcn;%ecc*(1 + eccNoise*(1/2)*(rand(1,1)-.5));%*xcn; 
                pixelX = degX*pixelsPerDegree + degX*eccNoise*(1/2)*(rand(1,1)-.5);                
                degY = degCenterY + ecc*ycn; %ecc*(1 + eccNoise*(1/2)*(rand(1,1)-.5));%ecc*ycn; 
                pixelY = degY*pixelsPerDegree + degY*eccNoise*(1/2)*(rand(1,1)-.5);

                % Get RGC RF diameter from LUT
                LUTshift = (szRows - size(rgcDiameterLUT,1)+1)/2;
                %disp(-LUTshift+round(pixelX));
                %disp(round(pixelY));
                if (-LUTshift+round(pixelX)) > 0
                    rgcDiameterDegrees = sRFmultFactor*rgcDiameterLUT(-LUTshift+round(pixelX),round(pixelY));
                else
                    rgcDiameterDegrees = sRFmultFactor*rgcDiameterLUT(1,round(pixelY));
                end
                rgcDiameterPixels = rgcDiameterDegrees*pixelsPerDegree;

                rgcRadiusPixels = rgcDiameterPixels/2;

                % Build spatial RF according to diameter
                % If less than one pixel, only use one pixel
                % if ecc>5/(scaleFactor*sqrt(2))
                if rgcDiameterDegrees > degreesPerPixel
                    [so, spatialRFonedim, magnitude1STD] = buildSpatialRF(rgcRadiusPixels);
                    so = so./max(abs(so(:)));
                else
                    so = 1;
                end

                % Add temporal resopnse
                sta = so(:)*irfMean;
                % sta3 = round(128+127* (reshape(sta,[size(so,1),size(so,2),size(irfMean,2)]))./max(abs(sta(:))) );
                sta3 = ( (reshape(sta,[size(so,1),size(so,2),size(irfMean,2)])) );

                % Get start and end position values for adding in STA
                % stimulus
                xcvecst = pixelX - ceil(size(sta3,1)/2) + 1;
                xcvecend = pixelX + floor(size(sta3,1)/2);
                ycvecst = pixelY - ceil(size(sta3,1)/2) + 1;
                ycvecend = pixelY + floor(size(sta3,1)/2);

                % Choose random starting time points for stimulus
                rstart = round((200-1)*rand(12,1))+1;

                if eccind < (length(eccArr) - 2)
                    tind = 1;%:2
                        movieBig(round(xcvecst:xcvecend),round(ycvecst:ycvecend),rstart:rstart+length(irfMean)-1) = ...
                            movieBig(round(xcvecst:xcvecend),round(ycvecst:ycvecend),rstart:rstart+length(irfMean)-1) + (sta3);
                end


%                     % Choose random starting time points for stimulus
%                 rstart = round((200-15)*rand(5,1))+1;
%                 
%                 for tind = 1:length(rstart)
%                     movieBig(round(xcvecst:xcvecend),round(ycvecst:ycvecend),rstart(tind):rstart(tind)+length(irfMean)-1) = ...
%                         movieBig(round(xcvecst:xcvecend),round(ycvecst:ycvecend),rstart(tind):rstart(tind)+length(irfMean)-1) + (sta3);
%                 end

            end
    %         end
    end
    
end

figure; imagesc(sum(movieBig,3)); colormap gray; axis equal

%%
moviePiece = movieBig(:,:,50+[1:50]);
maxMovie = max(abs(moviePiece(:)));
medMovie = median(moviePiece(:));
clear moviePiece

movieSmall = uint8(128 + 127*movieBig/maxMovie);

movieSmall(end,end,:) = 128; movieSmall(end-1,end,:) = 128;
figure; imagesc(sum(movieSmall,3)); colormap gray; axis equal

% clear movieBig

%% Show movie and save

% p.save = false;% 
p.save = true;
% p.vname = 'C:/Users/laha/Documents/GitHub/HLMaxFiring/test.avi';
p.vname = 'test.avi';
p.FrameRate = 90;
figure; 
% set(gcf,'position',[1000         157        1411        1181]);

%disp('test1')
%ieMovie(movieBig(:,:,1:100), p);
%ieMovie(movieBig, p);
%disp('test2')

% ieMovie(movieBig(:,:,1:100));
% ieMovie(movieSmall(:,:,1:100));

movieSmall(end,end,:) = 0; movieSmall(end-1,end,:) = 255;


if p.save
    vObj = VideoWriter(p.vname);
    vObj.FrameRate = p.FrameRate;
    vObj.Quality = 100;
    open(vObj);
end

for fr = 1:10
%     imagesc(movieBig(:,:,fr)); colormap gray    
    imagesc(movieSmall(:,:,fr)); colormap gray; axis image; set(gca,'xticklabel','','yticklabel','');
    if p.save,  F = getframe(gca,[0 0 szRows szRows]); writeVideo(vObj,F); end
%     drawnow;
end

% Write the video object if save is true
if p.save
    writeVideo(vObj,F);
    close(vObj);
end
