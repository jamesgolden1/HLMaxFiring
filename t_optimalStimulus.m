% t_optimalStimulus
% 
% Build the optimal stimulus for a small patch of RGCs and compute spike
% response.
% 
% Isetbio: git checkout 57f9775 (rgc spatial RF Q)
% Set extent to 5 in buildSpatialRFArray
% Add cm.os.timeStep = cm.integrationTime; to ieStimulusCMosaic

%%
degreesEccentricity = 1;
mmEcc = DegreesToRetinalEccentricityMM(degreesEccentricity);

coneParams.radius = mmEcc;
coneParams.meanLuminance = 100;
% movieInput = rand(10,10,10);
movieInput = zeros(10,10,100);
movieInput(:,:,40) = 255;
coneParams.integrationTime = .005;
iStim = ieStimulusMovieCMosaic(movieInput,coneParams);
cMosaic = iStim.cMosaic;
%% Bipolar

bpParams.cellType = 'offdiffuse';
bp = bipolar(cMosaic, bpParams);
bp.set('sRFcenter',1);
bp.set('sRFsurround',0);
bp.compute(cMosaic);

irfBpBig = RGB2XWFormat(bp.responseCenter)';
irfBp = irfBpBig(:,100); irfBp = irfBp(1:60)./max(irfBp);
figure; plot(.005*[1:60],irfBp)
% irfBp = irfBpBig(:,50); irfBp = irfBp(1:30)./max(irfBp);
% figure; plot(coneParams.integrationTime*[1:30],irfBp)


%% RGC
params.eyeRadius = mmEcc;
params.eyeAngle = 90;
innerRetina=ir(bp,params);
cellType = {'off parasol'};
innerRetina.mosaicCreate('type',cellType{1});

irf= innerRetina.mosaic{1}.tCenter;

innerRetina.compute(bp);

psth = innerRetina.mosaic{1}.get('psth');

% figure; ieMovie(psth(:,:,1:10:end));

%% Build image

rfCoords = vertcat(innerRetina.mosaic{1}.cellLocation{:});
rfMinR = min(rfCoords(:,1)); rfMaxR = max(rfCoords(:,1));
rfMinC = min(rfCoords(:,2)); rfMaxC = max(rfCoords(:,2));

rfSize = size(innerRetina.mosaic{1}.sRFcenter{1,1});

edgePadding = 4;
spStim = zeros(edgePadding+ceil(rfSize(1)/1)+ceil(rfMaxR-rfMinR),edgePadding+ceil(rfSize(2)/1)+ceil(rfMaxC-rfMinC));

for ri = 1:size(innerRetina.mosaic{1}.cellLocation,1)
    for ci = 1:size(innerRetina.mosaic{1}.cellLocation,2)
%         [ri ci]
        rvStart{ri,ci} = 1+ceil(innerRetina.mosaic{1}.cellLocation{ri,ci}(1) +ceil((rfMaxR-rfMinR)/2)+1);% - ceil(rfSize(1)/2)+1);
        rvEnd{ri,ci}   = 1+ceil(innerRetina.mosaic{1}.cellLocation{ri,ci}(1) +ceil((rfMaxR-rfMinR)/2)) + ceil(rfSize(1)/1);
        
        cvStart{ri,ci} = 1+ceil(innerRetina.mosaic{1}.cellLocation{ri,ci}(2) +ceil((rfMaxC-rfMinC)/2)+1);% - ceil(rfSize(2)/2)+1);
        cvEnd{ri,ci}   = 1+ceil(innerRetina.mosaic{1}.cellLocation{ri,ci}(2) +ceil((rfMaxC-rfMinC)/2) + ceil(rfSize(2)/1));
        
        spStim(rvStart{ri,ci}:rvEnd{ri,ci},cvStart{ri,ci}:cvEnd{ri,ci}) = ...
            spStim(rvStart{ri,ci}:rvEnd{ri,ci},cvStart{ri,ci}:cvEnd{ri,ci})+innerRetina.mosaic{1}.sRFcenter{ri,ci}-innerRetina.mosaic{1}.sRFsurround{ri,ci};
        
    end
end

figure; imagesc(spStim);

%% Build movie
movLen = 250;
movStim = zeros(edgePadding+ceil(rfSize(1)/1)+ceil(rfMaxR-rfMinR),edgePadding+ceil(rfSize(2)/1)+ceil(rfMaxC-rfMinC),movLen);

irf2 = irfBp; %irf(1:end);% NOT RIGHT IRF, subsampled here

for ri = 1:size(innerRetina.mosaic{1}.cellLocation,1)
    for ci = 1:size(innerRetina.mosaic{1}.cellLocation,2)
%         [ri ci]
        srf = RGB2XWFormat(innerRetina.mosaic{1}.sRFcenter{ri,ci}-innerRetina.mosaic{1}.sRFsurround{ri,ci});
        % irf= innerRetina.mosaic{1}.tCenter;
        
        sta = srf*irf2'; 

        for ii = 1:size(sta,2)
            staTemp = XW2RGBFormat(sta(:,ii),rfSize(1),rfSize(2));
            sta3(:,:,ii) = staTemp;%imresize(staTemp,[rs(eccind),rs(eccind)]);
        end
        
%         tStart = 1;%round(1+1*160*rand(7,1));
        
        if mod(ri,2)==0 && mod(ci,2)==0
            tStart = [1 40 80 120 160];
        elseif mod(ri,2)==0 && mod(ci,2)==1
            tStart = 20+[1 40 80 120 160];
        elseif mod(ri,2)==1 && mod(ci,2)==0
            tStart = 20+[1 40 80 120 160];
        elseif mod(ri,2)==1 && mod(ci,2)==1
            tStart = [1 40 80 120 160];
        end
        for ti = 1:length(tStart)
        movStim(rvStart{ri,ci}:rvEnd{ri,ci},cvStart{ri,ci}:cvEnd{ri,ci},tStart(ti):tStart(ti)+length(irf2)-1) = ...
            movStim(rvStart{ri,ci}:rvEnd{ri,ci},cvStart{ri,ci}:cvEnd{ri,ci},tStart(ti):tStart(ti)+length(irf2)-1)+sta3;
        end
    end
end

% figure; ieMovie(movStim);

%%
% movRZ = ieScale(movStim(:,:,1:40));
% movOptNorm = (movRZ./norm(movRZ(:)));
% movOptScaled = ieScale(movOptNorm,255);

frTotal = 20;%movLen;
movOptRed = movStim(:,:,1:frTotal);
movOptAbs = 1*movOptRed./sum(abs(movOptRed(:)));
movOptScaled = round(127*movOptAbs./max(movOptAbs(:))) + 128;
normScaled = norm(movOptScaled(:));
% meanRed = mean(movOptRed(:));
% movOptScaled = 

% figure; ieMovie(movOptScaled);
% movOptSum = sum(movOptRed(:));
% movOptNorm = movOptRed./movOptSum;
% movOptScaled = ieScale(movOptNorm,0,255);
% 
movOptScaled(1,1,:) = 1;
clear cMosaic
xofs = 1; yofs = 1;
iStim = ieStimulusMovieCMosaic(movOptScaled(xofs:end-xofs+1,yofs:end-yofs+1,:),coneParams);
cMosaic = iStim.cMosaic;

% cMosaic.integrationTime = .01;

figure; imagesc(cMosaic.absorptions(:,:,25));
% figure; plot(RGB2XWFormat(cMosaic.absorptions)')
% figure; plot(RGB2XWFormat(cMosaic.current)')
%% Bipolar
clear bp
bpParams.cellType = 'offdiffuse';
bp = bipolar(cMosaic, bpParams);
bp.set('sRFcenter',1);
bp.set('sRFsurround',0);
bp.compute(cMosaic);
% figure; plot(RGB2XWFormat(bp.responseCenter)')
%% RGC
% params.eyeRadius = mmEcc;
% params.eyeAngle = 90;
% innerRetina=ir(bp,params);
% cellType = {'off parasol'};
% innerRetina.mosaicCreate('type',cellType{1});
% 
% irf= innerRetina.mosaic{1}.tCenter;

innerRetina.compute(bp);

% figure; plot(RGB2XWFormat(innerRetina.mosaic{1}.responseLinear))

psth = innerRetina.mosaic{1}.get('psth');
% figure; ieMovie(psth(:,:,1:10:end));

spOptimal =  innerRetina.mosaic{1}.get('spikes');
szSpOptimal = size(spOptimal);
spAvg = sum(spOptimal(:))/(frTotal*coneParams.integrationTime*szSpOptimal(1)*szSpOptimal(2));
spOptimalCell = reshape(spOptimal,[szSpOptimal(1)*szSpOptimal(2) szSpOptimal(3)]);
sumSpOptimalCell = sum(spOptimalCell'); stdSpOptimalCell = std(sumSpOptimalCell);

sum(spOptimal(:))