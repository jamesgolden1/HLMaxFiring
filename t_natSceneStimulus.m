

load(['/Users/james/Documents/MATLAB/'...
    'akheitman/NSEM_mapPRJ/Stimuli/'...
    'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);

% % Needs this repo: https://github.com/isetbio/RemoteDataToolbox
% rdt = RdtClient('isetbio');
% rdt.crp('resources/data/rgc');
% data = rdt.readArtifact('testmovie_schemeA_8pix_Identity_8pix', 'type', 'mat');
% testmovie = data.testmovie;
                
%%

for fr = 1:movLen
    movNS(:,:,fr) = imresize(single(testmovie.matrix(:,:,fr)'),[92 121]);
end
    

% movNSNorm = (movNS./norm(movNS(:)));
% movScaled = ieScale(movNSNorm,255);


movNSRed = movNS(:,:,1:movLen);
movNSAbs = 1*movNSRed./sum(abs(movNSRed(:)));
movNSScaled = 127*movNSAbs./max(movNSAbs(:)) + 128;
movNSNormed = round(normScaled*movNSScaled./norm(movNSScaled(:)));
%%

coneParams.integrationTime = .005;
% iStim = ieStimulusMovieCMosaic(255*ieScale(single(movNS)),coneParams);
iStimNS = ieStimulusMovieCMosaic(movNSNormed,coneParams);
cMosaicNS = iStimNS.cMosaic;
% figure; plot(RGB2XWFormat(cMosaicNS.absorptions)')
% figure; plot(RGB2XWFormat(cMosaicNS.current)')
% figure; ieMovie(cMosaicNS.absorptions);
% figure; ieMovie(cMosaicNS.current);
%% Bipolar
clear bp
bpParams.cellType = 'offdiffuse';
bp = bipolar(cMosaicNS, bpParams);
bp.set('sRFcenter',1);
bp.set('sRFsurround',0);
bp.compute(cMosaicNS);
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

psth = innerRetina.mosaic{1}.get('psth');
% figure; ieMovie(psth(:,:,1:10:end));

% figure; plot(RGB2XWFormat(innerRetina.mosaic{1}.responseLinear))

spNS =  innerRetina.mosaic{1}.get('spikes');
sum(spNS(:))

szSpNS = size(spNS);

spNSAvg = sum(spNS(:))/(frTotal*coneParams.integrationTime*szSpNS(1)*szSpNS(2));
spNSCell = reshape(spNS,[szSpNS(1)*szSpNS(2) szSpNS(3)]);
sumSpNSCell = sum(spNSCell'); stdSpNSCell = std(sumSpNSCell);