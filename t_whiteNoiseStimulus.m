

% load(['/Users/james/Documents/MATLAB/'...
%     'akheitman/NSEM_mapPRJ/Stimuli/'...
%     'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);
% 
% 
% for fr = 1:200
%     frIm(:,:,fr) = imresize(single(testmovie.matrix(:,:,fr)'),[92 121]);
% end
    
movRand = (rand(92,121,movLen));
% movNorm = (movRand./norm(movRand(:)));
% movScaled = ieScale(movNorm,255);

movRandAbs = 1*movRand./sum(abs(movRand(:)));
movRandScaled = 127*movRandAbs./max(movRandAbs(:)) + 128;
movRandNormed = round(normScaled*movRandScaled./norm(movRandScaled(:)));
%%
iStimWN = ieStimulusMovieCMosaic(movRandNormed,coneParams);
cMosaicWN = iStimWN.cMosaic;

%% Bipolar

bpParams.cellType = 'offdiffuse';
bp = bipolar(cMosaicWN, bpParams);
bp.set('sRFcenter',1);
bp.set('sRFsurround',0);
bp.compute(cMosaicWN);

%% RGC
% params.eyeRadius = mmEcc;
% params.eyeAngle = 90;
% innerRetina=ir(bp,params);
% cellType = {'off parasol'};
% innerRetina.mosaicCreate('type',cellType{1});

innerRetina.compute(bp);

psth = innerRetina.mosaic{1}.get('psth');
% figure; ieMovie(psth(:,:,1:10:end));

spWN =  innerRetina.mosaic{1}.get('spikes');
sum(spWN(:))

szSpWN = size(spWN);

spWNAvg = sum(spWN(:))/(frTotal*coneParams.integrationTime*szSpWN(1)*szSpWN(2));
spWNCell = reshape(spWN,[szSpWN(1)*szSpWN(2) szSpWN(3)]);
sumSpWNCell = sum(spWNCell'); stdSpWNlCell = std(sumSpWNCell);

%%
figure; hold on;
bar([1 2 3],[spAvg,spWNAvg,spNSAvg]); %colormap(hot)
errorbar([1 2 3],[spAvg,spWNAvg,spNSAvg],...
    ([stdSpOptimalCell stdSpWNlCell stdSpNSCell]/2)./sqrt(szSpWN(1)*szSpWN(2)),...
    ([stdSpOptimalCell stdSpWNlCell stdSpNSCell]/2)./sqrt(szSpWN(1)*szSpWN(2)),'x','linewidth',4);
grid on
xtickcell{1} = ''; xtickcell{3} = ''; xtickcell{5} = ''; xtickcell{7} = '';
xtickcell{2} = 'OPTIMAL'; xtickcell{4} = 'WHITE NOISE'; xtickcell{6} = 'NATURAL SCENES';
set(gca,'XTickLabel',xtickcell);
xlabel('Category of Stimulus'); ylabel('Mean Rate Per Cell (spikes/sec)');
set(gca,'fontsize',14);
title('Effect of stimulus type on average firing rate');
