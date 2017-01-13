% t_maxFiring
% 
% This tutorial generates an RGC mosaic response to a movie stimulus.
%
% The response is generated using the isetbio framework. The general steps 
% are as follows:
% 
%   * display
%   * scene
%   * optical image [oi]
%   * cone mosaic   [absorptions and current]
%   * bipolar mosaic
%   * inner retina  [rgc mosaics]
% 
% The scene is set to be a natural image displayed on a particular 
% calibrated monitor.  The optical image structure computes the image that 
% falls on the retina. The cone mosaic object generates a cone mosaic
% of the appropriate size and computes cone absorptions and current.  The 
% bipolar object generates a bipolar mosaic and computes the bipolar cell 
% responses. The RGC responses are then computed from the bipolar mosaic.

% Need to revert to this commit for isetbio!
% git checkout 27993f6
%% Display

display = displayCreate('LCD-Apple');

%% Scene
% Generic scene
scene = sceneCreate;
figure; sceneShowImage(scene);

% Alternatively, load a custom scene
Ibig = imread('peppers.png');
I = rgb2gray(imresize(Ibig,0.5));
params.meanLuminance = 200;
scene = sceneFromFile(I, 'rgb', params.meanLuminance, display);
%% Optical image
oi  = oiCreate('wvf human');

% Set moving bar parameters
frameTotal = 80;

params.fov = 1.5;
params.barWidth = 10;
params.meanLuminance = 200;
sceneRGB = zeros([sceneGet(scene, 'size'), frameTotal, 3]);

%% Cone mosaic
cMosaic = coneMosaic;

params.fov    = 0.6;
params.radius = 12e-3;
params.theta  = 0;
params.side   = 'left';

% compute cone packing density
fLength = oiGet(oi, 'focal length');
eccMM = 2 * tand(params.radius/2) * fLength * 1e3;
coneD = coneDensity(eccMM, [params.theta], params.side);
coneSz(1) = sqrt(1./coneD) * 1e-3;  % avg cone size with gap in meters
coneSz(2) = coneSz(1);

% Set the cone aperture size
cMosaic.pigment.width  = coneSz(1); 
cMosaic.pigment.height = coneSz(2);

% Set cone mosaic field of view to match the scene
scene = sceneSet(scene,'fov', params.fov);
sceneFOV = [sceneGet(scene, 'h fov') sceneGet(scene, 'v fov')];
sceneDist = sceneGet(scene, 'distance');
cMosaic.setSizeToFOV(sceneFOV, 'sceneDist', sceneDist, 'focalLength', fLength);

% Set the exposure time for each step
cMosaic.integrationTime = .001;%cMosaic.os.timeStep;
%%

for frameNumber = 1:frameTotal
    
    barMovie = ones([sceneGet(scene, 'size'), 3])*0.5;
    
    % Bar at this time
    colStart = frameNumber + 1;
    colEnd   = colStart + params.barWidth - 1;
    % barMovie(:,t-startFrames + 1:(t-startFrames+1+params.barWidth-1),:) = 1;
    barMovie(:,colStart:colEnd,:) = 5;
    
    scene = sceneFromFile(barMovie, 'rgb', params.meanLuminance, display);
    
    scene = sceneSet(scene, 'h fov', params.fov);

    % Get scene RGB data
    sceneRGB(:,:,frameNumber,:) = sceneGet(scene,'rgb');
    
    % Compute optical image
    oi = oiCompute(oi, scene);
    
    % Compute absorptions and photocurrent
    % cMosaic.compute(oi, 'append', true, 'emPath', [0 0]);
    cMosaic.compute(oi, 'emPath', [0 0]);
    if frameNumber == 1; absorptionsMat = zeros([size(cMosaic.pattern) frameTotal]); end;
    absorptionsMat(:,:,frameNumber) = cMosaic.absorptions;
end
cMosaic.absorptions = absorptionsMat;
cMosaic.emPositions=zeros(frameTotal,2);
cMosaic.os.noiseFlag = 'none';
cMosaic.computeCurrent();

%% Bipolar

bpParams.cellType = 'offdiffuse';
bp = bipolar(cMosaic, bpParams);
bp.set('sRFcenter',1);
bp.set('sRFsurround',0);
bp.compute(cMosaic);

%% RGC
params.eyeRadius = 4;
params.eyeAngle = 90;
innerRetina=ir(bp,params);
cellType = {'on parasol'};
innerRetina.mosaicCreate('type',cellType{1});

innerRetina.compute(bp);

psth = innerRetina.mosaic{1}.get('psth');

figure; ieMovie(psth(:,:,1:10:end));

figure; sceneShowImage(scene);
