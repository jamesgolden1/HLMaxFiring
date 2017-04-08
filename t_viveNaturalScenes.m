% t_viveNaturalScenes
% 
% Build the natural scene stimulus for one eye of the HTC Vive

% Vive specs:
% 1080x1200 per eye, 2160x1200 total
% 2160 pixels = 110 degrees
% 19.6364 pixels per degree
% 0.0509 degrees per pixel

% For midgets, 1 pixel per degree under 10 degrees ecc

%% Add paths - isetbio, HLMaxFiring, RemoteDataToolbox
% addpath(genpath(isetbioRootPath))

%% Initialize movie

clear

% Set size of movie, based on Vive specs:
% 1080x1200 per eye, 2160x1200 total

zeroPad = 0;
szCols = 1080+zeroPad; szRows = 1200+zeroPad;
disp(['szCols:' num2str(round(szCols))]);
disp(['szRows:' num2str(round(szRows))]);

timeLength = 20; % seconds

fps = 30; % frames per second
frames = timeLength*fps;

% Binocular field of view (FOV)
fovCols = 110; % horizontal fov in degrees for vive: http://doc-ok.org/?p=1414
fovRows = 122;

% 2160/2 pixels = 110 degrees
% 9.82 pixels per degree
pixelsPerDegree = szCols/fovCols;
% 0.1019 degrees per pixel
degreesPerPixel = fovCols/szCols;

%% Add to big movie

% Initialize movie
% load(['/Users/james/Documents/MATLAB/'...
%     'akheitman/NSEM_mapPRJ/Stimuli/'...
%     'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);
load(['dat/testmovie_schemeA_8pix_Identity_8pix.mat']);
movieSmall = testmovie.matrix(:,:,1:frames);
movieSmall = permute(movieSmall,[2 1 3]);

eccind = 0;

% Set center - Vive display is 55 degrees per eye
degCenterX = fovRows/2; degCenterY = fovCols/2;

figure; imagesc(-szRows/2+1:szRows/2,-szCols/2+1:szCols/2,sum(movieSmall,3)); colormap gray; 
axis(2*[-100 100 -100 100]); axis equal

figure; imagesc(sum(movieSmall,3)); colormap gray; axis equal

% clear movieBig

%% Show movie and save
disp('creating movie now...');
% p.save = false;% 
p.save = true;
% p.vname = ['/Users/james/Documents/matlab/HLMaxFiring/testNS_April7_fps' num2str(fps) '.avi']
p.vname = ['C:\Users\laha\Documents\GitHub\regenInVR\media\testNS_April7_fps' num2str(fps) '.avi'];
% p.vname = 'test_April7.avi';
p.FrameRate = fps;
% figure; 
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
% end

for fr = 1:size(movieSmall,3)
%     imagesc(movieBig(:,:,fr)); colormap gray    
%     imagesc(movieSmall(:,:,fr)); colormap gray; axis image; set(gca,'xticklabel','','yticklabel','');
    if p.save  
%         F = getframe;%(gca);%,[0 0 szRows szRows]); 
        
        F = movieSmall(:,:,fr);
        writeVideo(vObj,F); 
    end
%     drawnow;
end

% Write the video object if save is true
% if p.save
%     writeVideo(vObj,F);
    close(vObj);
end
