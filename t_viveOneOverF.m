% t_viveOneOverF

% Build the 1/f spatiotemporal noise stimulus for the HTC Vive

% Vive specs:
% 1080x1200 per eye, 2160x1200 total
% 2160/2 pixels = 110 degrees
% 9.82 pixels per degree
% 0.1019 degrees per pixel
%% Choose cell type

clear;
close all;

%% Add paths - isetbio, HLMaxFiring, RemoteDataToolbox
% addpath(genpath(isetbioRootPath))

%% Initialize movie

% Set size of movie, based on Vive specs:
% 1080x1200 per eye, 2160x1200 total


eyeLateral = 'both';
% eyeLateral = 'left';
% eyeLateral = 'right';

zeroPad = 0;
% szCols = 1080+zeroPad; szRows = 1080+zeroPad;
% szCols = 900+zeroPad; szRows = 900+zeroPad;
szCols = 1080*(5/3)+zeroPad; szRows = 1080*(5/3)+zeroPad;
disp(['szCols:' num2str(round(szCols))]);
disp(['szRows:' num2str(round(szRows))]);

timeLength = 1; % seconds

fps = 15; % frames per second
frames = 4;%timeLength*fps;


% Binocular field of view (FOV)
fovCols = 110; % horizontal fov in degrees for vive: http://doc-ok.org/?p=1414
fovRows = 110;%122;

% % movieBig = 128*uint8(ones(szRows,szCols,frames));
% movieBig = single(zeros(szRows,szCols,frames));

% 2160/2 pixels = 110 degrees
% 9.82 pixels per degree
pixelsPerDegree = szCols/fovCols;
% 0.1019 degrees per pixel
degreesPerPixel = fovCols/szCols;

%% 

for wmov = 4%16


size1 = szCols%16^2;
len1 = frames;
rsp1 = rand(size1,size1,len1);
isp1 = i*rand(size1,size1,len1);
sp1 = rsp1 + isp1;
% alpha = wmov*(ceil(wmov/2))/2 - .5;
% beta = wmov*rem(wmov,2)/2;

% cuts are uncorrelated, therefore beta > 0

alph1 = [0 1 0 1];
alph2 = zeros(4);
beta1 = [0 0 1 1];
beta2 = zeros(4);

% alph1 = [0 .33 .66 1 0 .33 .66 1 0 .33 .66 1 0 .33 .66 1]-.5; 
% % alph2 = [0 0 0 0 .33 .33 .33 .33 .66 .66 .66 .66 1 1 1 1];
% alph2 = zeros(16);
% % beta1 = alph2;
% beta1 = [0 0 0 0 .33 .33 .33 .33 .66 .66 .66 .66 1 1 1 1]-.5;
% beta2 = alph2;

% alpha = alphaset(wmov);
% beta = betaset(wmov);

% Bandpass Filter
scaler1 = ones(len1);
% scaler1(1:7) = 1; scaler1(16:24) = 0;
scaler2 = ones(size1,size1);
% for i = 1:size1
%     for j = 1:size1
%         if (sqrt((i-size1/2)^2+(j-size1/2)^2)) > 8
%             scaler2(i,j) = 0;
%         end
%     end
% end

tdcw = ones(size1,size1,len1);
for i = 1:size1
    for j = 1:size1
        for k = 1:len1
            if ((i == size1/2) && (j == size1/2))
                tdcw(i,j,k) = scaler1(k)*(1/abs(k - len1/2))^(beta1(wmov)+ beta2(wmov)*(sqrt((i-1)^2+(j-1)^2))/(sqrt(2)*(size1)));
                % tdcw(i,j,k) = 1e-3;
                % tdcw(i,j) = 1-sqrt((i-200)^2+(j-200)^2)/(200*sqrt(2));
            elseif (k == len1/2)
                tdcw(i,j,k) = scaler2(i,j)*(1/(sqrt((i-size1/2)^2+(j-size1/2)^2))^(alph1(wmov)+ alph2(wmov)*(sqrt((i-1)^2+(j-1)^2))/(sqrt(2)*(size1))));
            else
                tdcw(i,j,k) = scaler2(i,j)*(1/(sqrt((i-size1/2)^2+(j-size1/2)^2))^(alph1(wmov) + alph2(wmov)*(k-1)/len1))*scaler1(k)*(1/abs(k - len1/2))^(beta1(wmov)+ beta2(wmov)*(sqrt((i-1)^2+(j-1)^2))/(sqrt(2)*(size1)));
            end
        end
    end
end
% How should you handle the values at (size1/2, size1/2)????

tdcw(size1/2,size1/2,:) = tdcw(size1/2+1,size1/2,:);%1;%sqrt(2*(size1/2));
% figure; imagesc(log(abs(tdcw(:,:,20)))); 
% colormap gray;
texspec1(1,:) = tdcw(23,2,:);
% figure; loglog(abs(texspec1));

spec1 = sp1.*tdcw;
% figure; surf(log(spec1(:,:,3))); colormap jet;
imf1 = ifftn(ifftshift(spec1));


% testfr = 15;
% mimf1 =abs(max(max(imf1(:,:,testfr))));
% figure; imagesc((abs((imf1(size1/8:size1-size1/8,size1/8:size1-size1/8,testfr)))./mimf1)); colormap gray; 

% Next
rsize1 = size1/5;
lside1 = length(rsize1:(size1-(rsize1)));

rimf1 = imag(imf1(rsize1:(size1-rsize1),rsize1:(size1-rsize1),:));
% rimf1 = real(imf1(50:350,50:350,:));
maxrimf1 = max(max(max(abs(rimf1))));
scfact1 = (2^8-1)/(2*maxrimf1);
wmov=1;
nimf1((lside1*(ceil(wmov/2)-1) + 1):(lside1*(ceil(wmov/2)-1) + lside1), (lside1*(rem((wmov-1),2))+1):lside1*(rem((wmov-1),2))+lside1, :) = scfact1.*rimf1 + 2^7 +1;

% (i + rd2*(ceil(n1/2)-1),j + cd2*(rem((n1-1),2)))



end %wmov

movieBig = nimf1;
%% Plot sum of frames over time
moviePiece = movieBig(:,:,1);
maxMovie = max(abs(moviePiece(:)));
% medMovie = median(moviePiece(:));
clear moviePiece

movieSmall = uint8(128 + 127*movieBig/maxMovie);

% movieSmall(end,end,:) = 128; movieSmall(end-1,end,:) = 128;
crossHairSize = 10;
rowSize = 1081.0;
colSize = 1081.0;
movieSmall(rowSize/2-crossHairSize:rowSize/2+crossHairSize, colSize/2-1:colSize/2+1,:) = 120;
movieSmall(rowSize/2-1:rowSize/2+1, colSize/2-crossHairSize:colSize/2+crossHairSize,:) = 120;

figure; imagesc(sum(abs(movieSmall),3)); colormap gray; axis equal

% clear movieBig
figure; plot(RGB2XWFormat(movieSmall(301:345,301:345,:))');
xlabel('frame'); ylabel('Black <--------------------> White');
title('Check for on vs. off center');

% Normalize intensity to gray background mean/median
for i = [2:4]; sfr = movieSmall(:,:,i); mfr(i-1) = mean(sfr(:)); end;
% movieSmall2 = 128 + (-128+movieSmall) * (128/mean(mfr(:)));
movieSmall = movieSmall * (128/mean(mfr(:)));

% Measure norm
% for i = 1:4; sfr = movieSmall(:,:,i); nfr(i) = norm(single(sfr(:))); end; figure; plot(nfr);
% sqrt(1080*1080*128^2) % comes out to about this
%% Show movie and save
disp('creating movie now...');
% p.save = false;% 
p.save = true;
% p.vname = ['C:/Users/laha/Documents/GitHub/HLMaxFiring/april18_' cellType '_fps' num2str(fps) '.avi']
% p.vname = ['C:\Users\laha\Documents\GitHub\regenInVR\media\test_fps' num2str(fps) '.avi'];
% p.vname = ['C:\Users\laha\Documents\GitHub\regenInVR\media\test2_oneOverF_April20_fps' num2str(fps) '.avi'];
p.vname = ['C:\Users\laha\Documents\GitHub\regenInVR\media\whiteNoise1.avi'];
% p.vname = ['/Users/james/Documents/matlab/isetbio/local/test2_oneOverF_May4_fps' num2str(fps) '.avi'];
p.FrameRate = fps;
% figure; 
% set(gcf,'position',[1000         157        1411        1181]);
% set(gcf,'position',[0 0   szRows   szCols]);

%disp('test1')
%ieMovie(movieBig(:,:,1:100), p);
%ieMovie(movieBig, p);
%disp('test2')
 
% ieMovie(movieBig(:,:,1:100));
% ieMovie(movieSmall(:,:,1:100));

movieSmall(end,end,:) = 0; movieSmall(end-1,end,:) = 255;


if p.save
    vObj = VideoWriter(p.vname);%,'Uncompressed AVI');
    vObj.FrameRate = p.FrameRate;
    vObj.Quality = 100;
    open(vObj);
% end

for fr = 1:size(movieSmall,3)
%     imagesc(movieBig(:,:,fr)); colormap gray    
%     imagesc(movieSmall(:,:,fr)); colormap gray; 
%     axis image; set(gca,'xticklabel','','yticklabel','');
    if p.save  
%         F = getframe(gca);%,[0 0 szRows szCols]); 
        F = movieSmall(:,:,fr);
        writeVideo(vObj,F); 
    end
%     drawnow;
end

% Write the video object if save is true
% if p.save
    close(vObj);
end

disp(sprintf('done! movie at\n%s',p.vname));