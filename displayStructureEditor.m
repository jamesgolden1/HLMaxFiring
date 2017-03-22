%pathToMatFileToAlter = '/Users/bireswarlaha/Box Sync/Neuro/github/isetbioModified/Vive-HTC-leftEye.mat';
pathToMatFileToAlter = '/Users/bireswarlaha/Box Sync/Neuro/github/isetbioModified/Vive-HTC-rightEye.mat';

inputData = load(pathToMatFileToAlter);
displayStructure = inputData.d;
disp(class(displayStructure))

names = fieldnames(displayStructure)

%fprintf(displayStructure.fieldnames());

type = displayStructure.type;
name = displayStructure.name;
wave = displayStructure.wave;
spd = displayStructure.spd;
gamma = displayStructure.gamma;
dpi = displayStructure.dpi;
dist = displayStructure.dist;
isEmissive = displayStructure.isEmissive;
dixel = displayStructure.dixel;
ambient = displayStructure.ambient;

%updating values below:
name = 'Vive-HTC-rightEye'
displayStructure.name = name;
dpi = 456 %source: https://steamcommunity.com/app/358040/discussions/0/365163686079554210/
displayStructure.dpi = dpi;
dist = 0.02 %the screens in VR headsets are like 2 cms away from the eye
displayStructure.dist = dist;

inputData.d = displayStructure;

%writing back in the output file
save '/Users/bireswarlaha/Box Sync/Neuro/github/isetbioModified/Vive-HTC-rightEye.mat' inputData;
