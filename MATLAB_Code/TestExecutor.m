clear;close all;clc
CSVData0 = csvread('BaselineData_8.5_min.csv');
CSVData1 = csvread('WalkingData_9min.csv');
CSVData2 = csvread('JoggingRunningData_3min.csv');
CSVData = [CSVData0; CSVData1; CSVData2];
rFB = 0; %Remove From Beginning
rFE = 0; %Remove From End
a1 = CSVData(1+rFB:end-rFE,1)./9.806; a2 = CSVData(1+rFB:end-rFE,2)./9.806; a3 = CSVData(1+rFB:end-rFE,3)./9.806;
g1 = CSVData(1+rFB:end-rFE,4).*57.2958; g2 = CSVData(1+rFB:end-rFE,5).*57.2958; g3 = CSVData(1+rFB:end-rFE,6).*57.2958;
ARES = (resultant([a1,a2,a3])-1);
GRES = resultant([g1,g2,g3]);
figure(1);hold on;
tH1 = 0.40;
plot(abs(resultant([a1,a2,a3])-1)),xlim([0 length(a1)]);
title('Accelerometer Acceleration Data'); ylabel('Acceleration (g)'); xlabel('Index (16.6Hz)');
title('Accelerometer Acceleration Training Data'); ylabel('Acceleration (g)'); xlabel('Index (16.6Hz)');
rl = refline(0,tH1); rl.Color = 'r';
refline(0,1.2);
rl2 = refline(0,1.2); rl2.Color = 'm';
%% Gyr
tH2 = 63;
figure(2);hold on;
plot(resultant([g1,g2,g3]));ylabel('Rotational Movement (Degrees per Second (?/s))'); xlabel('Index (16.6Hz)');
title('Gyroscope Rotational Movement Data'); xlim([0,length(g1)]);
rl2 = refline(0,tH2); rl2.Color = 'r';
refline(0,175);
refline(0,88);
title('Gyroscope Rotational Movement Training Data'); xlim([0,length(g1)]);
rl3 = refline(0,tH2); rl3.Color = 'r';
rl4 = refline(0,175); rl4.Color = 'k';
rl5 = refline(0,88); rl5.Color = 'm';
%% TEST CLASSIFICATION! %{ 
wL = 50;
for w=1:floor(length(a1)/50)
    start = 1 + (w-1)*wL;
    wend = start+wL-1;
    Y(w) = activityClassifier(CSVData(start:wend,1),CSVData(start:wend,2),CSVData(start:wend,3),CSVData(start:wend,4),CSVData(start:wend,5),CSVData(start:wend,6));
end

%}
%% Run Feature Extraction on Cleaned Data: 
%{
X = [a1 a2 a3];
G = [g1 g2 g3];
wL = 50; %Win length ; 50 ~= 3s
for w=1:floor(length(a1)/50)
    start = 1 + (w-1)*wL;
    wend = start+wL-1;
    FX1 = activityFeatureExtractionAcc(a1(start:wend));
    FX2 = activityFeatureExtractionAcc(a2(start:wend));
    FX3 = activityFeatureExtractionAcc(a3(start:wend));
    FXR = activityFeatureExtractionAcc(ARES(start:wend));
    FG1 = activityFeatureExtractionGyro(g1(start:wend));
    FG2 = activityFeatureExtractionGyro(g2(start:wend));
    FG3 = activityFeatureExtractionGyro(g3(start:wend));
    FGR = activityFeatureExtractionGyro(GRES(start:wend));
    tXA(w,:) = [FX1 FX2 FX3 FXR FG1 FG2 FG3 FGR];
end
commandwindow;
tYClass = input('Class?\n'); 
%{ 
CLASS 1: IDLE/NOTHING/TRANSITIONS
CLASS 2: WALKING [SLOW-FAST] NO JOGGING {i.e. feet always on ground}
CLASS 3: JOGGING - RUNNING

%}
for i=1:size(tXA,1)
    tYA(i,1) = tYClass;
end
clearvars -except tXA tYA
%% Train classifier %{
filename = 'dataTraining.mat';
tX = [];
tY = [];
if exist(filename,'file')==2
    load(filename);
end
tX = [tX;tXA];
tY = [tY;tYA];
save(filename,'tX','tY');
%}
%% plot 
%{
close all;
numfeatures = 12;
for i = 1:numfeatures
    figure(i); hold on;
    plot(tX(:,i));
    plot(tX(:,i+12));
    plot(tX(:,i+24));
    plot(tX(:,i+36));
end
% tXtY = [tX,tY];
%}