function [ F ] = activityFeatureExtractionGyro( sX )
%featureExtraction Summary of this function goes here if I ever feel like
%writing one up.
% sX = samples X input
tH1 = 51;
lowTh = 92;
highTh = 175;
sX = sX(:);
[Max,Imax] = findpeaks(diff(sX), 'SORTSTR','descend','NPEAKS',1);
[Min,Imin] = findpeaks(-diff(sX),'SORTSTR','descend','NPEAKS',1);
Amplitude = abs(Max-Min);
Velocity = Amplitude/(Imax - Imin);
Mean = (1/size(sX,2))*sum(sX);
T_stdv = std(sX);
T_max = max(sX);
T_min = min(sX);
T_Integrate = trapz(sX);
% peaks = [];
T_average_peak_distance = 0;
% T_findpeaks_distX=[];
[peaks, loc] = findpeaks(sX, 'MinPeakHeight', tH1);
if isempty(peaks)
    T_count_findpeaks = 0;
    T_findpeaks_distX = 0;
else
    T_count_findpeaks = length(peaks);
    lp = zeros(length(peaks),1);
    if length(peaks)>1
        for i = 2:length(peaks)
            lp(i-1) = loc(i)-loc(i-1);
        end
        T_average_peak_distance = mean(lp); 
        T_findpeaks_distX = loc(end) - loc(1); %TODO: TAKE AVG, NOT MAX-MIN
    else
        T_findpeaks_distX = 0;
    end
end
%%? Threshold Stuff:?
T_countmin_1 = sum(sX<(highTh) & sX>(lowTh));
T_countmin_2 = FcountMin(sX,lowTh);
T_countmax = FcountMax(sX,highTh);

F = horzcat(Amplitude, Velocity, Mean, T_stdv, T_max, T_min, T_Integrate, T_count_findpeaks, T_findpeaks_distX, T_average_peak_distance, T_countmax, T_countmin_1, T_countmin_2);
% hold off;
end

