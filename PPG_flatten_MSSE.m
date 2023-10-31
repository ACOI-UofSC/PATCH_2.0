%Import PPG file
epoch=0.5;       %epoch length in seconds for removal of motion artefacts


%PPG
delimiterIn = ',';
headerlinesIn =1;   
A = importdata(filename_Patch,delimiterIn,headerlinesIn);

import_line=2;                          %row that you want to start analyzing from or STARTING POINT


PPG_raw=transpose(A.data(1+round(trial_start*fs_PPG):end,1)); %raw PPG trace in column given by second number, 3
% bit_adjust; 
PPG=PPG_raw-movmean(PPG_raw, epoch*fs_PPG);
% PPG=movmean(PPG,round(fs_PPG*0.25));   %smooth out peak systole to make fundamental dominant

%dt=(t_PPG(end)-t_PPG(1))/length(t_PPG); %sampling frequency
samples_PPG=round((epoch)*fs_PPG);   %number of samples for the epoch window
t_PPG=(1/fs_PPG)*[1:1:length(PPG)];  %set up time

%PPG

% PPG_filt=highpass(PPG,0.8,fs_PPG,'ImpulseResponse','iir','Steepness',0.5);%high pass filter to eliminate motion artefacts and other slow processes 
PPG_filt=lowpass(PPG,3.5,fs_PPG,'ImpulseResponse','iir','Steepness',0.8);%low pass filter to eliminate motion artefacts and other slow processes 
% PPG_high=highpass(PPG,3.5,fs_PPG,'ImpulseResponse','iir','Steepness',0.9);%low pass filter to eliminate motion artefacts and other slow processes

% PPG_filt=PPG;

PPG_trend=movmax(abs(PPG_filt), samples_PPG);

PPG_flat=PPG_filt./PPG_trend;  %normalize and flatten

[m,I]=findpeaks(PPG_flat.*(PPG_flat>0.5));
% plot(t_PPG, PPG_trend, t_PPG, PPG_filt)
% plot(t_PPG, PPG_flat)
% hold on
% scatter(t_PPG(I), PPG_flat(I))  %plot peak
% ylim([-2 2])

%calculate time domain HR
HR_flat=(60./[diff(t_PPG(I)),1]);
HR_flat=filloutliers(HR_flat, 'previous','movmedian',40);
% HR_flat=filloutliers(HR_flat, 'previous');
% HR_flat=filloutliers(HR_flat, 'previous');
HR_flat=movmean(HR_flat,time_resolution);  

% figure
% plot(t_PPG(I),HR_flat)

% figure
% histfit(HR_flat)
% pd=fitdist(transpose(HR_flat),'Normal')  %fit the histogram of HR with normal distribution

%find breaths
% breath=HR_flat./movmean(HR_flat, 10);
% [height, time_index]=findpeaks(breath);
% 
% figure
% plot(t_PPG(I),(HR_flat./movmean(HR_flat, 10)))
% hold on
% scatter(t_PPG(I(time_index)),height)
    


