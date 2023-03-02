%CODE VERSION DATE 12/14/2022 Version 6
%In version 6, we now correlate acceleration signals first to determine a
%refined time shift to better align signals between actiheart and patch.
%This adds a new time shift t_shift from a cross-correlation analysis in
%accel_overlay_V6.m
%
%V6 also truncates actiheart signal to match patch length so that accuracy
%metric is correctly calculated
%
%separating participant specific modifications into separate script
%participant_parameters_V5.m
%This script compares Patch PPG with ACtiheart ECG, and aligns the two
%using participant_parameters_V5. The datasets are aligned in time  by
%interpolating down to the smaller vector, in this case HRPPG
%Additional change, PPG_flatten time domain HR_flat is modified to be tied 
%to time_resolution in partcipant_parameters_V5 


%HR_corr_accel_V5 can be used to remove motion artefacts with an external
%optimization loop or manually within the script.
%Can specify starting HR based on participant's resting HR
%Can specify time resolution for spectrogram ONLY if needed


%Make sure version numbers are same as top level script.
clear
accel_overlay_V9;

PPG_flatten_V8; %time domain PPG HR from Systolic-Systolic spacing. 
%Make sure version number is same as top level script.

%Imports slow motion artefact and DC-offset removed PPG raw signal from
%PPG_flatten of same version number as this
% PPG_conv=PPG_raw-mean(PPG_raw);
PPG_conv=PPG_flat;


%Convert PPG data into spectrogram
% p is the power of the spectrogram
% f is the corresponding frequencies to the spectral estimates p
% t is a vector of time instants corresponding to the centers of the 
% windowed segments used to compute short-time power spectrum estimates.
[p,f,t]=pspectrum(PPG_conv, fs_PPG, 'spectrogram','FrequencyLimits',[0.7 3.5],'TimeResolution',time_resolution, "Reassign", true);
% 
% accel;  %import accelerometer data 
% p=p./p_x; %divide out motion artefacts from accelerometry


%initialize time variables, indices etc. and plot initial spectrogram
%without data
x=1:1:length(t);    
p(:,x)=p(:,x)./max(p(:,x)); % Normalize to between 0 and 1
% s=surf(t,f*60,p, 'FaceAlpha',0.4);
% s.EdgeColor = 'none';
% view(0,90);
% hold on
% ylim([0 600]);

% [pks, idx]=max(p(:,x));

idx=x./x;
width=x./x;
[pk1, idx1]=max(p(:,1));
idx(1)=idx1*(specify_start==0)+specify_start*(specify_start>0); 
%if specify_start=0, uses extracted first HR

%run fourier peak detect and decision tree for each time.
for i=2:1:length(t)
[X, iX, wX]=findpeaks(p(:,i).*(p(:,i)>0.5)); 
% X is a vector with local maxima
% iX is a vector of the location of peaks
% wX is a vector of the width of each peak

[Xmin,iXmin]=min(abs(iX-idx(i-1)));
if isempty(iXmin)
    idx(i)=idx(i-1);
    width(i)=width(i-1);
else
idx(i)=iX(iXmin);
%if jump is too big, keep previous value
idx(i)=iX(iXmin)*(abs(f(idx(i))-f(idx(i-1)))<=0.9)+idx(i-1)*(abs(f(idx(i))-f(idx(i-1)))>0.9);
width(i)=wX(iXmin);
end

end

HR_sp=60*f(idx);  %derivation of heart rate from peak frequency

%HR spectrogram calculation, with outliers removed and replaced
HR_sp_smooth=movmean(filloutliers(HR_sp,'spline','movmedian',20),1)';
HR_sp_smooth=interp1(t, HR_sp_smooth, t_PPG(I));  %HR_sp_smooth upsample to t_PPG x-axiss

%taking only from frequency domain with smooth changes
HRPPG=HR_sp_smooth;

% Import and Plotting of Actiheart Data
fid=fopen(filename_Actiheart, 'rt');
ACTIHEART = textscan(fid, '%T%f%f%s%s%s','headerlines', 7,'delimiter','\t'); %from ACtiheart per second file
fclose(fid);
t_Actiheart=86400*datenum(ACTIHEART{1,1});
t_Actiheart=t_Actiheart-t_Actiheart(1);  %sets first value of time to zero to get elapsed time only
t_Actiheart=t_Actiheart+t_Actiheart_shift-trial_start+t_shift;  %aligns with PATCH PPG time from participant parameters of the same version
HR_Actiheart=ACTIHEART{1,2};  %HR second by second from imported Actiheart data
HR_Actiheart=HR_Actiheart./(HR_Actiheart>0);  %dividing by 0's to eliminate 0's

%truncate Actiheart to patch length
[start_Actiheart, I_start]=min(abs(t_Actiheart-0));  %find starting
[end_Actiheart, I_end]=min(abs(t_Actiheart-time(end)));      %find ending
t_Actiheart=t_Actiheart(I_start:I_end);
HR_Actiheart=HR_Actiheart(I_start:I_end);

%plot HR_sp_smooth from PPG and overlay with Polar
Spect_Graph = figure;
Spect_Graph_n = strcat(participant_num,"spect.png");
s=surf(t,f*60,p, 'FaceAlpha',0.4);
s.EdgeColor = 'none';
view(0,90)
hold on
plot(t_PPG(I), HR_sp_smooth)
ylim([0 250])
xlim([0 trial_length])
plot(t_Actiheart, HR_Actiheart)
legend('PPG Spectrogram', 'HR PPG-sp','HR_Actiheart_per_second')
saveas(Spect_Graph, Spect_Graph_n)

%plot self_consistency vs. time averaged over num_beats beats
figure
num_beats=60;
plot(t_PPG(I), movmean(abs(HR_flat-HR_sp_smooth)<10,num_beats))

% Intialize an array to contain the self consistency per second
scs_per_sec = zeros(length(HR_sp_smooth)-30, 3);
k = 1; % Used to iterate through scs_per_sec
width_interp=interp1(t, width, t_PPG(I)); % Interpret Width to ppg_time

% The window in which self consistency is being looked at is 30 seconds.
% Need to find the index of the first 30 seconds.
ppg_time = t_PPG(I);    % Store the time of trial into a variable
i = 1;                   % intialize end of window
j = 1;                   % intialize begining of window
while ppg_time(i) <= 30
    i = i + 1;
end
% Calculate the self consistency per second. 
% To do this I use a moving window of 30 seconds and move it by 1 second
% each calculation.
while i < length(HR_sp_smooth)
    % Select 30 seconds of HR_sp_smooth
    smooth = HR_sp_smooth(j:i);

    %  Select 30 seconds of HR_flat
    flat = HR_flat(j:i);
    % Calculate the Self Consistency
    self_consist = sum((abs(flat-smooth)<10)/length(smooth));

    % Store the self consistency
    scs_per_sec(k, 2) = self_consist;

    % Select 30 seconds of width
    std_width = std((10/1024)*width_interp(j:i)*60);

    % Store the self consistency
    scs_per_sec(k, 3) = std_width;
    % Calculate time stamp
    scs_per_sec(k, 1) = ppg_time(floor((i + j)/2));

    k = k+1;
    j = j + 1;
    i = i + 1;
end


%Bland Atlman
% figure
% scatter(HR_polar_interp, (HR_sp_smooth-HR_polar_interp))

% NICK'S ADDITION. PLEASE DOUBLE CHECK
%Resample HR_flat time domain Actiheart HR down to t_Actiheart for metrics
HR_smooth_interp=interp1(t_PPG(I), HR_sp_smooth,t_Actiheart);

% Printing all metrics
Actiheart_accuracy=sum(abs(HR_Actiheart-HR_smooth_interp)<5)/length(HR_smooth_interp);
PPG_self_consistency=sum(abs(HR_flat-HR_sp_smooth)<10)/length(HR_sp_smooth);
Actiheart_RMSE=nanstd((HR_smooth_interp-HR_Actiheart));
Actiheart_MAE_BPM=nanmean(abs(HR_smooth_interp-HR_Actiheart));
Actiheart_MAE_percent=nanmean(abs(HR_smooth_interp-HR_Actiheart)./HR_Actiheart);
time_PPG=(length(PPG)/fs_PPG)/60;  %time of PPG in minutes

% HR_Actiheart_interp = interp1(t_Actiheart, HR_Actiheart, t_PPG(i));
% 
% Actiheart_accuracy=sum(abs(HR_Actiheart_interp-HR_sp_smooth)<5)/length(HR_sp_smooth);
% PPG_self_consistency=sum(abs(HR_flat-HR_sp_smooth)<10)/length(HR_sp_smooth);
% Actiheart_RMSE=nanstd((HR_sp_smooth-HR_Actiheart_interp));
% Actiheart_MAE_BPM=nanmean(abs(HR_sp_smooth-HR_Actiheart_interp));
% Actiheart_MAE_percent=nanmean(abs(HR_sp_smooth-HR_Actiheart_interp)./HR_Actiheart_interp);
% time_PPG=(length(PPG)/fs_PPG)/60;  %time of PPG in minutes

mean_width=mean((10/1024)*width*60);
std_width=std((10/1024)*width*60);
median_width=median((10/1024)*width*60);
max_width=max((10/1024)*width*60);


% Need to output the metrics and Heart rate
metrics = table(Actiheart_accuracy, PPG_self_consistency, Actiheart_RMSE, Actiheart_MAE_BPM, Actiheart_MAE_percent, mean_width, std_width, median_width, max_width);
metrics_n = participant_num + "_metrics.xlsx";
writetable(metrics,metrics_n);

% Call activity_lables.m to get timestamped activity lables 
activity_lables;
heart_lables = strings(length(t_Actiheart), 1);
% Iterate through all of the times that correspond with a HR reading
for i = 1:1:length(t_Actiheart)
    % Iterate through each label in acti_times
    for j = 1:1:length(acti_times)
        % If the heart rate occurs during an activity label it
        if (t_Actiheart(i) >= acti_times(j,1)) && (t_Actiheart(i) <= acti_times(j,2))
            heart_lables(i) = labels(j);
        end
    end
    % If the heart rate wasn't labeled with an activity label it as a
    % transition
    if heart_lables(i) == ""
        heart_lables(i) = "Transition";
    end
end

% Check if there are flags. If there are flags add them to the data
heart_flags = strings(length(t_Actiheart), 1);
if iscell(flag_times) == 0
    % iterate through each time corresponding to a heart rate
    for i = 1:1:length(t_Actiheart)
        % Iterate through each flag time
        for j=1:1:length(flag_times)
            % If the heart rate time falls in the flag time, label it
            if (t_Actiheart(i) >= flag_times(j,1)) && (t_Actiheart(i) <= flag_times(j,2))
                heart_flags(i) = flag(j);
            end
        end
    end
end
% Output Heart Rate to a file
%HR_Output ={[{'Time'}; t_PPG(I)'], [{'HR_Patch'}; HR_sp_smooth']};
HR_Output = table(["Time", t_Actiheart.'; "Activity", heart_lables.'; "Flags", heart_flags.'; "HR_PPG" , HR_smooth_interp.'; "HR_Actiheart", HR_Actiheart.'].');

HR_Output_n = participant_num + "_PATCH_HR.csv";
writetable(HR_Output, HR_Output_n)

% Label and Flag ENMO data in similar way that HR data was done.
enmo_lables = strings(length(time_actiheart), 1);
% Iterate through all of the times that correspond with a ENMO reading
for i = 1:1:length(time_actiheart)
    % Iterate through each label in acti_times
    for j = 1:1:length(acti_times)
        % If the enmo was measured during an activity label it
        if (time_actiheart(i) >= acti_times(j,1)) && (time_actiheart(i) <= acti_times(j,2))
            enmo_lables(i) = labels(j);
        end
    end
    % If the ENMO wasn't labeled with an activity, label it as a
    % transition
    if enmo_lables(i) == ""
        enmo_lables(i) = "Transition";
    end
end
% Check if there are flags. If there are flags add them to the data
enmo_flags = strings(length(time_actiheart), 1);
if iscell(flag_times) == 0
    % iterate through each time corresponding to a heart rate
    for i = 1:1:length(time_actiheart)
        % Iterate through each flag time
        for j=1:1:length(flag_times)
            % If the heart rate time falls in the flag time, label it
            if (time_actiheart(i) >= flag_times(j,1)) && (time_actiheart(i) <= flag_times(j,2))
                enmo_flags(i) = flag(j);
            end
        end
    end
end
% Output ENMO to a file
enmo_table = table(["Time", time_actiheart; "Activity", enmo_lables.'; "Flag", enmo_flags.'; "PATCH_ENMO", ENMO.'; "ActiHeart_ENMO", ENMO_actiheart.'].');
enmo_table_name = participant_num + "_PATCH_ENMO.csv";
writetable(enmo_table, enmo_table_name)

% Output Self Consistency Per Second

var_names = ["Time", "Self Consistency", "STD Width"];
scs_output = array2table(scs_per_sec);
scs_output.Properties.VariableNames(1:3) = {'Time', 'Self Consistency', 'STD Width'};
writetable(scs_output, participant_num + "_Self_Consistency.csv")
