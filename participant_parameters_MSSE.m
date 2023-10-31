            %SPECIFY save this file in the same folder as script
participant_num = '1106';
filename_Patch=strcat(participant_num,'_TEST.TXT');    %SPECIFY save this file in the same folder as script
filename_Actiheart= strcat(participant_num,'_hr.txt');
filename_Actiheart_accl=strcat(participant_num,'_accel.txt');

trial_start=9*60+10;                  %SPECIFY trial start time in seconds, minimum is 1s, can be changed
trial_length=43*60+35;             %SPECIFY trial length time in seconds 22mins total

fs_PPG=86.8;            %PPG and accel sampling rate. SHOULD NOT CHANGE

%specify starting HR for resting HR. Ideally specify someting.
%If HR_start=0, will use calculated starting peak.

HR_start=80;    %SPECIFY starting HR in BPM. Use 0 to leave blank 
specify_start=round((HR_start/600)*1024); %convert to index for f vector in spectrogram

t_Actiheart_shift=(2*60+50);     %offset in seconds between PATCH and ACtiheart from protocol tracking sheet SPECIFY EVERYTIME

%specify time resolution resolution for spectrogram.
%Suggest default value of 30s. Should not change unless required. 
%Speak to PI about this if needed.
time_resolution=20;  %SPECIFY ONLY IF NEEDED. ASK PI.
average_time=5;      %time to average rms acceleration over in seconds SPECIFY. ASK PI

