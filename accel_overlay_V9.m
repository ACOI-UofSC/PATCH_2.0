%CODE VERSION DATA 7_27_2021
%Incorporates two different spectrogam power thresholds for steps vs.
%motion artefact
%accelerometer processing
clear
participant_parameters_V6;
fs_accel=fs_PPG;
A = importdata(filename_Patch,',',1);  %importdata

% Import data. Comment out if already imported.
x_a=A.data(1+round(trial_start*fs_PPG):end,2);
y_a=A.data(1+round(trial_start*fs_PPG):end,3);
z_a=A.data(1+round(trial_start*fs_PPG):end,4);

%convert acceleration in bits to G's
x_a=(x_a/1024)*25.2-16;
y_a=(y_a/1024)*25.2-16;
z_a=(z_a/1024)*25.2-16;


% PPG=accelerometer.data(1+round(trial_start*fs_PPG):round(trial_length*fs_PPG),5);
index=1:1:length(x_a);
time=index./fs_accel;
mag=sqrt(x_a.^2+y_a.^2+z_a.^2);  %compute magnitude of accel 
jerk=[diff(mag);0];


% mag_filt=mag_filt/maximum_accel; %normalize to maximum value of acceleration

%import actiheart Accel as sampled to 50Hz by Actiheart
AA=importdata(filename_Actiheart_accl,'\t',6);
x_actiheart=AA.data(1+round(trial_start*50):end,1);
y_actiheart=AA.data(1+round(trial_start*50):end,2);
z_actiheart=AA.data(1+round(trial_start*50):end,3);
mag_actiheart=sqrt(x_actiheart.^2+y_actiheart.^2+z_actiheart.^2)/9.81;

%compute jerk and set up time vector
jerk_actiheart=[diff(mag_actiheart);0];
time_actiheart=(1/50)*[0:1:length(mag_actiheart)-1];
time_actiheart=time_actiheart+t_Actiheart_shift;  %offset to be varied

%truncate actiheart data to match patch time
[start_actiheart, i_start]=min(abs(time_actiheart-time(1)));  %find starting
[end_actiheart, i_end]=min(abs(time_actiheart-time(end)));      %find ending
time_actiheart=time_actiheart(i_start:i_end);
jerk_actiheart=jerk_actiheart(i_start:i_end);
% 
% yyaxis left
% plot(time, jerk)
% yyaxis right
% plot(time_actiheart, jerk_actiheart)

%downsample PATCH jerk down to 50Hz Actiheart rate
jerk_interp=interp1(time, jerk, time_actiheart)';
jerk_interp=fillmissing(jerk_interp,'spline');

%optimize time_shift by getting cross-correlation between actiheart and
%patch jerk signals max lead/lag samples is 3rd argument in xcorr
[correl, lags]=xcorr(jerk_interp, jerk_actiheart,50*300);
[xx, yy]=max(correl);
t_shift=lags(yy)/50;  %convert lag in vector index to second using 50Hz

%remake actiheart jerk based on refined delay t_shift
jerk_actiheart=[diff(mag_actiheart);0];
time_actiheart=(1/50)*[0:1:length(mag_actiheart)-1];
time_actiheart=time_actiheart+t_Actiheart_shift+t_shift;  %offset determined by cross-correlation
%truncate again to PATCH time
[start_actiheart, i_start]=min(abs(time_actiheart-time(1)));  %find starting
[end_actiheart, i_end]=min(abs(time_actiheart-time(end)));      %find ending
time_actiheart=time_actiheart(i_start:i_end);
jerk_actiheart=jerk_actiheart(i_start:i_end);
jerk_interp=interp1(time, jerk, time_actiheart)';
jerk_interp=fillmissing(jerk_interp,'spline');


%calculate scaling factor from actiheart to patch
% average_time=5;        %time to average over in seconds becasue 50Hz rate
scatter(sqrt(movmean(jerk_actiheart.^2,50*average_time)),sqrt(movmean(jerk_interp.^2,50*average_time)));
factor=sqrt(movmean(jerk_actiheart.^2,50*average_time))\sqrt(movmean(jerk_interp.^2,50*average_time));

ENMO=cumsum(jerk_interp)-movmean(cumsum(jerk_interp),50*4); %integrate and remove slow moving artefacts
ENMO_actiheart=cumsum(jerk_actiheart)-movmean(cumsum(jerk_actiheart),50*4);

plot(time_actiheart,ENMO_actiheart, time_actiheart, (1/factor)*ENMO);
figure
scatter(sqrt(movmean(ENMO_actiheart.^2,50*average_time)),(1/factor)*sqrt(movmean(ENMO.^2,50*average_time)));
hold on

%linear regression to RMS ENMO over average time to get rsq
x=sqrt(movmean(ENMO_actiheart.^2,50*average_time));
% x=fillmissing(x,'spline');
y=(1/factor)*sqrt(movmean(ENMO.^2,50*average_time));
% y=fillmissing(y,'spline');
scaling=x\y
yfit = scaling*x;   %make vector for linear regression
plot(x, yfit);          %overlay regression fit
%residuals analysis
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal
