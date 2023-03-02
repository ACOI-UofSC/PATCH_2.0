HR_spect_smooth_width_V5;
accel_V5;
close all

%correct using accelerometry
threshold_HR_accel=200; %set threshold for HR motion artefacts with accelerometry
p_x(:,L)=p_a(:,L).*(p_a(:,L)>threshold_HR_accel); 
p_x(:,L)=p_x(:,L)+1*(p_x(:,L)==0); %replace all 0's with 1's for division  

% p_x(:,L)=(p_a(:,L)<threshold_HR_accel); 
% p_x(:,L)=p_x(:,L)+100*(p_x(:,L)==0); %replace all 0's with 1's for division  


p_corr=p./p_x;

x=1:1:length(t);    

p_corr(:,x)=p_corr(:,x)./max(p_corr(:,x));

idx=x./x;
width=x./x;

[pk1, idx1]=max(p_corr(:,1));
% idx(1)=idx1*(specify_start==0)+specify_start*(specify_start>0); 
%if specify_start=0, uses extracted first HR

%run fourier peak detect and decision tree for each time.
for i=2:1:length(t)
[X, iX, wX]=findpeaks(p_corr(:,i).*(p_corr(:,i)>0.3));

[Xmin,iXmin]=min(abs(iX-idx(i-1)));
idx(i)=iX(iXmin);
%if jump is too big, keep previous value
% idx(i)=iX(iXmin)*(abs(f(idx(i))-f(idx(i-1)))<=0.8)+idx(i-1)*(abs(f(idx(i))-f(idx(i-1)))>0.8);
width(i)=wX(iXmin);
end

HR_corr=60*f(idx);  %derivation of heart rate from peak frequency

%HR spectrogram calculation, with outliers removed and replaced
HR_corr_smooth=movmean(filloutliers(HR_corr,'spline','movmedian',20),1)';
HR_corr_smooth=interp1(t, HR_corr_smooth, t_PPG(I));  %HR_sp_smooth upsample to t_PPG x-axiss

s=surf(t,f*60,p_corr, 'FaceAlpha',0.4);
s.EdgeColor = 'none';
view(0,90)
hold on
plot(t_PPG(I), HR_corr_smooth)
ylim([0 600])
xlabel('time(s)')
ylabel('HR/Motion in BPM/RPM')
%resampled HR_polar to PPG data
plot(t_PPG(I), HR_polar_interp)
% scatter(t_a, 60*step_speed) %motion artefacts
% legend('PPG Spectrogram', 'HR PPG-sp', 'HR Polar','motion artefact from accelerometry')

%Bland Atlman
figure
scatter(HR_polar_interp, (HR_corr_smooth-HR_polar_interp))

%Printing all metrics
Polar_accuracy=sum(abs(HR_polar_interp-HR_corr_smooth)<5)/length(HR_corr_smooth)
PPG_self_consistency=sum(abs(HR_flat-HR_corr_smooth)<10)/length(HR_corr_smooth)
Polar_RMSE=nanstd((HR_corr_smooth-HR_polar_interp))
Polar_MAE_BPM=nanmean(abs(HR_corr_smooth-HR_polar_interp))
Polar_MAE_percent=nanmean(abs(HR_corr_smooth-HR_polar_interp)./HR_polar_interp)
time_PPG=(length(PPG)/fs_PPG)/60  %time of PPG in minutes

mean_width=mean((10/1024)*width*60)
std_width=std((10/1024)*width*60)
median_width=median((10/1024)*width*60)
max_width=max((10/1024)*width*60)

figure
width_confidence;
