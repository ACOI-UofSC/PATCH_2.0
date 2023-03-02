%CODE VERSION DATA 7_27_2021
%Incorporates two different spectrogam power thresholds for steps vs.
%motion artefact
%accelerometer processing
fs_accel=fs_PPG;
participant_parameters_V5;


% Import data. Comment out if already imported.
x_a=A.data(1+round(trial_start*fs_PPG):round(trial_length*fs_PPG),2);
y_a=A.data(1+round(trial_start*fs_PPG):round(trial_length*fs_PPG),3);
z_a=A.data(1+round(trial_start*fs_PPG):round(trial_length*fs_PPG),4);

% PPG=accelerometer.data(1+round(trial_start*fs_PPG):round(trial_length*fs_PPG),5);

index=1:1:length(x_a);

time=index./fs_accel;

mag=sqrt(x_a.^2+y_a.^2+z_a.^2);  %compute magnitude of accel 

mag_filt=mag-mean(mag);

maximum_accel=max(abs(mag_filt));

% mag_filt=mag_filt/maximum_accel; %normalize to maximum value of acceleration

[p_a,f_a,t_a]=pspectrum(mag_filt,fs_PPG, 'spectrogram','FrequencyLimits',[0 10],'TimeResolution',time_resolution);

L=1:1:length(t_a);    
% p_a(:,L)=p_a(:,L)./max(p_a(:,L));

threshold_steps=10; %set motion threshold for counting steps with accelerometry
p_steps(:,L)=p_a(:,L).*(p_a(:,L)>threshold_steps); %for counting steps
%for eliminating motion artefacts for HR. May need different thresholds
%from that previously used for steps
threshold_HR_accel=10; %set threshold for HR motion artefacts with accelerometry
p_x(:,L)=p_a(:,L).*(p_a(:,L)>threshold_HR_accel); 
p_x(:,L)=p_x(:,L)+1*(p_x(:,L)==0); %replace all 0's with 1's for division     
 


step_speed=L; %initialize step_speed vector

for j=1:1:length(L)
    [m_a,I_a]=max(p_steps(:,j));
    step_speed(j)=f_a(I_a);
end

steps=cumsum(step_speed*mean(diff(t_a)));


%comment out plotting as needed
% p_a(:,L)=p_a(:,L)./max(p_a(:,L));  % normalize each epoch
% p_steps(:,L)=p_steps(:,L)./max(p_steps(:,L));  % normalize each epoch
figure
s=surf(t_a,f_a*60,p_steps, 'FaceAlpha',0.4);
% s=surf(t_a, f_a*60, p_a(:,L)./max(p_a(:,L)),  'FaceAlpha',0.4);
s.EdgeColor = 'none';
view(0,90)

hold on
plot(t_a, 60*step_speed)
xlabel('time(s)')
ylabel('motion BPM/RPM')

yyaxis right
plot(t_a, steps, 'o r')
ylabel('steps')
% plot(t_CARS+t_delay, CARS_vigor)