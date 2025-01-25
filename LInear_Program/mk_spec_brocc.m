clear all

%case1='rand_phase';
case1='zero_phase';

% ------- get data ----
spec_data1=load('brocch_data.txt');  % amplitude of each cell
spec_data=spec_data1'; % rotate

freq=load('brocch_freq.txt');
dir=load('brocch_dir.txt');


[numfreq,numdir]=size(spec_data);
% -------- check errors
if numfreq ~= length(freq)
disp('freq size does not match spec size, stop!')
return
end
if numdir ~= length(dir)
disp('dir size does not match spec size, stop!')
return
end


% PHASES

if case1(1:3) == 'zeo' 
phase=zeros(size(spec_data));
elseif case1(1:3) == 'ran'
phase=rand(numfreq,numdir)*2*pi;
else
phase=ones(size(spec_data));
end

Dep_Ser = 10.0;
g=9.81;

m=400; % grid
n=400;
dx=3.0;
dy=3.0;
TIME=[0:100]; 

% ------ energy

Amp=spec_data;
En=0.5*Amp.^2; % energy in each cell

En_Dir=sum(En); % dir distribution
En_Freq=sum(En,2); % freq distribution


% --------- plots

fmin=min(freq);
fmax=max(freq);
f_range=[fmin fmax];
d_range=[-50 50];

figure(1)
subplot(2,2,1)
plot(En_Freq,freq)
xlabel('S(f)/Hz')
ylabel('Freq(Hz)')
grid
axis([0 0.005 f_range(1) f_range(2)])


subplot(2,2,4)
plot(dir,En_Dir)
ylabel('S(\theta)/deg')
xlabel('Dire(deg)')
grid
axis([d_range(1) d_range(2) 0 0.002]);

subplot(2,2,2)
v=[0.005 0.01:0.01:0.1];
%c=contour(dir,freq,En,v);
%clabel(c,v)
contourf(dir,freq,En)
xlabel('Dire(deg)')
ylabel('Freq(Hz)')
grid
axis([d_range(1) d_range(2) f_range(1) f_range(2)]);
title('data, m^2/Hz/Deg')

fname=[case1 '.jpg']
print('-djpeg',fname)

totalE=sum(sum(En));
Hrms=sqrt(totalE*8);

% ------------ wave surface

fname=[case1 '.mp4'];
% define movie file and parameters
myVideo = VideoWriter(fname,'MPEG-4');
myVideo.FrameRate = 10;  
myVideo.Quality = 100;
%vidHeight = 576; %this is the value in which it should reproduce
%vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);

fig=figure(2)
colormap jet
wid=4;
len=7;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);

Sigma_Ser=2.0*pi*freq;

for k=1:length(TIME)

BB=cos(Sigma_Ser*TIME(k)); 
CC=sin(Sigma_Ser*TIME(k));

x_maker=[0:m-1]*dx;
y_maker=[0:n-1]*dy;
[X_maker Y_maker]=meshgrid(x_maker,y_maker);

Ein2D=zeros([n,m]);

for kf=1:length(freq)

Cm_eta=zeros([n,m]);
Sm_eta=zeros([n,m]);

for kd=1:length(dir)
WaveNumber = wvnum_omvec(Dep_Ser,freq(kf)*2*pi,g);
WaveNumberX = WaveNumber*cos(dir(kd)*pi/180.0);
WaveNumberY = WaveNumber*sin(dir(kd)*pi/180.0);
Cm_eta = Cm_eta +Amp(kf,kd)*cos(WaveNumberX*X_maker + WaveNumberY*Y_maker+phase(kf,kd));
Sm_eta = Sm_eta +Amp(kf,kd)*sin(WaveNumberX*X_maker + WaveNumberY*Y_maker+phase(kf,kd));
end  % kd

Ein2D=Ein2D+Cm_eta*BB(kf)+Sm_eta*CC(kf);

end  % kf

clf
pcolor(X_maker,Y_maker,Ein2D),shading flat
colorbar
caxis([-1 1])
title(['Time = ' num2str(TIME(k)) ' sec'])
xlabel('x(m)')
ylabel('y(m)')
pause(0.1)

% save image
F = print('-RGBImage','-r300');
%J = imresize(F,[vidHeight vidWidth]);
mov(k).cdata = F;

writeVideo(myVideo,mov(k).cdata);

end % time
close(myVideo)

return

% truncate
[v1,ind1]=max(max(Amp));
[peakf,peakd]=find(Amp==v1);
MaxFreqNum = 50;
MaxDireNum = 30;
NFreq = min(NumFreq,MaxFreqNum);
NDir  = min(NumDir,MaxDireNum);
n_peakf=[max(1,peakf-floor(NFreq/2)):min(NumFreq,peakf+floor(NFreq/2))];
n_peakd=[max(1,peakd-floor(NDir/2)):min(NumDir,peakd+floor(NDir/2))];

Freq_model=Freq(n_peakf);
Dire_model=Dire(n_peakd);
Amp_model=Amp(n_peakf,n_peakd);

PeakFreq = Freq(peakf)

Amp_Dir_model=sum(Amp_model)*dfreq;
Amp_Frq_model=sum(Amp_model,2)*ddire;


figure(2)



subplot(2,2,1)
plot(Amp_Frq_model,Freq_model)
ylabel('S(f)/Hz')
xlabel('Freq(Hz)')
grid
axis([0 8 f_range(1) f_range(2)])


subplot(2,2,4)
plot(Dire_model,Amp_Dir_model)
ylabel('S(\theta)/deg')
xlabel('Dire(deg)')
grid
axis([d_range(1) d_range(2) 0 0.008]);

subplot(2,2,2)
v=[0.005 0.01:0.01:0.1];
c=contour(Dire_model,Freq_model,Amp_model,v);
clabel(c,v)
xlabel('Dire(deg)')
ylabel('Freq(Hz)')
grid
axis([d_range(1) d_range(2) f_range(1) f_range(2)]);
title('truncated spec, m^2/Hz/Deg')

totalE_model=sum(sum(Amp_model*dfreq*ddire));
Hrms_model=sqrt(totalE_model*8);

% write out

Amp_input=sqrt(Amp_model*dfreq*ddire*8.0)/2.0;
Amp_input=Amp_input';
Phase_input=Amp_input*0.0;

PeakPeriod = 1.0/PeakFreq;

% write 2D
fname='wave2d_frf.txt';

% write data
fid=fopen(fname,'w');
fprintf(fid,'%5i %5i   - NumFreq NumDir \n',length(Freq_model),length(Dire_model));
fprintf(fid,'%10.3f   - PeakPeriod  \n',PeakPeriod);
fprintf(fid,'%10.3f   - Freq \n',Freq_model');
fprintf(fid,'%10.3f   - Dire \n',Dire_model');
dlmwrite(fname,Amp_input,'delimiter','\t','-append','precision',5);
dlmwrite(fname,Phase_input,'delimiter','\t','-append','precision',5);

fclose(fid)

% write 1D

fname='wave1d_frf.txt';

amp_frq_bin=sqrt(Amp_Frq_model*dfreq*8.0)/2.0;

fid=fopen(fname,'w');
fprintf(fid,'%5i %5i   - NumFreq NumDir \n',length(Freq_model),1);
fprintf(fid,'%10.3f   - PeakPeriod  \n',PeakPeriod);
fprintf(fid,'%10.3f   - Freq \n',Freq_model');
fprintf(fid,'%10.3f   - Dire \n',0.0');
dlmwrite(fname,amp_frq_bin,'delimiter','\t','-append','precision',5);

fclose(fid)


