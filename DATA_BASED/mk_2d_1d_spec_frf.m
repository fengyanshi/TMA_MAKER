clear all

% get data
spec_data1=load('../DATA/spectra_frf.txt');
spec_data=spec_data1'; % rotate

[NumFreq,NumDir]=size(spec_data);

dfreq=0.0075;
ddire=5.;
angle_start=1;
angle_end=72;

fmin=0.04;
fmax=fmin+(NumFreq-1)*dfreq;
thmin=((angle_start-1)*ddire-71.8); % rotate 0 deg
thmax=thmin+(NumDir-1)*ddire;
Freq=[fmin:dfreq:fmax];
Dire=[thmin:ddire:thmax];



Amp=spec_data(1:end,1:end);
    %dfreq=(Freq(end)-Freq(1))/(length(Freq)-1); % Hz
    %ddire=(Dire(end)-Dire(1))/(length(Dire)-1); % degree

Amp_Dir=sum(Amp)*dfreq;
Amp_Frq=sum(Amp,2)*ddire;

f_range=[fmin fmax];
d_range=[-90 90];

figure(1)
subplot(2,2,1)
plot(Amp_Frq,Freq)
xlabel('S(f)/Hz')
ylabel('Freq(Hz)')
grid
axis([0 8 f_range(1) f_range(2)])


subplot(2,2,4)
plot(Dire,Amp_Dir)
ylabel('S(\theta)/deg')
xlabel('Dire(deg)')
grid
axis([d_range(1) d_range(2) 0 0.008]);

subplot(2,2,2)
v=[0.005 0.01:0.01:0.1];
c=contour(Dire,Freq,Amp,v);
clabel(c,v)
xlabel('Dire(deg)')
ylabel('Freq(Hz)')
grid
axis([d_range(1) d_range(2) f_range(1) f_range(2)]);
title('raw data, m^2/Hz/Deg')

totalE=sum(sum(Amp*dfreq*ddire));
Hrms=sqrt(totalE*8);



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


