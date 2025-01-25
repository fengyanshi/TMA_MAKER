clear all
fdir='../TMA/';
fre=load([fdir 'SPC_frq.txt']);
angle=load([fdir 'SPC_angle.txt']);
spc=load([fdir 'SPC_HMO.txt']);

pcolor(fre,angle,spc),shading flat

Etotal=sum(sum(spc.^2))/8.0;
Hrms=sqrt(8*Etotal);

E=spc.^2/8.0;

dangle=diff(angle)*pi/180.;
dfre=diff(fre);
[n,m]=size(spc);

for kf=1:m-1
    for ktheta=1:n-1
        Em(ktheta,kf)=0.25*(E(ktheta,kf)...
            +E(ktheta+1,kf)+E(ktheta,kf+1)+E(ktheta+1,kf+1));
        Dens(kf,ktheta)=Em(ktheta,kf)/dangle(ktheta)/dfre(kf);
        F(kf)=0.5*(fre(kf)+fre(kf+1));
        A(ktheta)=0.5*(angle(ktheta)+angle(ktheta+1));
    end
end

E_a=dangle*0.0;
E_f=dfre*0.0;

for kf=1:m-1
for ktheta=1:n-1
    E_f(kf)=E_f(kf)+Em(ktheta,kf);
end
    D_f(kf)=E_f(kf)/dfre(kf);
end

for ktheta=1:n-1
for kf=1:m-1
    E_a(ktheta)=E_a(ktheta)+Em(ktheta,kf);
end
    D_a(ktheta)=E_a(ktheta)/dangle(ktheta);
end


Etotal_1=sum(sum(Em));


d_range=[-50 50];
f_range=[0.05 0.25];

% plot
figure(1)
clf
subplot(2,2,1)
plot(D_f,F)
xlabel('S(f)/Hz')
ylabel('Freq(Hz)')
grid
axis([0 max(D_f) f_range(1) f_range(2)])


subplot(2,2,4)
plot(A,D_a)
ylabel('S(\theta)/ Rad')
xlabel('Dire(deg)')
grid
axis([d_range(1) d_range(2) 0 max(D_a)]);

subplot(2,2,2)
v=[0:10:120];
%c=contour(Dire,Freq,Amp,v);
contourf(A,F,Dens,40,'Edge','none')
%clabel(c,v)
xlabel('Dire(deg)')
ylabel('Freq(Hz)')
grid
axis([d_range(1) d_range(2) f_range(1) f_range(2)]);
title('m^2/Hz/Rad')
%colorbar

%totalE=sum(Dens.*dangle.*dfre);
Hrms=sqrt(Etotal_1*8);

subplot(2,2,3)
text(0.2,0.5,['Hrms = ' num2str(Hrms)  ' m'])
axis off

