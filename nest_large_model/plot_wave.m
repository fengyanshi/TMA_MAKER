clear all
fdir='./output/';

SLP=0.05;
Xslp = 300.0;
DEPTH_FLAT = 8.0;

eta=load([fdir 'eta_00001']);

[n,m]=size(eta);
dx=2.0;
dy=2.0;
x=[0:m-1]*dx;
y=[0:n-1]*dy;

for i=1:m
if x(i)<Xslp
dep(i)=-DEPTH_FLAT;
else
dep(i)=-DEPTH_FLAT+SLP*(x(i)-Xslp);
end
end

nfile=[1:1:99];

myVideo = VideoWriter('videoOut.mp4','MPEG-4');
myVideo.FrameRate = 10;  
myVideo.Quality = 100;
vidHeight = 576; %this is the value in which it should reproduce
vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);

colormap jet
wid=8;
len=5;
set(gcf,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[0 0 wid len]);



for num=1:length(nfile)
    
fnum=sprintf('%.5d',nfile(num));
eta=load([fdir 'eta_' fnum]);
mask=load([fdir 'mask_' fnum]);

eta(mask==0)=NaN;


clf
subplot(1,2,1)

mesh(x,y,eta)
axis([0 x(end) 0 y(end) -2 4])
view(17,16)

subplot(1,2,2)
plot(x,eta(floor(n/2),:))
hold on
plot(x,dep,'k-')
axis([0 x(end) -2 4])
grid

if num==1
ylabel(' y (m) ')
end

xlabel(' x (m) ')
%cbar=colorbar;
%set(get(cbar,'ylabel'),'String','\eta (m) ')
title(['time = ' num2str(num*2) ' sec'])

set(gcf,'Renderer','zbuffer')

pause(0.1)

F = print('-RGBImage','-r300');
J = imresize(F,[vidHeight vidWidth]);
mov(num).cdata = J;

writeVideo(myVideo,mov(num).cdata);

end
close(myVideo)



%print -djpeg eta_inlet_shoal_irr.jpg