clear all
fdir='../output/';

dep=load([fdir 'dep.out']);
dx=2.0;
dy=2.0;
[n,m]=size(dep);
x=[0:m-1]*dx;
y=[0:n-1]*dy;
[xx,yy]=meshgrid(x,y);


fig=figure(1);
wid=14;
len=8;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[1 1 wid len]);

ax=[0 600 0 800];

files=[1:3];

% define movie file and parameters
myVideo = VideoWriter('uvmean.mp4','MPEG-4');
myVideo.FrameRate = 10;  
myVideo.Quality = 100;
%vidHeight = 576; %this is the value in which it should reproduce
%vidWidth = 1024; %this is the value in which it should reproduce
open(myVideo);

for num=1:length(files)

fnum=sprintf('%.5d',files(num));

u=load([fdir 'umean_' fnum]);
v=load([fdir 'vmean_' fnum]);
eta=load([fdir 'etamean_' fnum]);
mask=load([fdir 'mask_' fnum]);
dep1=dep;
eta(eta+dep<0.01)=NaN;
dep1(mask<1)=NaN;

fnum1=sprintf('%.5d',files(num)*100+99);
etam=load([fdir 'eta_map_' fnum1]);

uu=sqrt(u.^2+v.^2);
[vort vort1]=curl(xx,yy,u,v);

time=num2str(files(num)*100+99,'%0.1f');

clf

colormap(jet)

subplot(1,2,1)

pcolor(x,y,eta),shading flat
caxis([-0.5 2.0])
xlabel('x (m)')
ylabel('y (m)')
title(['wave average, time = ' time ' s'])
colorbar
hold on
sk=10;
sc=20;
quiver(xx(1:sk:end,1:sk:end),yy(1:sk:end,1:sk:end),u(1:sk:end,1:sk:end)*sc,v(1:sk:end,1:sk:end)*sc,0,'Color','k')

axis(ax)

subplot(1,2,2)

pcolor(x,y,etam),shading flat
caxis([-0.5 2.0])
xlabel('x (m)')
ylabel('y (m)')
title(['Large-scale data,  time = ' time ' s'])
colorbar

axis(ax)

pause(0.1)


% save image
F = print('-RGBImage','-r300');
%J = imresize(F,[vidHeight vidWidth]);
mov(num).cdata = F;

writeVideo(myVideo,mov(num).cdata);

end
close(myVideo)






