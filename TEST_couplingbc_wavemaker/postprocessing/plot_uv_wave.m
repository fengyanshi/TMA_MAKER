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
wid=8;
len=10;
set(fig,'units','inches','paperunits','inches','papersize', [wid len],'position',[1 1 wid len],'paperposition',[1 1 wid len]);

ax=[0 400 0 400];

files=[30:30];

for num=1:length(files)

fnum=sprintf('%.5d',files(num));

u=load([fdir 'u_' fnum]);
v=load([fdir 'v_' fnum]);
eta=load([fdir 'eta_' fnum]);
mask=load([fdir 'mask_' fnum]);
dep1=dep;
eta(mask<1)=NaN;
dep1(mask<1)=NaN;

%u(mask<1)=NaN;
%v(mask<1)=NaN;

uu=sqrt(u.^2+v.^2);
[vort vort1]=curl(xx,yy,u,v);

time=num2str(files(num)*100,'%0.1f');

clf

colormap(jet)

subplot(6,1,[1])

plot(xx(1,:),-dep(1,:),'k-','LineWidth',1)
hold on
plot([-1000 1000],[0 0],'k--')
grid
axis([0 ax(2) -10 1])
xlabel('x (m)')
ylabel('z (m)')
title(['time = ' time ' s'])

subplot(6,1,[2])

plot(xx(1,:),v(floor(n/2),:),'k-','LineWidth',1)
hold on
axis([0 ax(2) 0 1.2])
grid
xlabel('x (m)')
ylabel('v at middle row(m/s)')

subplot(6,1,[3 4])
pcolor(xx,yy,eta),shading flat
%caxis([-1 1])

h_bar=colorbar('location','EastOutside');
set(get(h_bar,'xlabel'),'string','wave surface (m)' )

hold on
sc=10;
sk=5;
quiver(xx(1:sk:end,1:sk:end),yy(1:sk:end,1:sk:end),u(1:sk:end,1:sk:end)*sc,v(1:sk:end,1:sk:end)*sc,0,'Color','k')

xlabel('x (m)')
ylabel('y (m)')

subplot(6,1,[5 6])
% 
hp=pcolor(xx,yy,vort);shading interp
%caxis([0.0 1.0])
hold on
%contour(xx,yy,dep,[-5:0.2:0],'Color','k')
quiver(xx(1:sk:end,1:sk:end),yy(1:sk:end,1:sk:end),u(1:sk:end,1:sk:end)*sc,v(1:sk:end,1:sk:end)*sc,0,'Color','k')

h_bar=colorbar('location','EastOutside');
set(get(h_bar,'xlabel'),'string','vorticity (1/s)' )

xlabel('x (m)')
ylabel('y (m)')


axis(ax)


end
print -djpeg case_slope_wave_15deg.jpg






