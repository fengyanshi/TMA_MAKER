clear all
fdir='output/';

eta=load([fdir 'eta_00024']);
tmp=load([fdir 'tmp_00024']);

figure(1)
clf
colormap jet
subplot(121)
pcolor(eta),shading flat
xlabel('x')
ylabel('y')
caxis([-0.4 0.4])
title(['numerical T=150s'])
subplot(122)
pcolor(tmp),shading flat
xlabel('x')
ylabel('y')
caxis([-0.4 0.4])
title(['linear theory T=150s'])

print('-djpeg','theory_funwave.jpg')