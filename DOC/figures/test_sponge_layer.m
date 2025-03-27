clear all
R=0.85;
width=40;
i=[1:width];
A=100.0;
ri=R.^(50*i/width);
spp=A.^ri;
spp=1./spp;

hold on
plot(spp)
