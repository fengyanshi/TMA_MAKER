clear all
close all

% input
period=3600;
data_time=[0:100:7200];

% output file name
coupling_filename='../coupling.txt';
logfile='log.txt';      

npoints(1)=-1; % east
npoints(2)=-1; % west
npoints(3)=200; % south
npoints(4)=200; % north

% initial
eastfine=[];
westfine=[];
southfine=[];
northfine=[];

% construct data, i.g., eastfine(npoints,var,ntime)

for ti=1:length(data_time);

for ipoint=1:npoints(1)
% u
eastfine(ipoint,1,ti)=0.0; %;
% v
eastfine(ipoint,2,ti)=0;
% eta
eastfine(ipoint,3,ti)=0.0;
end % ip

for ipoint=1:npoints(2)
westfine(ipoint,1,ti)=0.0; %;
westfine(ipoint,2,ti)=0;
westfine(ipoint,3,ti)=0.0;
end % ip

for ipoint=1:npoints(3)
southfine(ipoint,1,ti)=0.0;
southfine(ipoint,2,ti)=1.0*(npoints(3)-ipoint)/(npoints(3)-1);
southfine(ipoint,3,ti)=0.0; 
end % ip

for ipoint=1:npoints(4)
northfine(ipoint,1,ti)=0.0; 
northfine(ipoint,2,ti)=1.0*(npoints(4)-ipoint)/(npoints(4)-1);
northfine(ipoint,3,ti)=0.0; 
end % ip

end % itime

filename=coupling_filename;  

FIN = fopen(filename,'w'); 

% log file
Finfo = fopen(logfile,'w');  

fprintf(Finfo,'coupling data\nboundary info: num of points, start point');
fprintf(Finfo,'\nEAST\n\t%d\t\t%d',size(eastfine,1),1);
fprintf(Finfo,'\nWEST\n\t%d\t\t%d',size(westfine,1),1);
fprintf(Finfo,'\nSOUTH\n\t%d\t\t%d',size(southfine,1),1);
fprintf(Finfo,'\nNORTH\n\t%d\t\t%d',size(northfine,1),1);
% end log file

fprintf(FIN,'coupling data\nboundary info: num of points, start point');
fprintf(FIN,'\nEAST\n\t%d\t\t%d',size(eastfine,1),1);
fprintf(FIN,'\nWEST\n\t%d\t\t%d',size(westfine,1),1);
fprintf(FIN,'\nSOUTH\n\t%d\t\t%d',size(southfine,1),1);
fprintf(FIN,'\nNORTH\n\t%d\t\t%d',size(northfine,1),1);
fprintf(FIN,'\nTIME SERIES');
for t = 1:length(data_time)
    disp(sprintf('Writing Time Step No. %d    of   %d',t,length(data_time) ))
    fprintf(FIN,'\n\t%f',data_time(t));
    printside(FIN,'EAST',eastfine,t)
    printside(FIN,'WEST',westfine,t)
    printside(FIN,'SOUTH',southfine,t)
    printside(FIN,'NORTH',northfine,t)
end
fclose(FIN);
disp('Finished!')


fclose(Finfo);



