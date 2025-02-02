function printside(FIN, DIR, sidevar, t)
% Organizes the format of data printing into coupling.txt

fprintf(FIN,['\n' DIR ' SIDE']);
if isempty(sidevar)==0
for kk = 1:3
    fprintf(FIN,'\n%16.6E%16.6E%16.6E%16.6E%16.6E',sidevar(:,kk,t));
end
end
end