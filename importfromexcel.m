clear; close all; clc

deltaT = 10;

filenameimport = input('Enter file name (the file has to be in the same directory): ');
filenameexport = input('Enter name of the EXPORT file: ');

try
    Tabl = readtable(filenameimport, "VariableNamingRule","preserve");
catch
    Tabl = readtable(filenameimport); % For older MATLAB versions
end

i = 1;
while ~isnan(Tabl{i,2})
    IDh(i) = Tabl{i,1};
    Tinh(i) = Tabl{i,2};
    Touth(i) = Tabl{i,3};
    cph(i) = Tabl{i,4};
    i = i+1;
end


IDc = Tabl{i+1:end,1}';
Tinc = Tabl{i+1:end,2}';
Toutc = Tabl{i+1:end,3}';
cpc = Tabl{i+1:end,4}';

Tablexport = networkdesign(IDh, Tinh, Touth, cph, IDc, Tinc, Toutc, cpc, deltaT)

writetable(Tablexport,filenameexport)