% Read the OD data in sheet 1, remove time data and save to .mat file
t = readtable("pramis_data.xlsx","Sheet",1,"VariableNamingRule","preserve","ReadVariableNames",1);
t.("Time (min)") = [];

od = table2array(t);

% Code was written with od indexed the other way
od = transpose(od);
save("od.mat","od")
