t = readtable("pramis_data.xlsx","Sheet",2,"VariableNamingRule","preserve");

% Handle missing data
t.("Fibrinogen (g/L)")(isnan(t.("Fibrinogen (g/L)"))) = mean(t.("Fibrinogen (g/L)")(~isnan(t.("Fibrinogen (g/L)"))));
t.("TFPI (% of mean control)")(isnan(t.("TFPI (% of mean control)"))) = 1.0;
% Other factors are percentage of pooled plasma
for i=4:11
    t{isnan(t{:,i}),i} = 100;
end

IC = zeros(333,11);
% TF, II, V, VII, VIII, IX, X, XI, AT, TFPI, Fbg
IC(:,1) = t.("TF (pM)") * 1e-12;
IC(:,2) = t.("II (% of pooled plasma)")/100 * 1.4e-6;
IC(:,3) = t.("V (% of pooled plasma)")/100 * 2e-8;
IC(:,4) = t.("VII (% of pooled plasma)")/100 * 1e-8;
IC(:,5) = t.("VIII (% of pooled plasma)")/100 * 7e-10;
IC(:,6) = t.("IX (% of pooled plasma)")/100 * 9e-8;
IC(:,7) = t.("X (% of pooled plasma)")/100 * 1.6e-7;
IC(:,8) = t.("XI (% of pooled plasma)")/100 * 3e-8;
IC(:,9) = t.("AT (% of pooled plasma)")/100 * 3.4e-6;
IC(:,10) = t.("TFPI (% of mean control)")/100 * 2.5e-9; % ~1.09 is the mean activity of the controls, used to normalise to a percentage
IC(:,11) = t.("Fibrinogen (g/L)") / 340000; % g/L to mol/L through 340kDa weight

save('IC.mat', 'IC')
