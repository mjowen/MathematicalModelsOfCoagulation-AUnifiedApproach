function kM = prior2param(kB)
%Converts are parameter vector from the Prior/ Bayesian form (Kd/k+ and
%Km/kcat) into the mass action law parameter vector form (k+/k- and k+/k-/kcat)

%kB - Prior (Bayesian) form parameters
%kM - Mass Action Law form parameters
kM = zeros(109,1);
% Indexes to swap between the different rate parameter forms
%Kd,k+ -> k+,k-
KdinkB = [1,3,32,59,63,65];
KdinkM = [1,3,45,81,85,87];

%Km,kcat -> k+,k-,kcat
KminkB = [5:2:17,20:2:30,36,38,41,43,46:2:54,85];
KminkM = [5:3:23,27:3:42,49,52,56,59,63:3:75,107];

% k+ -> k+
kinkB = [19,34,35,40,45,58,61,62,67:84];
kinkM = [26,47,48,55,62,80,83,84,89:106];

VIIIa_break_in_kB = 56;
VIIIa_break_in_kM = 78;

% Assign all rates of the form (Kd,k+) into (k+,k-)
for i = 1:length(KdinkB)
    kM(KdinkM(i)) = kB(KdinkB(i)+1);%k+
    kM(KdinkM(i)+1) = kB(KdinkB(i))*kB(KdinkB(i)+1);%k-=Kd*k+
end

% Assign all rates of the form (Km,kcat) into (k+,k-,kcat)
for i = 1:length(KminkB)
    kPlus = 1e8; % Ensure k- > 0 for final rates, if not then increase k+
    while kM(KminkM(i)+1)<=0 %Always run atleast once as kM starts as all 0's
        kM(KminkM(i)) = kPlus; %k+ 1e8
        kM(KminkM(i)+2) = kB(KminkB(i)+1);%kcat
        kM(KminkM(i)+1) = kB(KminkB(i))*kPlus - kB(KminkB(i)+1);%k- = Km*k+ - kcat
        kPlus = kPlus*10;
    end
end

% Assign all of the extra rates that aren't changed (but do shift in index)
for i = 1:length(kinkB)
    kM(kinkM(i)) = kB(kinkB(i));
end

% FVIIIa break is more awkward (defined in reverse to the others, forwards rate is 1st order and reverse rate is 2nd order)
kM(VIIIa_break_in_kM) = kB(VIIIa_break_in_kB+1); %k_+ fowards rate (not association rate)
kM(VIIIa_break_in_kM+1) = kB(VIIIa_break_in_kB+1)/kB(VIIIa_break_in_kB); %k_- = k_+/Kd


if sum(kM<0)>0
    disp("Negative parameters given")
    disp("Indexes at")
    find(kM<0)
end


