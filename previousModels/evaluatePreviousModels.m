load('previousModelsRates.mat')
load('od.mat')
load('IC.mat')

maxT = 20*60;
hockinFromIC = @(ic)(hockin(hockinK,ic,maxT));
danforthFromIC = @(ic)(danforth(danforthK,ic,maxT));
brummelFromIC = @(ic)(brummel(brummelK,ic,maxT));
bungayFromIC = @(ic)(bungay(bungayK,ic,maxT));
tyurinFromIC = @(ic)(tyurin(tyurinK,ic,maxT));
lakshmananFromIC = @(ic)(lakshmanan(lakshmananK,ic,maxT));
chatterjeeFromIC = @(ic)(chatterjee(chatterjeeK,ic,maxT));
panteleevFromIC = @(ic)(panteleev(panteleevK,ic,maxT));

% Set the baseline initial conditions: Not needed for all the models (some,
% like Hockin, override all the baseline ICs), but is needed for some (like C1-inh in Chatterjee)
hockinBaselineIC = zeros(34,1);
% TF VII VIIa X IX II VIII V TFPI AT
hockinBaselineIC([1,2,4,8,11,14,15,21,26,29])=[15e-12,1.0e-8,1.0e-10,1.6e-7,9.0e-8,1.4e-6,7.0e-10,2.0e-8,2.5e-9,3.4e-6];

danforthBaselineIC = zeros(34,1);
% TF VII VIIa X IX II VIII V TFPI AT
danforthBaselineIC([1,2,4,8,11,14,15,21,26,29])=[15e-12,1.0e-8,1.0e-10,1.6e-7,9.0e-8,1.4e-6,7.0e-10,2.0e-8,2.5e-9,3.4e-6];

brummelBaselineIC = zeros(58,1);
% TF VII VIIa X IX II VIII V TFPI AT
brummelBaselineIC([1,2,4,8,11,14,15,21,26,29])=[15e-12,1.0e-8,1.0e-10,1.6e-7,9.0e-8,1.4e-6,7.0e-10,2.0e-8,2.5e-9,3.4e-6];

bungayBaselineIC = zeros(74,1);
% TFL IIf Vf VIIf VIIaf VIIIf IXf Xf PSf XIf ATf LIPID
bungayBaselineIC([35,1,5,9,11,13,17,21,27,52,58,73])=[15e-12,1.4e-6 2e-8 1e-8 1e-10 7e-10 9e-8 1.6e-7 0 3e-8 3.4e-6 6.79e-3];

tyurinBaselineIC = zeros(33,1);
% TF VIIa VII X V XI IX VIII II PC TFPI AT PAI1 a1AT a2AP a2M
tyurinBaselineIC([1,4,5,7,9,12,14,16,19,23,24,26,28,29,30,31])=[15e-12,1e-10 1e-8 1.6e-7 2e-8 3e-8 9e-8 7e-10 1.4e-6 0 2.5e-9 3.4e-6 4.6e-10 4e-5 9.5e-7 3.25e-6];

lakshmananBaselineIC = zeros(42,1);
% TF VII VIIa X IX II VIII V TFPI AT XI
lakshmananBaselineIC([1,2,4,8,11,14,15,21,26,29,35])=[15e-12,1.0e-8,1.0e-10,1.6e-7,9.0e-8,1.4e-6,7.0e-10,2.0e-8,2.5e-9,3.4e-6,3e-8];

chatterjeeBaselineIC = zeros(82,1);
% TF VII VIIa X IX II VIII V TFPI AT PK C1inh XI a1AT a2AP Fbg eps
chatterjeeBaselineIC([1,2,4,8,11,14,15,21,26,29,42,48,51,57,59,66,81])=[15e-12,1.0e-8,1.0e-10,1.6e-7,9.0e-8,1.4e-6,7.0e-10,2.0e-8,2.5e-9,3.4e-6,0,2.1e-6,3e-8,4e-5,9.5e-7,0,0.01];

panteleevBaselineIC = zeros(24,1);
% TF VIIa VII IX X II I VIII VIII V XI AT TFPI
panteleevBaselineIC([3,4,5,7,9,11,13,15,17,19,20,21])=[15e-12,1e-10,1e-8,9e-8,1.6e-7,1.4e-6,0,7e-10,2e-8,3e-8,3.4e-6,2.5e-9];

% Calculate pooled ETP
[t,thr] = hockinFromIC(hockinBaselineIC);
pooledOD = thr2od(t,thr);
hockinPooledETP = pooledOD(end);

[t,thr] = danforthFromIC(danforthBaselineIC);
pooledOD = thr2od(t,thr);
danforthPooledETP = pooledOD(end);

[t,thr] = brummelFromIC(brummelBaselineIC);
pooledOD = thr2od(t,thr);
brummelPooledETP = pooledOD(end);

[t,thr] = bungayFromIC(bungayBaselineIC);
pooledOD = thr2od(t,thr);
bungayPooledETP = pooledOD(end);

[t,thr] = tyurinFromIC(tyurinBaselineIC);
pooledOD = thr2od(t,thr);
tyurinPooledETP = pooledOD(end);

[t,thr] = lakshmananFromIC(lakshmananBaselineIC);
pooledOD = thr2od(t,thr);
lakshmananPooledETP = pooledOD(end);

[t,thr] = panteleevFromIC(panteleevBaselineIC);
pooledOD = thr2od(t,thr);
panteleevPooledETP = pooledOD(end);

[t,thr] = chatterjeeFromIC(chatterjeeBaselineIC);
pooledOD = thr2od(t,thr);
chatterjeePooledETP = pooledOD(end);

% Calculate OD curves for each individual and find score/cost and ETP
scores = zeros(1,8);
etp = zeros(333,8);
trueETP = zeros(333,1);
for i = 1:333
    % Get the true OD data
    trueOD = od(i,:);
    trueETP(i) = trueOD(end);
    
    %Hockin
    hockinIC = hockinBaselineIC;
    % TF
    hockinIC(1) = IC(i,1);
    % II, V, VII, VIII, IX, X, AT
    hockinIC([14,21,2,15,11,8,29]) = [IC(i,2);IC(i,3);IC(i,4);IC(i,5);IC(i,6);IC(i,7);IC(i,9)];
    % TFPI
    hockinIC(26) = IC(i,10);
    % VIIa
    hockinIC(4) = hockinIC(2)/100;

    [t,thr] = hockinFromIC(hockinIC);
    OD = thr2od(t,thr);
    scores(1) = scores(1) + cost(t,OD/hockinPooledETP*100,trueOD);
    etp(i,1) = OD(end);
    
    %Danforth
    danforthIC = danforthBaselineIC;
    % TF
    danforthIC(1) = IC(i,1);
    % II, V, VII, VIII, IX, X, AT
    danforthIC([14,21,2,15,11,8,29]) = [IC(i,2);IC(i,3);IC(i,4);IC(i,5);IC(i,6);IC(i,7);IC(i,9)];
    % TFPI
    danforthIC(26) = IC(i,10);
    % VIIa
    danforthIC(4) = danforthIC(2)/100;
    
    [t,thr] = danforthFromIC(danforthIC);
    OD = thr2od(t,thr);
    scores(2) = scores(2) + cost(t,OD/danforthPooledETP*100,trueOD);
    etp(i,2) = OD(end);
    
    %Brummel
    brummelIC = brummelBaselineIC;
    % TF
    brummelIC(1) = IC(i,1);
    % II, V, VII, VIII, IX, X, AT
    brummelIC([14,21,2,15,11,8,29]) = [IC(i,2);IC(i,3);IC(i,4);IC(i,5);IC(i,6);IC(i,7);IC(i,9)];
    % TFPI
    brummelIC(26) = IC(i,10);
    % VIIa
    brummelIC(4) = brummelIC(2)/100;
    
    [t,thr] = brummelFromIC(brummelIC);
    OD = thr2od(t,thr);
    scores(3) = scores(3) + cost(t,OD/brummelPooledETP*100,trueOD);
    etp(i,3) = OD(end);
    
    %Bungay
    bungayIC = bungayBaselineIC;
    % TF
    bungayIC(35) = IC(i,1);
    % II, V, VII, VIII, IX, X, XI, AT
    bungayIC([1,5,9,13,17,21,52,58]) = [IC(i,2);IC(i,3);IC(i,4);IC(i,5);IC(i,6);IC(i,7);IC(i,8);IC(i,9)];
    % TFPI
    bungayIC(57) = IC(i,10);
    % VIIa
    bungayIC(11) = bungayIC(9)/100;
    
    [t,thr] = bungayFromIC(bungayIC);
    OD = thr2od(t,thr);
    scores(4) = scores(4) + cost(t,OD/bungayPooledETP*100,trueOD);
    etp(i,4) = OD(end);
    
    %Tyurin
    tyurinIC = tyurinBaselineIC;
    % TF
    tyurinIC(1) = IC(i,1);
    % II, V, VII, VIII, IX, X, XI, AT
    tyurinIC([19,9,5,16,14,7,12,26]) = [IC(i,2);IC(i,3);IC(i,4);IC(i,5);IC(i,6);IC(i,7);IC(i,8);IC(i,9)];
    % TFPI
    tyurinIC(24) = IC(i,10);
    % VIIa
    tyurinIC(4) = tyurinIC(5)/100;
    
    [t,thr] = tyurinFromIC(tyurinIC);
    OD = thr2od(t,thr);
    scores(5) = scores(5) + cost(t,OD/tyurinPooledETP*100,trueOD);
    etp(i,5) = OD(end);

    %Lakshmanan
    lakshmananIC = lakshmananBaselineIC;
    % TF
    lakshmananIC(1) = IC(i,1);
    % II, V, VII, VIII, IX, X, XI, AT
    lakshmananIC([14,21,2,15,11,8,35,29]) = [IC(i,2);IC(i,3);IC(i,4);IC(i,5);IC(i,6);IC(i,7);IC(i,8);IC(i,9)];
    % TFPI
    lakshmananIC(26) = IC(i,10);
    %VIIa
    lakshmananIC(4) = lakshmananIC(2)/100;

    [t,thr] = lakshmananFromIC(lakshmananIC);
    OD = thr2od(t,thr);
    scores(6) = scores(6) + cost(t,OD/lakshmananPooledETP*100,trueOD);
    etp(i,6) = OD(end);

    %Chatterjee
    chatterjeeIC = chatterjeeBaselineIC;
    % TF
    chatterjeeIC(1) = IC(i,1);
    % II, V, VII, VIII, IX, X, XI, AT
    chatterjeeIC([14,21,2,15,11,8,51,29]) = [IC(i,2);IC(i,3);IC(i,4);IC(i,5);IC(i,6);IC(i,7);IC(i,8);IC(i,9)];
    %TFPI
    chatterjeeIC(26) = IC(i,10);
    % VIIa
    chatterjeeIC(4) = chatterjeeIC(2)/100;
    
    [t,thr] = chatterjeeFromIC(chatterjeeIC);
    OD = thr2od(t,thr);
    scores(7) = scores(7) + cost(t,OD/chatterjeePooledETP*100,trueOD);
    etp(i,7) = OD(end);
    
    %Panteleev
    panteleevIC = panteleevBaselineIC;
    % TF
    panteleevIC(3) = IC(i,1);
    % II, V, VII, VIII, IX, X, XI, AT
    panteleevIC([11,17,5,15,7,9,19,20]) = [IC(i,2);IC(i,3);IC(i,4);IC(i,5);IC(i,6);IC(i,7);IC(i,8);IC(i,9)];
    % TFPI
    panteleevIC(21) = IC(i,10);
    % VIIa
    panteleevIC(4) = panteleevIC(5)/100;
    
    [t,thr] = panteleevFromIC(panteleevIC);
    OD = thr2od(t,thr);
    scores(8) = scores(8) + cost(t,OD/panteleevPooledETP*100,trueOD);
    etp(i,8) = OD(end);
    
    disp(strcat('Completed person: ',num2str(i)))
end

scores = sqrt(scores./333);

% Produce summary statistics
names = ["Hockin", "Danforth", "Brummel", "Bungay", "Tyurin", "Lakshmanan", "Chatterjee", "Panteleev"];
for i = 1:8
    disp("============")
    disp(names(i))
    disp("Cost")
    disp(scores(i))
    disp("RMSE")
    a = fitlm(etp(:,i), trueETP, intercept=false);
    disp(a.RMSE)
    disp("R^2")
    a = fitlm(etp(:,i), trueETP);
    disp(a.Rsquared.Ordinary)
    disp("p-value")
    disp(a.Coefficients.pValue(2))
end

function od = thr2od(t,thr)
od = zeros(size(thr));
for i = 2:length(od)
    od(i) = trapz(t(1:i),thr(1:i));
end
end

function score = cost(t,od,trueOD)
timepoints = 0:30:1200;
interpOD = interp1(t,od,timepoints);
score = sum((interpOD-trueOD).^2);
end