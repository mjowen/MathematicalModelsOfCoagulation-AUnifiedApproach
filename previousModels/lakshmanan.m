% TF
% VII
% TF==VII
% VIIa
% TF==VIIa
% Xa
% IIa
% X
% TF==VIIa==X
% TF==VIIa==Xa
% IX
% TF==VIIa==IX
% IXa
% II
% VIII
% VIIIa
% IXa==VIIIa
% IXa==VIIIa==X
% VIIIa1-L
% VIIIa2
% V
% Va
% Xa==Va
% Xa==Va==II
% mIIa
% TFPI
% Xa==TFPI
% TF==VIIa==Xa==TFPI
% ATIII
% Xa==ATIII
% mIIa==ATIII
% IXa==ATIII
% IIa==ATIII
% TF==VIIa==ATIII
% XI
% XI==IIa
% XIa
% XIa==IX
% XIa==ATIII
% C1-inh
% XIa==C1-inh
% IXa==X

function [t,thr,totThr,sol]=main(k,c0,maxt)
% k - Parameters
% c0 - Initial concentrations
options = odeset('AbsTol',1e-16, 'NonNegative',1);
sol = ode23tb(@(t,y)f(y,k),[0,maxt],c0,options);
t = sol.x;
y = sol.y;
thr = y(7,:);
totThr = y(7,:)+1.2*y(25,:);
end
function dy = f(y,k)
dy = zeros(size(y));

    % TF + VII < TF==VII   k1
    rate = k(2)*y(3); %TF==VII
    dy([3,1,2]) = dy([3,1,2]) + [-rate;rate;rate];

    % TF + VII > TF==VII   k2
    rate = k(1)*y(1)*y(2); %TF + VII
    dy([1,2,3]) = dy([1,2,3]) + [-rate;-rate;rate];

    % TF + VIIa < TF==VIIa   k3
    rate = k(4)*y(5); %TF==VIIa
    dy([5,1,4]) = dy([5,1,4]) + [-rate;rate;rate];

    % TF + VIIa > TF==VIIa   k4
    rate = k(3)*y(1)*y(4); %TF + VIIa
    dy([1,4,5]) = dy([1,4,5]) + [-rate;-rate;rate];

    % TF==VIIa + VII > TF==VIIa + VIIa   k5
    rate = k(5)*y(5)*y(2); %TF==VIIa + VII
    dy([2,4]) = dy([2,4]) + [-rate;rate];

    % TF==VIIa + TF==VII > TF==VIIa + TF==VIIa
    rate = k(6)*y(5)*y(3);
    dy([3,5]) = dy([3,5]) + [-rate;rate];

    % Xa + VII > Xa + VIIa   k6
    rate = k(7)*y(6)*y(2); %Xa + VII
    dy([2,4]) = dy([2,4]) + [-rate;rate];

    % Xa + TF==VII > Xa + TF==VIIa
    rate = k(8)*y(6)*y(3); %Xa + VII
    dy([3,5]) = dy([3,5]) + [-rate;rate];

    % IIa + VII > IIa  + VIIa  k7
    rate = k(9)*y(7)*y(2); %IIa + VII
    dy([2,4]) = dy([2,4]) + [-rate;rate];

    % IIa + TF==VII > IIa  + TF==VIIa
    rate = k(10)*y(7)*y(3); %IIa + VII
    dy([3,5]) = dy([3,5]) + [-rate;rate];

    % TF==VIIa + X < TF==VIIa==X   k8
    rate = k(12)*y(9); %TF==VIIa==X
    dy([9,5,8]) = dy([9,5,8]) + [-rate;rate;rate];

    % TF==VIIa + X > TF==VIIa==X   k9
    rate = k(11)*y(8)*y(5); %X + TF==VIIa
    dy([5,8,9]) = dy([5,8,9]) + [-rate;-rate;rate];

    % TF==VIIa==X > TF==VIIa==Xa  k10
    rate = k(13)*y(9); %TF==VIIa==X
    dy([9,10]) = dy([9,10]) + [-rate;rate];

    % TF==VIIa + Xa < TF==VIIa==Xa   k11
    rate = k(15)*y(10); %TF==VIIa==Xa
    dy([10,5,6]) = dy([10,5,6]) + [-rate;rate;rate];

    % TF==VIIa + Xa > TF==VIIa==Xa   k12
    rate = k(14)*y(5)*y(6); %TF==VIIa + Xa
    dy([5,6,10]) = dy([5,6,10]) + [-rate;-rate;rate];

    % TF==VIIa + IX < TF==VIIa==IX   k13
    rate = k(17)*y(12); %TF==VIIa==IX
    dy([12,5,11]) = dy([12,5,11]) + [-rate;rate;rate];

    % TF==VIIa + IX > TF==VIIa==IX   k14
    rate = k(16)*y(5)*y(11); %TF==VIIa + IX
    dy([5,11,12]) = dy([5,11,12]) + [-rate;-rate;rate];

    % TF==VIIa==IX > TF==VIIa + IXa  k15
    rate = k(18)*y(12); %TF==VIIa==IX
    dy([12,5,13]) = dy([12,5,13]) + [-rate;rate;rate];

    % Xa + II > Xa + IIa   k16
    rate = k(19)*y(6)*y(14); %Xa + II
    dy([14,7]) = dy([14,7]) + [-rate;rate];

    % IIa + VIII > IIa + VIIIa   k17
    rate = k(20)*y(7)*y(15); %IIa + VIII
    dy([15,16]) = dy([15,16]) + [-rate;rate];

    % VIIIa + IXa < IXa==VIIIa   k18
    rate = k(22)*y(17); %IXa==VIIIa
    dy([17,13,16]) = dy([17,13,16]) + [-rate;rate;rate];

    % VIIIa + IXa > IXa==VIIIa   k19
    rate = k(21)*y(16)*y(13); %VIIIa + IXa
    dy([13,16,17]) = dy([13,16,17]) + [-rate;-rate;rate];

    % IXa==VIIIa + X < IXa==VIIIa==X   k20
    rate = k(24)*y(18); %IXa==VIIIa==X
    dy([18,8,17]) = dy([18,8,17]) + [-rate;rate;rate];

    % IXa==VIIIa + X > IXa==VIIIa==X   k21
    rate = k(23)*y(8)*y(17); %X + IXa==VIIIa
    dy([8,17,18]) = dy([8,17,18]) + [-rate;-rate;rate];

    % IXa==VIIIa==X > IXa==VIIIa + Xa   k22
    rate = k(25)*y(18); %IXa==VIIIa==X
    dy([18,6,17]) = dy([18,6,17]) + [-rate;rate;rate];

    % VIIIa < VIIIa1-L + VIIIa2   k23
    rate = k(27)*y(19)*y(20); %VIIIa1-L + VIIIa2
    dy([19,20,16]) = dy([19,20,16]) + [-rate;-rate;rate];

    % VIIIa > VIIIa1-L + VIIIa2   k24
    rate = k(26)*y(16); %VIIIa
    dy([16,19,20]) = dy([16,19,20]) + [-rate;rate;rate];

    % IXa==VIIIa==X > VIIIa1-L + VIIIa2 + X + IXa   k25
    rate = k(28)*y(18); %IXa==VIIIa==X
    dy([18,8,13,19,20]) = dy([18,8,13,19,20]) + [-rate;rate;rate;rate;rate];

    % IXa==VIIIa > VIIIa1-L + VIIIa2 + IXa   k26
    rate = k(29)*y(17); %IXa==VIIIa
    dy([17,13,19,20]) = dy([17,13,19,20]) + [-rate;rate;rate;rate];

    % IIa + V > IIa + Va   k27
    rate = k(30)*y(7)*y(21); %IIa + V
    dy([21,22]) = dy([21,22]) + [-rate;rate];

    % Xa + Va < Xa==Va   k28
    rate = k(32)*y(23); %Xa==Va
    dy([23,6,22]) = dy([23,6,22]) + [-rate;rate;rate];

    % Xa + Va > Xa==Va   k29
    rate = k(31)*y(6)*y(22); %Xa + Va
    dy([6,22,23]) = dy([6,22,23]) + [-rate;-rate;rate];

    % Xa==Va + II < Xa==Va==II   k30
    rate = k(34)*y(24); %Xa==Va==II
    dy([24,14,23]) = dy([24,14,23]) + [-rate;rate;rate];

    % Xa==Va + II > Xa==Va==II   k31
    rate = k(33)*y(23)*y(14); %Xa==Va + II
    dy([14,23,24]) = dy([14,23,24]) + [-rate;-rate;rate];

    % Xa==Va==II > Xa==Va + mIIa   k32
    rate = k(35)*y(24); %Xa==Va==II
    dy([24,23,25]) = dy([24,23,25]) + [-rate;rate;rate];

    % mIIa + Xa==Va > IIa + Xa==Va   k33
    rate = k(36)*y(25)*y(23); %mIIa + Xa==Va
    dy([25,7]) = dy([25,7]) + [-rate;rate];

    % Xa + TFPI < Xa==TFPI   k34
    rate = k(38)*y(27); %Xa==TFPI
    dy([27,6,26]) = dy([27,6,26]) + [-rate;rate;rate];

    % Xa + TFPI > Xa==TFPI   k35
    rate = k(37)*y(6)*y(26); %Xa + TFPI
    dy([6,26,27]) = dy([6,26,27]) + [-rate;-rate;rate];

    % TF==VIIa==Xa + TFPI < TF==VIIa==Xa==TFPI   k36
    rate = k(40)*y(28); %TF==VIIa==Xa==TFPI
    dy([28,10,26]) = dy([28,10,26]) + [-rate;rate;rate];

    % TF==VIIa==Xa + TFPI > TF==VIIa==Xa==TFPI   k37
    rate = k(39)*y(10)*y(26); %TF==VIIa==Xa + TFPI
    dy([10,26,28]) = dy([10,26,28]) + [-rate;-rate;rate];

    % TF==VIIa + Xa==TFPI > TF==VIIa==Xa==TFPI   k38
    rate = k(41)*y(5)*y(27); %TF==VIIa + Xa==TFPI
    dy([5,27,28]) = dy([5,27,28]) + [-rate;-rate;rate];

    % Xa + ATIII > Xa==ATIII   k39
    rate = k(42)*y(29)*y(6); %ATIII + Xa
    dy([6,29,30]) = dy([6,29,30]) + [-rate;-rate;rate];

    % mIIa + ATIII > mIIa==ATIII   k40
    rate = k(43)*y(25)*y(29); %mIIa + ATIII
    dy([25,29,31]) = dy([25,29,31]) + [-rate;-rate;rate];

    % IXa + ATIII > IXa==ATIII   k41
    rate = k(44)*y(13)*y(29); %IXa + ATIII
    dy([13,29,32]) = dy([13,29,32]) + [-rate;-rate;rate];

    % IIa + ATIII > IIa==ATIII   k42
    rate = k(45)*y(7)*y(29); %IIa + ATIII
    dy([7,29,33]) = dy([7,29,33]) + [-rate;-rate;rate];
    
    % TF==VIIa + ATIII > TF==VIIa==ATIII   k43
    rate = k(46)*y(5)*y(29); %TF==VIIa + ATIII
    dy([5,29,34]) = dy([5,29,34]) + [-rate;-rate;rate];

    % XI + IIa > XI==IIa
    rate = k(47)*y(35)*y(7);
    dy([35,7,36]) = dy([35,7,36]) + [-rate;-rate;rate];

    % XI + IIa < XI==IIa
    rate = k(48)*y(36);
    dy([36,35,7]) = dy([36,35,7]) + [-rate;rate;rate];
    
    % XI==IIa > XIa + IIa
    rate = k(49)*y(36);
    dy([36,37,7]) = dy([36,37,7]) + [-rate;rate;rate];
    
    % XIa + IX > XIa==IX
    rate = k(50)*y(37)*y(11);
    dy([37,11,38]) = dy([37,11,38]) + [-rate;-rate;rate];
    
    % XIa + IX < XIa==IX
    rate = k(51)*y(38);
    dy([38,37,11]) = dy([38,37,11]) + [-rate;rate;rate];
    
    % XIa==IX > XIa + IXa
    rate = k(52)*y(38);
    dy([38,37,13]) = dy([38,37,13]) + [-rate;rate;rate];
    
    % XIa + ATIII > XIa==ATIII
    rate = k(53)*y(37)*y(29);
    dy([37,29,39]) = dy([37,29,39]) + [-rate;-rate;rate];
    
    % XIa + C1-inh > XIa==C1-inh
    rate = k(54)*y(37)*y(40);
    dy([37,40,41]) = dy([37,40,41]) + [-rate;-rate;rate];
    
    % IXa + X > IXa==X
    rate = k(55)*y(13)*y(8);
    dy([13,8,42]) = dy([13,8,42]) + [-rate;-rate;rate];
    
    % IXa + X < IXa==X
    rate = k(56)*y(42);
    dy([42,13,8]) = dy([42,13,8]) + [-rate;rate;rate];
    
    % IXa==X > IXa + Xa
    rate = k(57)*y(42);
    dy([42,13,6]) = dy([42,13,6]) + [-rate;rate;rate];
    
end
