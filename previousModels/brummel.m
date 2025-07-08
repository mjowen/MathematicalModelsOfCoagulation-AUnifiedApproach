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
% TM
% TM==IIa
% PC
% TM==IIa==PC
% APC
% APC==Va
% Va5
% Va3
% APC==Va5
% APC==Va3
% Va53
% HCF
% LCA1
% APC==LCA1
% TM==IIa==APC
% Xa==Va5
% Xa==Va3
% Xa==Va5==II
% Xa==Va3==II
% TM==mIIa
% TM==mIIa==PC
% Xa==Va53
% Xa==Va53==II
% II==Va

function [t,thr,totThr,sol]=main(k,c0,maxt)
% k - Parameters
% c0 - Initial concentrations
options = odeset('AbsTol',1e-15);
sol = ode23tb(@(t,y)f(y,k),[0,maxt],c0,options);
t = sol.x;
y = sol.y;
thr = y(7,:);
totThr = y(7,:)+1.2*y(25,:);
end
function dy = f(y1,k)
dy = zeros(size(y1));
    % TF + VII < TF==VII   k1
    rate = k(1)*y1(3); %TF==VII
    dy([3,1,2]) = dy([3,1,2]) + [-rate;rate;rate];
    % TF + VII > TF==VII   k2
    rate = k(2)*y1(1)*y1(2); %TF + VII
    dy([1,2,3]) = dy([1,2,3]) + [-rate;-rate;rate];
    % TF + VIIa < TF==VIIa   k3
    rate = k(3)*y1(5); %TF==VIIa
    dy([5,1,4]) = dy([5,1,4]) + [-rate;rate;rate];
    % TF + VIIa > TF==VIIa   k4
    rate = k(4)*y1(1)*y1(4); %TF + VIIa
    dy([1,4,5]) = dy([1,4,5]) + [-rate;-rate;rate];
    % TF==VIIa + VII > TF==VIIa + VIIa   k5
    rate = k(5)*y1(5)*y1(2); %TF==VIIa + VII
    dy([2,4]) = dy([2,4]) + [-rate;rate];
    % Xa + VII > Xa + VIIa   k6
    rate = k(6)*y1(6)*y1(2); %Xa + VII
    dy([2,4]) = dy([2,4]) + [-rate;rate];
    % IIa + VII > IIa  + VIIa  k7
    rate = k(7)*y1(7)*y1(2); %IIa + VII
    dy([2,4]) = dy([2,4]) + [-rate;rate];
    % TF==VIIa + X < TF==VIIa==X   k8
    rate = k(8)*y1(9); %TF==VIIa==X
    dy([9,5,8]) = dy([9,5,8]) + [-rate;rate;rate];
    % TF==VIIa + X > TF==VIIa==X   k9
    rate = k(9)*y1(8)*y1(5); %X + TF==VIIa
    dy([5,8,9]) = dy([5,8,9]) + [-rate;-rate;rate];
    % TF==VIIa==X > TF==VIIa==Xa  k10
    rate = k(10)*y1(9); %TF==VIIa==X
    dy([9,10]) = dy([9,10]) + [-rate;rate];
    % TF==VIIa + Xa < TF==VIIa==Xa   k11
    rate = k(11)*y1(10); %TF==VIIa==Xa
    dy([10,5,6]) = dy([10,5,6]) + [-rate;rate;rate];
    % TF==VIIa + Xa > TF==VIIa==Xa   k12
    rate = k(12)*y1(5)*y1(6); %TF==VIIa + Xa
    dy([5,6,10]) = dy([5,6,10]) + [-rate;-rate;rate];
    % TF==VIIa + IX < TF==VIIa==IX   k13
    rate = k(13)*y1(12); %TF==VIIa==IX
    dy([12,5,11]) = dy([12,5,11]) + [-rate;rate;rate];
    % TF==VIIa + IX > TF==VIIa==IX   k14
    rate = k(14)*y1(5)*y1(11); %TF==VIIa + IX
    dy([5,11,12]) = dy([5,11,12]) + [-rate;-rate;rate];
    % TF==VIIa==IX > TF==VIIa + IXa  k15
    rate = k(15)*y1(12); %TF==VIIa==IX
    dy([12,5,13]) = dy([12,5,13]) + [-rate;rate;rate];
    % Xa + II > Xa + IIa   k16
    rate = k(16)*y1(6)*y1(14); %Xa + II
    dy([14,7]) = dy([14,7]) + [-rate;rate];
    % IIa + VIII > IIa + VIIIa   k17
    rate = k(17)*y1(7)*y1(15); %IIa + VIII
    dy([15,16]) = dy([15,16]) + [-rate;rate];
    % VIIIa + IXa < IXa==VIIIa   k18
    rate = k(18)*y1(17); %IXa==VIIIa
    dy([17,13,16]) = dy([17,13,16]) + [-rate;rate;rate];
    % VIIIa + IXa > IXa==VIIIa   k19
    rate = k(19)*y1(16)*y1(13); %VIIIa + IXa
    dy([13,16,17]) = dy([13,16,17]) + [-rate;-rate;rate];
    % IXa==VIIIa + X < IXa==VIIIa==X   k20
    rate = k(20)*y1(18); %IXa==VIIIa==X
    dy([18,8,17]) = dy([18,8,17]) + [-rate;rate;rate];
    % IXa==VIIIa + X > IXa==VIIIa==X   k21
    rate = k(21)*y1(8)*y1(17); %X + IXa==VIIIa
    dy([8,17,18]) = dy([8,17,18]) + [-rate;-rate;rate];
    % IXa==VIIIa==X > IXa==VIIIa + Xa   k22
    rate = k(22)*y1(18); %IXa==VIIIa==X
    dy([18,6,17]) = dy([18,6,17]) + [-rate;rate;rate];
    % VIIIa < VIIIa1-L + VIIIa2   k23
    rate = k(23)*y1(19)*y1(20); %VIIIa1-L + VIIIa2
    dy([19,20,16]) = dy([19,20,16]) + [-rate;-rate;rate];
    % VIIIa > VIIIa1-L + VIIIa2   k24
    rate = k(24)*y1(16); %VIIIa
    dy([16,19,20]) = dy([16,19,20]) + [-rate;rate;rate];
    % IXa==VIIIa==X > VIIIa1-L + VIIIa2 + X + IXa   k25
    rate = k(25)*y1(18); %IXa==VIIIa==X
    dy([18,8,13,19,20]) = dy([18,8,13,19,20]) + [-rate;rate;rate;rate;rate];
    % IXa==VIIIa > VIIIa1-L + VIIIa2 + IXa   k25
    rate = k(25)*y1(17); %IXa==VIIIa
    dy([17,13,19,20]) = dy([17,13,19,20]) + [-rate;rate;rate;rate];
    % IIa + V > IIa + Va   k26
    rate = k(26)*y1(7)*y1(21); %IIa + V
    dy([21,22]) = dy([21,22]) + [-rate;rate];
    % Xa + Va < Xa==Va   k27
    rate = k(27)*y1(23); %Xa==Va
    dy([23,6,22]) = dy([23,6,22]) + [-rate;rate;rate];
    % Xa + Va > Xa==Va   k28
    rate = k(28)*y1(6)*y1(22); %Xa + Va
    dy([6,22,23]) = dy([6,22,23]) + [-rate;-rate;rate];
    % Xa==Va + II < Xa==Va==II   k29
    rate = k(29)*y1(24); %Xa==Va==II
    dy([24,14,23]) = dy([24,14,23]) + [-rate;rate;rate];
    % Xa==Va + II > Xa==Va==II   k30
    rate = k(30)*y1(23)*y1(14); %Xa==Va + II
    dy([14,23,24]) = dy([14,23,24]) + [-rate;-rate;rate];
    % Xa==Va==II > Xa==Va + mIIa   k31
    rate = k(31)*y1(24); %Xa==Va==II
    dy([24,23,25]) = dy([24,23,25]) + [-rate;rate;rate];
    % mIIa + Xa==Va > IIa + Xa==Va   k32
    rate = k(32)*y1(25)*y1(23); %mIIa + Xa==Va
    dy([25,7]) = dy([25,7]) + [-rate;rate];
    % Xa + TFPI < Xa==TFPI   k33
    rate = k(33)*y1(27); %Xa==TFPI
    dy([27,6,26]) = dy([27,6,26]) + [-rate;rate;rate];
    % Xa + TFPI > Xa==TFPI   k34
    rate = k(34)*y1(6)*y1(26); %Xa + TFPI
    dy([6,26,27]) = dy([6,26,27]) + [-rate;-rate;rate];
    % TF==VIIa==Xa + TFPI < TF==VIIa==Xa==TFPI   k35
    rate = k(35)*y1(28); %TF==VIIa==Xa==TFPI
    dy([28,10,26]) = dy([28,10,26]) + [-rate;rate;rate];
    % TF==VIIa==Xa + TFPI > TF==VIIa==Xa==TFPI   k36
    rate = k(36)*y1(10)*y1(26); %TF==VIIa==Xa + TFPI
    dy([10,26,28]) = dy([10,26,28]) + [-rate;-rate;rate];
    % TF==VIIa + Xa==TFPI > TF==VIIa==Xa==TFPI   k37
    rate = k(37)*y1(5)*y1(27); %TF==VIIa + Xa==TFPI
    dy([5,27,28]) = dy([5,27,28]) + [-rate;-rate;rate];
    % Xa + ATIII > Xa==ATIII   k38
    rate = k(38)*y1(29)*y1(6); %ATIII + Xa
    dy([6,29,30]) = dy([6,29,30]) + [-rate;-rate;rate];
    % mIIa + ATIII > mIIa==ATIII   k39
    rate = k(39)*y1(25)*y1(29); %mIIa + ATIII
    dy([25,29,31]) = dy([25,29,31]) + [-rate;-rate;rate];
    % IXa + ATIII > IXa==ATIII   k40
    rate = k(40)*y1(13)*y1(29); %IXa + ATIII
    dy([13,29,32]) = dy([13,29,32]) + [-rate;-rate;rate];
    % IIa + ATIII > IIa==ATIII   k41
    rate = k(41)*y1(7)*y1(29); %IIa + ATIII
    dy([7,29,33]) = dy([7,29,33]) + [-rate;-rate;rate];
    % TF==VIIa + ATIII > TF==VIIa==ATIII   k42
    rate = k(42)*y1(5)*y1(29); %TF==VIIa + ATIII
    dy([5,29,34]) = dy([5,29,34]) + [-rate;-rate;rate];

    %Brummel Reactions
    % TM + IIa < TM==IIa
    rate = k(43)*y1(36);
    dy([36,35,7]) = dy([36,35,7]) + [-rate;rate;rate];

    % TM + IIa > TM==IIa
    rate = k(44)*y1(35)*y1(7);
    dy([35,7,36]) = dy([35,7,36]) + [-rate;-rate;rate];

    % TM==IIa + PC < TM==IIa==PC
    rate = k(45)*y1(38);
    dy([38,36,37]) = dy([38,36,37]) + [-rate;rate;rate];

    % TM==IIa + PC > TM==IIa==PC
    rate = k(46)*y1(36)*y1(37);
    dy([36,37,38]) = dy([36,37,38]) + [-rate;-rate;rate];

    % TM==IIa==PC > TM==IIa + APC
    rate = k(47)*y1(38);
    dy([38,36,39]) = dy([38,36,39]) + [-rate;rate;rate];

    % TM==IIa + ATIII > IIa==ATIII + TM
    rate = k(48)*y1(36)*y1(29);
    dy([29,36,33,35]) = dy([29,36,33,35]) + [-rate;-rate;rate;rate];

    % APC + Va < APC==Va
    rate = k(49)*y1(40);
    dy([40,39,22]) = dy([40,39,22]) + [-rate;rate;rate];

    % APC + Va > APC==Va
    rate = k(50)*y1(22)*y1(39);
    dy([22,39,40]) = dy([22,39,40]) + [-rate;-rate;rate];

    % APC==Va > APC + Va5
    rate = k(51)*y1(40);
    dy([40,39,41]) = dy([40,39,41]) + [-rate;rate;rate];

    % APC==Va > APC + Va3
    rate = k(52)*y1(40);
    dy([40,39,42]) = dy([40,39,42]) + [-rate;rate;rate];

    % APC + Va5 < APC==Va5
    rate = k(49)*y1(43);
    dy([43,39,41]) = dy([43,39,41]) + [-rate;rate;rate];

    % APC + Va5 > APC==Va5
    rate = k(50)*y1(39)*y1(41);
    dy([39,41,43]) = dy([39,41,43]) + [-rate;-rate;rate];

    % APC + Va3 < APC==Va3
    rate = k(49)*y1(44);
    dy([44,39,42]) = dy([44,39,42]) + [-rate;rate;rate];

    % APC + Va3 > APC==Va3
    rate = k(50)*y1(39)*y1(42);
    dy([39,42,44]) = dy([39,42,44]) + [-rate;-rate;rate];

    % APC==Va3 > APC + Va53
    rate = k(51)*y1(44);
    dy([44,39,45]) = dy([44,39,45]) + [-rate;rate;rate];

    % APC==Va5 > APC + Va53
    rate = k(52)*y1(43);
    dy([43,39,45]) = dy([43,39,45]) + [-rate;rate;rate];

    % Va3 > HCF + LCA1
    rate = k(53)*y1(42);
    dy([42,46,47]) = dy([42,46,47]) + [-rate;rate;rate];

    % Va53 > HCF + LCA1
    rate = k(53)*y1(45);
    dy([45,46,47]) = dy([45,46,47]) + [-rate;rate;rate];

    % APC + LCA1 < APC==LCA1
    rate = k(49)*y1(48);
    dy([48,39,47]) = dy([48,39,47]) + [-rate;rate;rate];

    % APC + LCA1 > APC==LCA1
    rate = k(50)*y1(39)*y1(47);
    dy([39,47,48]) = dy([39,47,48]) + [-rate;-rate;rate];

    % APC + TM==IIa < TM==IIa==APC
    rate = k(45)*y1(49);
    dy([49,39,36]) = dy([49,39,36]) + [-rate;rate;rate];

    % APC + TM==IIa > TM==IIa==APC
    rate = k(46)*y1(39)*y1(36);
    dy([36,39,49]) = dy([36,39,49]) + [-rate;-rate;rate];

    % Xa + Va5 < Xa==Va5
    rate = k(54)*y1(50);
    dy([50,6,41]) = dy([50,6,41]) + [-rate;rate;rate];

    % Xa + Va5 > Xa==Va5
    rate = k(28)*y1(6)*y1(41);
    dy([6,41,50]) = dy([6,41,50]) + [-rate;-rate;rate];

    % Xa + Va3 < Xa==Va3
    rate = k(54)*y1(51);
    dy([51,6,42]) = dy([51,6,42]) + [-rate;rate;rate];

    % Xa + Va3 > Xa==Va3
    rate = k(28)*y1(6)*y1(42);
    dy([6,42,51]) = dy([6,42,51]) + [-rate;-rate;rate];

    % Xa==Va5 + II < Xa==Va5==II
    rate = k(29)*y1(52);
    dy([52,50,14]) = dy([52,50,14]) + [-rate;rate;rate];

    % Xa==Va5 + II > Xa==Va5==II
    rate = k(30)*y1(50)*y1(14);
    dy([50,14,52]) = dy([50,14,52]) + [-rate;-rate;rate];

    % Xa==Va5==II > Xa==Va5 + mIIa
    rate = k(55)*y1(52);
    dy([52,50,25]) = dy([52,50,25]) + [-rate;rate;rate];

    % Xa==Va3 + II < Xa==Va3==II
    rate = k(29)*y1(53);
    dy([53,14,51]) = dy([53,14,51]) + [-rate;rate;rate];

    % Xa==Va3 + II > Xa==Va3==II
    rate = k(30)*y1(14)*y1(51);
    dy([14,51,53]) = dy([14,51,53]) + [-rate;-rate;rate];

    % Xa==Va3==II > Xa==Va3 + mIIa
    rate = k(56)*y1(53);
    dy([53,51,25]) = dy([53,51,25]) + [-rate;rate;rate];

    % Xa==Va5 + mIIa > IIa + Xa==Va5
    rate = k(57)*y1(50)*y1(25);
    dy([25,7]) = dy([25,7]) + [-rate;rate];

    % Xa==Va3 + mIIa > IIa + Xa==Va3
    rate = k(58)*y1(51)*y1(25);
    dy([25,7]) = dy([25,7]) + [-rate;rate];

    % Xa==Va3 > HCF + LCA1 + Xa
    rate = k(59)*y1(51);
    dy([51,46,47,6]) = dy([51,46,47,6]) + [-rate;rate;rate;rate];

    % Xa==Va3==II > HCF + LCA1 + Xa + II
    rate = k(59)*y1(53);
    dy([53,46,47,6,14]) = dy([53,46,47,6,14]) + [-rate;rate;rate;rate;rate];

    % IXa + X > IXa + Xa
    rate = k(60)*y1(13)*y1(8);
    dy([8,6]) = dy([8,6]) + [-rate;rate];

    % mIIa + V > mIIa + Va
    rate = k(61)*y1(25)*y1(21);
    dy([21,22]) = dy([21,22]) + [-rate;rate];

    % TM + mIIa < TM==mIIa
    rate = k(43)*y1(54);
    dy([54,35,25]) = dy([54,35,25]) + [-rate;rate;rate];

    % TM + mIIa > TM==mIIa
    rate = k(44)*y1(35)*y1(25);
    dy([25,35,54]) = dy([25,35,54]) + [-rate;-rate;rate];

    % TM==mIIa + PC < TM==mIIa==PC
    rate = k(45)*y1(55);
    dy([55,54,37]) = dy([55,54,37]) + [-rate;rate;rate];

    % TM==mIIa + PC > TM==mIIa==PC
    rate = k(46)*y1(54)*y1(37);
    dy([54,37,55]) = dy([54,37,55]) + [-rate;-rate;rate];

    % TM==mIIa==PC > TM==mIIa + APC
    rate = k(47)*y1(55);
    dy([55,54,39]) = dy([55,54,39]) + [-rate;rate;rate];

    % TM==mIIa + ATIII > mIIa==ATIII + TM
    rate = k(48)*y1(54)*y1(29);
    dy([54,29,31,35]) = dy([54,29,31,35]) + [-rate;-rate;rate;rate];

    % Xa + Va53 < Xa==Va53
    rate = k(54)*y1(56);
    dy([56,6,45]) = dy([56,6,45]) + [-rate;rate;rate];

    % Xa + Va53 > Xa==Va53
    rate = k(28)*y1(6)*y1(45);
    dy([6,45,56]) = dy([6,45,56]) + [-rate;-rate;rate];

    % Xa==Va53 + II < Xa==Va53==II
    rate = k(29)*y1(57);
    dy([57,56,14]) = dy([57,56,14]) + [-rate;rate;rate];

    % Xa==Va53 + II > Xa==Va53==II
    rate = k(30)*y1(14)*y1(56);
    dy([14,56,57]) = dy([14,56,57]) + [-rate;-rate;rate];

    % Xa==Va53==II > Xa==Va53 + mIIa
    rate = k(56)*y1(57);
    dy([57,56,25]) = dy([57,56,25]) + [-rate;rate;rate];

    % Xa==Va53 + mIIa > Xa==Va53 + IIa
    rate = k(58)*y1(56)*y1(25);
    dy([25,7]) = dy([25,7]) + [-rate;rate];

    % Xa==Va53 > HCF + LCA1 + Xa
    rate = k(59)*y1(56);
    dy([56,46,47,6]) = dy([56,46,47,6]) + [-rate;rate;rate;rate];

    % Xa==Va53==II > HCF + LCA1 + Xa + II
    rate = k(59)*y1(57);
    dy([57,46,47,6,14]) = dy([57,46,47,6,14]) + [-rate;rate;rate;rate;rate];

    % II + Va < II==Va
    rate = k(62)*y1(58);
    dy([58,14,22]) = dy([58,14,22]) + [-rate;rate;rate];

    % II + Va > II==Va
    rate = k(63)*y1(14)*y1(22);
    dy([14,22,58]) = dy([14,22,58]) + [-rate;-rate;rate];
    
    % Xa==Va5 + APC > Xa==Va53 + APC
    rate = k(64)*y1(50)*y1(39);
    dy([50,56]) = dy([50,56]) + [-rate;rate];
end