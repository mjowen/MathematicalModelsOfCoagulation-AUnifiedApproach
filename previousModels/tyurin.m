% TF
% TFVII
% TFVIIa
% VIIa
% VII
% Xa
% X
% Va
% V
% XaVa
% XIa
% XI
% IXa
% IX
% VIIIa
% VIII
% IXaVIIIa
% IIa
% II
% TM
% TMIIa
% PCa
% PC
% TFPI
% TFPIXa
% ATIII
% C1
% PAI1
% alpha1AT
% alpha2AP
% alpha2M
% PCI
% XIIa

function [t,thr,y,sol]=main(k,c0,maxt)
% k - Parameters
% c0 - Initial concentrations
options = odeset('AbsTol',1e-10);
sol = ode23tb(@(t,y)f(y,k),[0,maxt/60],c0/10^-6,options);
t = sol.x;
y = sol.y;
thr = y(18,:)*10^-6;
t = t*60;
y = y*10^-6;
end

function dy = f(y,k)
dy = zeros(size(y));

%Enzyme Reactions
% 1 XIIa-XI
rate = k(29)*y(12)*y(33)/( k(30) +y(12) );
dy([12,11]) = dy([12,11]) + [-rate;rate];
% 2 IIa-XI
rate = k(31)*y(12)*y(18)/( k(32)* (1+y(9)/k(54)+y(16)/k(62)) +y(12) );
dy([12,11]) = dy([12,11]) + [-rate;rate];
% 3 XIa-XI
rate = k(33)*y(12)*y(11)/( k(34)* (1+y(14)/k(36)) +y(12) );
dy([12,11]) = dy([12,11]) + [-rate;rate];
% 4 XIa-IX
rate = k(35)*y(11)*y(14)/( k(36)* (1+y(12)/k(34)) +y(14) );
dy([14,13]) = dy([14,13]) + [-rate;rate];
% 5 VIIa-IX
rate = k(37)*y(4)*y(14)/( k(38)* (1+y(7)/k(46)) +y(14) );
dy([14,13]) = dy([14,13]) + [-rate;rate];
% 6 TFVIIa-IX
rate = k(39)*y(3)*y(14)/( k(40)* (1+y(7)/k(48)) +y(14) );
dy([14,13]) = dy([14,13]) + [-rate;rate];
% 7 IXa-X
rate = k(41)*y(7)*y(13)/( k(42) +y(7) );
dy([7,6]) = dy([7,6]) + [-rate;rate];
% 8 IXaVIIIa-X
rate = k(43)*y(7)*y(17)/( k(44) +y(7) );
dy([7,6]) = dy([7,6]) + [-rate;rate];
% 9 VIIa-X
rate = k(45)*y(4)*y(7)/( k(46)* (1+y(14)/k(38)) +y(7) );
dy([7,6]) = dy([7,6]) + [-rate;rate];
% 10 TFVIIa-X
rate = k(47)*y(3)*y(7)/( k(48)* (1+y(14)/k(40)) +y(7) );
dy([7,6]) = dy([7,6]) + [-rate;rate];
% 11 Xa-II
rate = k(49)*y(6)*y(19)/( k(50)* (1+y(9)/k(56) +y(5)/k(58) +y(2)/k(60)) +y(19) );
dy([19,18]) = dy([19,18]) + [-rate;rate];
% 12 XaVa-II
rate = k(51)*y(10)*y(19)/( k(52) +y(19) );
dy([19,18]) = dy([19,18]) + [-rate;rate];
% 13 IIa-V
rate = k(53)*y(18)*y(9)/( k(54)* (1+y(12)/k(32)+y(16)/k(62)) +y(9) );
dy([9,8]) = dy([9,8]) + [-rate;rate];
% 14 Xa-V
rate = k(55)*y(6)*y(9)/( k(56)* (1+y(19)/k(50)+y(5)/k(58)+y(2)/k(60)) +y(9) );
dy([9,8]) = dy([9,8]) + [-rate;rate];
% 15 Xa-VII
rate = k(57)*y(6)*y(5)/( k(58)* (1+y(9)/k(56)+y(19)/k(50)+y(2)/k(60)) +y(5) );
dy([5,4]) = dy([5,4]) + [-rate;rate];
% 16 Xa-TFVII
rate = k(59)*y(2)*y(6)/( k(60)* (1+y(9)/k(56)+y(19)/k(50)+y(5)/k(58)) +y(2) );
dy([2,3]) = dy([2,3]) + [-rate;rate];
% 17 IIa-VIII
rate = k(61)*y(16)*y(18)/( k(62)* (1+y(12)/k(32) +y(9)/k(54)) +y(16) );
dy([16,15]) = dy([16,15]) + [-rate;rate];
% 18 TMIIa-PC
rate = k(63)*y(21)*y(23)/( k(64) +y(23) );
dy([23,22]) = dy([23,22]) + [-rate;rate];
% 19 PC (Va)
rate = k(65)*y(22)*y(8)/( k(66) +y(15) +y(8) +y(10) +y(17) );
dy(8) = dy(8) -rate;
% 20 PC (VIIIa)
rate = k(65)*y(22)*y(15)/( k(66) +y(15) +y(8) +y(10) +y(17) );
dy(15) = dy(15) -rate;
% 21 PC (IXaVIIIa)
rate = k(65)*y(22)*y(17)/( k(66) +y(15) +y(8) +y(10) +y(17) );
dy([13,17]) = dy([13,17]) + [rate;-rate];
% 22 PC (XaVa)
rate = k(65)*y(22)*y(10)/( k(66) +y(15)+y(8)+y(10)+y(17) );
dy([10,6]) = dy([10,6]) + [-rate;rate];

%Va + Xa > Xa==Va
rate = k(1)*y(8)*y(6);
dy([6,8,10]) = dy([6,8,10]) + [-rate;-rate;rate];
%VIIIa + IXa > VIIIa==IXa
rate = k(2)*y(15)*y(13); %IXa instead of XIa
dy([13,15,17]) = dy([13,15,17]) + [-rate;-rate;rate];
%VIIa + TF > TF==VIIa
rate = k(3)*y(4)*y(1);
dy([4,1,3]) = dy([4,1,3]) + [-rate;-rate;rate];
%VII + TF > TF==VII
rate = k(4)*y(5)*y(1);
dy([5,1,2]) = dy([5,1,2]) + [-rate;-rate;rate];
%TF==VIIa + TFPI==Xa > TF==VIIa==TFPI==Xa
rate = k(5)*y(3)*y(25);
dy([3,25]) = dy([3,25]) + [-rate;-rate];
%TF==VIIa + ATIII > TF==VIIa==ATIII
rate = k(6)*y(3)*y(26);
dy([3,26]) = dy([3,26]) + [-rate;-rate];
%IIa + ATIII > IIa==ATIII
rate = k(7)*y(18)*y(26);
dy([18,26]) = dy([18,26]) + [-rate;-rate];
%IIa + a1AT > IIa==a1AT
rate = k(8)*y(18)*y(29);
dy([18,29]) = dy([18,29]) + [-rate;-rate];
%IIa + a2M > IIa==a2M
rate = k(9)*y(18)*y(31);
dy([18,31]) = dy([18,31]) + [-rate;-rate];
%IIa + PCI > IIa==PCI
rate = k(10)*y(18)*y(32);
dy([18,32]) = dy([18,32]) + [-rate;-rate];
%Xa + ATIII > Xa==ATIII
rate = k(11)*y(6)*y(26);
dy([6,26]) = dy([6,26]) + [-rate;-rate];
%Xa + a1AT > Xa==a1AT
rate = k(12)*y(6)*y(29);
dy([6,29]) = dy([6,29]) + [-rate;-rate];
%Xa + TFPI > Xa==TFPI
rate = k(13)*y(6)*y(24);
dy([6,24,25]) = dy([6,24,25]) + [-rate;-rate;rate];
%XaVa + a1AT > Xa==a1AT + Va
rate = k(14)*y(10)*y(29);
dy([10,29,8]) = dy([10,29,8]) + [-rate;-rate;rate];
%XaVa + ATIII > Xa==ATIII + Va
rate = k(15)*y(10)*y(26);
dy([10,26,8]) = dy([10,26,8]) + [-rate;-rate;rate];
%IXa + ATIII > IXa==ATIII
rate = k(16)*y(13)*y(26);
dy([13,26]) = dy([13,26]) + [-rate;-rate];
%VIIIa==IXa + ATIII > IXa==ATIII + VIIIa
rate = k(17)*y(17)*y(26);
dy([17,26,15]) = dy([17,26,15]) + [-rate;-rate;rate];
%XIa + C1-Inh > XI==C1-Inh
rate = k(18)*y(11)*y(27);
dy([11,27]) = dy([11,27]) + [-rate;-rate];
%XIa + a1AT > XI==a1AT
rate = k(19)*y(11)*y(29);
dy([11,29]) = dy([11,29]) + [-rate;-rate];
%XIa + ATIII > XIa==ATIII
rate = k(20)*y(11)*y(26);
dy([11,26]) = dy([11,26]) + [-rate;-rate];
%XIa + a2AP > XIa==a2AP
rate = k(21)*y(11)*y(30);
dy([11,30]) = dy([11,30]) + [-rate;-rate];
%XIa + PAI1 > XIa==PAI1
rate = k(22)*y(11)*y(28);
dy([11,28]) = dy([11,28]) + [-rate;-rate];
%IIa + TM > IIa==TM
rate = k(23)*y(18)*y(20);
dy([18,20,21]) = dy([18,20,21]) + [-rate;-rate;rate];
%IIa==TM + PCI > IIa==TM==PCI
rate = k(24)*y(21)*y(32);
dy([21,32]) = dy([21,32]) + [-rate;-rate];
%PCa + PCI > PCa==PCI
rate = k(25)*y(22)*y(32);
dy([22,32]) = dy([22,32]) + [-rate;-rate];
%PCa + a1AT > PCa==a1AT
rate = k(26)*y(22)*y(29);
dy([22,29]) = dy([22,29]) + [-rate;-rate];
%TF==VIIa > TF + VIIa
rate = k(27)*y(3);
dy([3,1,4]) = dy([3,1,4]) + [-rate;rate;rate];
%TF==VII > TF + VII
rate = k(28)*y(2);
dy([2,1,5]) = dy([2,1,5]) + [-rate;rate;rate];

end
