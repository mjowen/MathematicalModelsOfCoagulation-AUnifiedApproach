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
% Boc-VPR-MCA
% Boc-VPR-MCA-IIa
% Boc-VPR
% AMC
% XII
% XIIa
% XIIa==XII
% PK
% XIIa==PK
% K
% XII==K
% Kinh
% CTI
% C1inh
% XIIa==C1inh
% XIIa==ATIII
% XI
% XI-IIa
% XIa
% XIIa==XI
% XIa==ATIII
% XIa==C1inh
% a1AT
% XIa==a1AT
% a2AP
% XIa==a2AP
% XIa==IX
% IXa==X
% Xa==VIII
% VIIa==IX
% VIIa==X
% Fbg
% Fbg==IIa
% Fbn1
% FPA
% Fbn2
% Fbn12
% Fbn22
% FPB
% Fbn2==IIa
% Fbn12==IIa
% Fbn12==IIa==ATIII
% Fbn1==IIa
% Fbn1==IIa==ATIII
% Fbn2==IIa==ATIII
% XIIa==CTI
% eps

function [t,thr,totThr,sol,eps,y]=main(k,c0,maxt)
% k - Parameters
% c0 - Initial concentrations
% maxt - End of time interval
global IIaMax
IIaMax = 0;
options = odeset('AbsTol',1e-16);
sol = ode23tb(@(t,y)f(y,k),[0,maxt],c0,options);
t = sol.x;
y = sol.y;
thr = y(7,:);
totThr = y(7,:)+1.2*y(25,:);
eps = y(81,:);
end

function dy = f(y1,k)
global IIaMax
IIaMax = max(IIaMax,y1(7));
fIIaMax = IIaMax^1.6123/(IIaMax^1.6123+(2.4279e-9)^1.6123);
epsMax = k(104)+(1-k(104))*fIIaMax;
dy = zeros(size(y1));
dy(81) = k(105)*(epsMax-y1(81));

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
rate = k(18)*y1(17)/y1(81); %IXa==VIIIa
dy([17,13,16]) = dy([17,13,16]) + [-rate;rate;rate];
% VIIIa + IXa > IXa==VIIIa   k19
rate = k(19)*y1(16)*y1(13); %VIIIa + IXa
dy([13,16,17]) = dy([13,16,17]) + [-rate;-rate;rate];
% IXa==VIIIa + X < IXa==VIIIa==X   k20
rate = k(20)*y1(18)/y1(81); %IXa==VIIIa==X
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
rate = k(24)*y1(16)/y1(81); %VIIIa
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
rate = k(27)*y1(23)/y1(81); %Xa==Va
dy([23,6,22]) = dy([23,6,22]) + [-rate;rate;rate];
% Xa + Va > Xa==Va   k28
rate = k(28)*y1(6)*y1(22); %Xa + Va
dy([6,22,23]) = dy([6,22,23]) + [-rate;-rate;rate];
% Xa==Va + II < Xa==Va==II   k29
rate = k(29)*y1(24)/y1(81); %Xa==Va==II
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
rate = k(33)*y1(27)/y1(81); %Xa==TFPI
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


% Boc-VPR-MCA + IIa > Boc-VPR-MCA==IIa
rate = k(43)*y1(35)*y1(7);
dy([7,35,36]) = dy([7,35,36]) + [-rate;-rate;rate];
% Boc-VPR-MCA + IIa < Boc-VPR-MCA==IIa
rate = k(44)*y1(36);
dy([36,7,35]) = dy([36,7,35]) + [-rate;rate;rate];
% Boc-VPR-MCA==IIa > Boc-VPR + AMC + IIa
rate = k(45)*y1(36);
dy([36,37,38,7]) = dy([36,37,38,7]) + [-rate;rate;rate;rate];
% XII > XIIa
rate = k(46)*y1(39);
dy([39,40]) = dy([39,40]) + [-rate;rate];
% XIIa + XII > XIIa==XII
rate = k(47)*y1(39)*y1(40);
dy([39,40,41]) = dy([39,40,41]) + [-rate;-rate;rate];
% XIIa + XII < XIIa==XII
rate = k(48)*y1(41)/y1(81);
dy([41,39,40]) = dy([41,39,40]) + [-rate;rate;rate];
% XIIa==XII > 2 XIIa
rate = k(49)*y1(41);
dy([41,40]) = dy([41,40]) + [-rate;2*rate];
%XIIa + PK > XIIa==PK
rate = k(50)*y1(40)*y1(42);
dy([40,42,43]) = dy([40,42,43]) + [-rate;-rate;rate];
%XIIa + PK < XIIa==PK
rate = k(51)*y1(43)/y1(81);
dy([43,40,42]) = dy([43,40,42]) + [-rate;rate;rate];
% XIIa==PK > XIIa + K
rate = k(52)*y1(43);
dy([43,40,44]) = dy([43,40,44]) + [-rate;rate;rate];
% XII + K > XII==K
rate = k(53)*y1(39)*y1(44);
dy([39,44,45]) = dy([39,44,45]) + [-rate;-rate;rate];
%XII + K < XII==K
rate = k(54)*y1(45)/y1(81);
dy([45,39,44]) = dy([45,39,44]) + [-rate;rate;rate];
%XII==K > XIIa + K
rate = k(55)*y1(45);
dy([45,40,44]) = dy([45,40,44]) + [-rate;rate;rate];
% PK + K > K + K
rate = k(56)*y1(42)*y1(44);
dy([42,44]) = dy([42,44]) + [-rate;rate];
% K > Kinh
rate = k(57)*y1(44);
dy([44,46]) = dy([44,46]) + [-rate;rate];
% XIIa + CTI > XIIa==CTI
rate = k(58)*y1(40)*y1(47);
dy([40,47,80]) = dy([40,47,80]) + [-rate;-rate;rate];
% XIIa + CTI < XIIa==CTI
rate = k(59)*y1(80);
dy([80,40,47]) = dy([80,40,47]) + [-rate;rate;rate];
% XIIa + C1inh > XIIa==C1inh
rate = k(60)*y1(40)*y1(48);
dy([40,48,49]) = dy([40,48,49]) + [-rate;-rate;rate];
% XIIa + ATIII > XIIa==ATIII
rate = k(61)*y1(40)*y1(29);
dy([40,29,50]) = dy([40,29,50]) + [-rate;-rate;rate];
% XI + IIa > XI==IIa
rate = k(62)*y1(51)*y1(7);
dy([51,7,52]) = dy([51,7,52]) + [-rate;-rate;rate];
% XI + IIa < XI==IIa
rate = k(63)*y1(52);
dy([52,7,51]) = dy([52,7,51]) + [-rate;rate;rate];
% XI==IIa > XIa + IIa
rate = k(64)*y1(52);
dy([52,7,53]) = dy([52,7,53]) + [-rate;rate;rate];
% XIIa + XI > XIIa==XI
rate = k(65)*y1(40)*y1(51);
dy([40,51,54]) = dy([40,51,54]) + [-rate;-rate;rate];
% XIIa + XI < XIIa==XI
rate = k(66)*y1(54)/y1(81);
dy([54,40,51]) = dy([54,40,51]) + [-rate;rate;rate];
% XIIa==XI > XIIa + XIa
rate = k(67)*y1(54);
dy([54,40,53]) = dy([54,40,53]) + [-rate;rate;rate];
% XIa + XI > 2 XIa
rate = k(68)*y1(53)*y1(51);
dy([51,53]) = dy([51,53]) + [-rate;rate];
%XIa + ATIII > XIa==ATIII
rate = k(69)*y1(53)*y1(29);
dy([53,29,55]) = dy([53,29,55]) + [-rate;-rate;rate];
% XIa + C1inh > XIa==C1inh
rate = k(70)*y1(53)*y1(48);
dy([53,48,56]) = dy([53,48,56]) + [-rate;-rate;rate];
% XIa + a1AT > XIa==a1AT
rate = k(71)*y1(53)*y1(57);
dy([53,57,58]) = dy([53,57,58]) + [-rate;-rate;rate];
%XIa + a2AP > XIa==a2AP
rate = k(72)*y1(53)*y1(59);
dy([53,59,60]) = dy([53,59,60]) + [-rate;-rate;rate];
% XIa + IX > XIa==IX
rate = k(73)*y1(53)*y1(11); %Should be y1(11) not y1(51)
dy([53,11,61]) = dy([53,11,61]) + [-rate;-rate;rate];
% XIa + IX < XIa==IX
rate = k(74)*y1(61)/y1(81);
dy([61,53,11]) = dy([61,53,11]) + [-rate;rate;rate];
% XIa==IX > XIa + IXa
rate = k(75)*y1(61);
dy([61,53,13]) = dy([61,53,13]) + [-rate;rate;rate];
% IXa + X > IXa==X
rate = k(76)*y1(13)*y1(8);
dy([13,8,62]) = dy([13,8,62]) + [-rate;-rate;rate];
% IXa + X < IXa==X
rate = k(77)*y1(62)/y1(81);
dy([62,8,13]) = dy([62,8,13]) + [-rate;rate;rate];
% IXa==X > IXa + Xa
rate = k(78)*y1(62);
dy([62,13,6]) = dy([62,13,6]) + [-rate;rate;rate];
% Xa + VIII > Xa==VIII
rate = k(79)*y1(6)*y1(15);
dy([6,15,63]) = dy([6,15,63]) + [-rate;-rate;rate];
% Xa + VIII < Xa==VIII
rate = k(80)*y1(63)/y1(81);
dy([63,6,15]) = dy([63,6,15]) + [-rate;rate;rate];
% Xa==VIII > Xa + VIIIa
rate = k(81)*y1(63);
dy([63,6,16]) = dy([63,6,16]) + [-rate;rate;rate];
% VIIa + IX > VIIa==IX
rate = k(82)*y1(4)*y1(11);
dy([4,11,64]) = dy([4,11,64]) + [-rate;-rate;rate];
% VIIa + IX < VIIa==IX
rate = k(83)*y1(64);
dy([64,11,4]) = dy([64,11,4]) + [-rate;rate;rate];
% VIIa==IX > VIIa + IXa
rate = k(84)*y1(64);
dy([64,4,13]) = dy([64,4,13]) + [-rate;rate;rate];
% VIIa + X > VIIa==X
rate = k(85)*y1(4)*y1(8);
dy([4,8,65]) = dy([4,8,65]) + [-rate;-rate;rate];
% VIIa + X < VIIa==X
rate = k(86)*y1(65)/y1(81);
dy([65,4,8]) = dy([65,4,8]) + [-rate;rate;rate];
% VIIa==X > VIIa + Xa
rate = k(87)*y1(65);
dy([65,4,6]) = dy([65,4,6]) + [-rate;rate;rate];
% Fbg + IIa > Fbg==IIa
rate = k(88)*y1(66)*y1(7);
dy([66,7,67]) = dy([66,7,67]) + [-rate;-rate;rate];
% Fbg + IIa < Fbg==IIa
rate = k(89)*y1(67);
dy([67,7,66]) = dy([67,7,66]) + [-rate;rate;rate];
% Fbg==IIa > Fbn1 + IIa + FPA
rate = k(90)*y1(67);
dy([67,7,68,69]) = dy([67,7,68,69]) + [-rate;rate;rate;rate];
% Fbn1 + IIa > Fbn1==IIa
rate = k(91)*y1(7)*y1(68);
dy([68,7,77]) = dy([68,7,77]) + [-rate;-rate;rate];
% Fbn1 + IIa < Fbn1==IIa
rate = k(92)*y1(77);
dy([77,7,68]) = dy([77,7,68]) + [-rate;rate;rate];
% Fbn1==IIa > Fbn2 + IIa + FPB
rate = k(93)*y1(77);
dy([77,7,70,73]) = dy([77,7,70,73]) + [-rate;rate;rate;rate];
% 2 Fbn1 > Fbn12
rate = k(94)*y1(68)^2;
dy([68,71]) = dy([68,71]) + [-2*rate;rate];
% 2 Fbn1 < Fbn12
rate = k(95)*y1(71);
dy([71,68]) = dy([71,68]) + [-rate;2*rate];
% Fbn12 + IIa > Fbn12==IIa
rate = k(96)*y1(71)*y1(7);
dy([71,7,75]) = dy([71,7,75]) + [-rate;-rate;rate];
% Fbn12 + IIa < Fbn12==IIa
rate = k(97)*y1(75);
dy([75,7,71]) = dy([75,7,71]) + [-rate;rate;rate];
% Fbn12==IIa > Fbn22 + IIa + FPB
rate = k(98)*y1(75);
dy([75,7,72,73]) = dy([75,7,72,73]) + [-rate;rate;rate;rate];
% Fbn2 + IIa > Fbn2==IIa
rate = k(99)*y1(7)*y1(70);
dy([70,7,74]) = dy([70,7,74]) + [-rate;-rate;rate];
% Fbn2 + IIa < Fbn2==IIa
rate = k(100)*y1(74);
dy([74,7,70]) = dy([74,7,70]) + [-rate;rate;rate];
% Fbn12==IIa + ATIII > Fbn12==IIa==ATIII
rate = k(101)*y1(75)*y1(29);
dy([75,29,76]) = dy([75,29,76]) + [-rate;-rate;rate];
% Fbn1==IIa + ATIII > Fbn1==IIa==ATIII
rate = k(102)*y1(77)*y1(29);
dy([77,29,78]) = dy([77,29,78]) + [-rate;-rate;rate];
% Fbn2==IIa + ATIII > Fbn2==IIa==ATIII
rate = k(103)*y1(74)*y1(29);
dy([74,29,79]) = dy([74,29,79]) + [-rate;-rate;rate];

end
