% IIf
% IIL
% mIIaf
% mIIaL
% Vf
% VL
% Vaf
% VaL
% VIIf
% VIIL
% VIIaf
% VIIaL
% VIIIf
% VIIIL
% VIIIaf
% VIIIaL
% IXf
% IXL
% IXaf
% IXaL
% Xf
% XL
% Xaf
% XaL
% APCf
% APCL
% PSf
% PSL
% VIIIaif
% VIIIaiL
% Vaif
% VaiL
% PCf
% PCL
% TFL
% TFVIIaL
% TFVIIL
% TFVIIaIXL
% TFVIIaXL
% TFVIIaXaL
% TFVIIXaL
% IXaVIIIaL
% XaVaL
% IXaVIIIaXL
% VXaL
% VIIIXaL
% IIaf
% V_IIaL
% VIII_IIaL
% XaVaIIL
% XaVamIIaL
% XIf
% XI_IIaf
% XIaf
% APCPSL
% APCPSVIIIaL
% TFPIf
% ATf
% IIaATf
% TFPIXaf
% TFPIXaTFVIIaL
% APCPSVaL
% IXaATf
% XaATf
% VIIXaL
% VmIIaL
% VIIImIIaL
% TML
% IIaTML
% IIaTMPCL
% mIIaATL
% XIaIXL
% LIPID

function [t,thr,totThr,sol]=bungay(k,c0,maxt)
% k - Parameters
% c0 - Initial concentrations
options = odeset('AbsTol',1e-16);
sol = ode23tb(@(t,y)f(y,k),[0,maxt],c0/10^-9,options);
t = sol.x;
y = sol.y;
thr = (y(47,:))*10^-9;
totThr = (y(47,:)+1.2*y(3,:)+1.2*y(4,:))*10^-9;
end



function dy=f(y,k)
kon = k(76:92);
koff = k(93:end-1);
nva = k(end);
dy = zeros(size(y));
lipidChange = 0;
% II lipid binding
rate=(kon(1)*y(1)*y(73)/nva-koff(1)*y(2));
dy([1,2]) = dy([1,2]) + [-rate;rate];
lipidChange = lipidChange + rate;
% mIIa lipid binding
rate=(kon(2)*y(3)*y(73)/nva-koff(2)*y(4));
dy([3,4]) = dy([3,4]) + [-rate;rate];
lipidChange = lipidChange + rate;
% V lipid binding
rate=(kon(3)*y(5)*y(73)/nva-koff(3)*y(6));
dy([5,6]) = dy([5,6]) + [-rate;rate];
lipidChange = lipidChange + rate;
% Va lipid binding
rate=(kon(4)*y(7)*y(73)/nva-koff(4)*y(8));
dy([7,8]) = dy([7,8]) + [-rate;rate];
lipidChange = lipidChange + rate;
% VII lipid binding
rate=(kon(5)*y(9)*y(73)/nva-koff(5)*y(10));
dy([9,10]) = dy([9,10]) + [-rate;rate];
lipidChange = lipidChange + rate;
% VIIa lipid binding
rate=(kon(6)*y(11)*y(73)/nva-koff(6)*y(12));
dy([11,12]) = dy([11,12]) + [-rate;rate];
lipidChange = lipidChange + rate;
% VIII lipid binding
rate=(kon(7)*y(13)*y(73)/nva-koff(7)*y(14));
dy([13,14]) = dy([13,14]) + [-rate;rate];
lipidChange = lipidChange + rate;
% VIIIa lipid binding
rate=(kon(8)*y(15)*y(73)/nva-koff(8)*y(16));
dy([15,16]) = dy([15,16]) + [-rate;rate];
lipidChange = lipidChange + rate;
% IX lipid binding
rate=(kon(9)*y(17)*y(73)/nva-koff(9)*y(18));
dy([17,18]) = dy([17,18]) + [-rate;rate];
lipidChange = lipidChange + rate;
% IXa lipid binding
rate=(kon(10)*y(19)*y(73)/nva-koff(10)*y(20));
dy([19,20]) = dy([19,20]) + [-rate;rate];
lipidChange = lipidChange + rate;
% X lipid binding
rate=(kon(11)*y(21)*y(73)/nva-koff(11)*y(22));
dy([21,22]) = dy([21,22]) + [-rate;rate];
lipidChange = lipidChange + rate;
% Xa lipid binding
rate=kon(12)*y(23)*y(73)/nva-koff(12)*y(24);
dy([23,24]) = dy([23,24]) + [-rate;rate];
lipidChange = lipidChange + rate;
% APC lipid binding
rate=(kon(13)*y(25)*y(73)/nva-koff(13)*y(26));
dy([25,26]) = dy([25,26]) + [-rate;rate];
lipidChange = lipidChange + rate;
% PS lipid binding
rate=(kon(14)*y(27)*y(73)/nva-koff(14)*y(28));
dy([27,28]) = dy([27,28]) + [-rate;rate];
lipidChange = lipidChange + rate;
% VIIIai lipid binding
rate=(kon(15)*y(29)*y(73)/nva-koff(15)*y(30));
dy([29,30]) = dy([29,30]) + [-rate;rate];
lipidChange = lipidChange + rate;
% Vai lipid binding
rate=(kon(16)*y(31)*y(73)/nva-koff(16)*y(32));
dy([31,32]) = dy([31,32]) + [-rate;rate];
lipidChange = lipidChange + rate;
% PC lipid binding
rate=(kon(17)*y(33)*y(73)/nva-koff(17)*y(34));
dy([33,34]) = dy([33,34]) + [-rate;rate];
lipidChange = lipidChange + rate;


% TF + VIIa <-> TF:VIIa
rate=(k(1)*y(35)*y(12)-k(2)*y(36));
dy([35,12,36]) = dy([35,12,36]) + [-rate;-rate;rate];
% TF + VII <-> TF:VII
rate=(k(3)*y(35)*y(10)-k(4)*y(37));
dy([10,35,37]) = dy([10,35,37]) + [-rate;-rate;rate];
% IX + TF:VIIa <-> IX:TF:VIIa
rate=(k(5)*y(36)*y(18)-k(6)*y(38));
dy([18,36,38]) = dy([18,36,38]) + [-rate;-rate;rate];
% IX:TF:VIIa -> IXa + TF:VIIa
rate=k(7)*y(38);
dy([38,36,20]) = dy([38,36,20]) + [-rate;rate;rate];
% X + TF:VIIa <-> X:TF:VIIa
rate=(k(8)*y(36)*y(22)-k(9)*y(39));
dy([22,36,39]) = dy([22,36,39]) + [-rate;-rate;rate];
% X:TF:VIIa -> Xa:TF:VIIa
rate=k(10)*y(39);
dy([39,40]) = dy([39,40]) + [-rate;rate];
% Xa:TF:VIIa -> Xa + TF:VIIa
rate=k(75)*y(40);
dy([40,36,24]) = dy([40,36,24]) + [-rate;rate;rate];
% TF:VII + Xa <-> Xa:TF:VII
rate=(k(11)*y(37)*y(24)-k(12)*y(41));
dy([24,37,41]) = dy([24,37,41]) + [-rate;-rate;rate];
% Xa:TF:VII -> Xa + TF:VIIa
rate=k(13)*y(41);
dy([41,24,36]) = dy([41,24,36]) + [-rate;rate;rate];
% IXa + VIIIa <-> IXa:VIIIa
rate=(k(14)*y(20)*y(16)-k(15)*y(42));
dy([16,20,42]) = dy([16,20,42]) + [-rate;-rate;rate];
% Xa + Va <-> Xa:Va
rate=(k(16)*y(24)*y(8)-k(17)*y(43));
dy([8,24,43]) = dy([8,24,43]) + [-rate;-rate;rate];
% X + IXa:VIIIa <-> X:IXa:VIIIa
rate=(k(18)*y(42)*y(22)-k(19)*y(44));
dy([22,42,44]) = dy([22,42,44]) + [-rate;-rate;rate];
% X:IXa:VIIIa -> Xa + IXa:VIIIa
rate=k(20)*y(44);
dy([44,24,42]) = dy([44,24,42]) + [-rate;rate;rate];
% V + Xa <-> V:Xa
rate=(k(21)*y(6)*y(24)-k(22)*y(45));
dy([6,24,45]) = dy([6,24,45]) + [-rate;-rate;rate];
% V:Xa -> Va + Xa
rate=k(23)*y(45);
dy([45,8,24]) = dy([45,8,24]) + [-rate;rate;rate];
% Xa + VIII <-> Xa:VIII
rate=k(24)*y(14)*y(24)-k(25)*y(46);
dy([14,24,46]) = dy([14,24,46]) + [-rate;-rate;rate];
% Xa:VIII -> Xa + VIIIa
rate=k(26)*y(46);
dy([46,16,24]) = dy([46,16,24]) + [-rate;rate;rate];
% V + IIa <-> V:IIa
rate=(k(27)*y(6)*y(47)-k(28)*y(48));
dy([6,47,48]) = dy([6,47,48]) + [-rate;-rate;rate];
% V:IIa -> Va:IIa
rate=k(29)*y(48);
dy([48,8,47]) = dy([48,8,47]) + [-rate;rate;rate];
% VIII + IIa <-> VIII:IIa
rate=(k(30)*y(14)*y(47)-k(31)*y(49));
dy([14,47,49]) = dy([14,47,49]) + [-rate;-rate;rate];
% VIII:IIa -> VIIIa + IIa
rate=k(32)*y(49);
dy([49,16,47]) = dy([49,16,47]) + [-rate;rate;rate];
% Xa:Va + II <-> Xa:Va:II
rate=(k(33)*y(43)*y(2)-k(34)*y(50));
dy([2,43,50]) = dy([2,43,50]) + [-rate;-rate;rate];
% Xa:Va + mIIa <-> Xa:Va:mIIa
rate=(k(35)*y(43)*y(4)-k(36)*y(51));
dy([4,43,51]) = dy([4,43,51]) + [-rate;-rate;rate];
% Xa:Va:II -> Xa:Va:mIIa
rate=k(37)*y(50);
dy([50,51]) = dy([50,51]) + [-rate;rate];
% Xa:Va:mIIa -> Xa:Va + IIa
rate=k(38)*y(51);
dy([51,43,47]) = dy([51,43,47]) + [-rate;rate;rate];
% Release lipid from the above reaction since IIa is not lipid bound
lipidChange = lipidChange - rate;
% VII + Xa <-> VII:Xa
rate=(k(39)*y(10)*y(24)-k(40)*y(65));
dy([10,24,65]) = dy([10,24,65]) + [-rate;-rate;rate];
% VII:Xa -> VIIa + Xa
rate=k(41)*y(65);
dy([65,12,24]) = dy([65,12,24]) + [-rate;rate;rate];
% XI + IIa <-> XI:IIa
rate=(k(42)*y(52)*y(47)-k(43)*y(53));
dy([47,52,53]) = dy([47,52,53]) + [-rate;-rate;rate];
% XI:IIa -> XIa + IIa
rate=k(44)*y(53);
dy([53,47,54]) = dy([53,47,54]) + [-rate;rate;rate];
% PCa:PS + VIIIa <-> PCa:PS:VIIIa
rate=(k(45)*y(55)*y(16)-k(46)*y(56));
dy([16,55,56]) = dy([16,55,56]) + [-rate;-rate;rate];
% PCa:PS:VIIIa -> PCa:PS + VIIIai
rate=k(47)*y(56);
dy([56,30,55]) = dy([56,30,55]) + [-rate;rate;rate];
% PCa:PS + Va <-> PCa:PS:Va
rate=(k(48)*y(55)*y(8)-k(49)*y(62));
dy([8,55,62]) = dy([8,55,62]) + [-rate;-rate;rate];
% PCa:PS:Va -> PCa:PS + Vai
rate=k(50)*y(62);
dy([62,32,55]) = dy([62,32,55]) + [-rate;rate;rate];
% TFPI + Xa <-> Xa:TFPI
rate=(k(51)*y(57)*y(23)-k(52)*y(60));
dy([23,57,60]) = dy([23,57,60]) + [-rate;-rate;rate];
% Xa:TFPI + TF:VIIa <-> Xa:TFPI:TF:VIIa
rate=(k(53)*y(60)*y(36)-k(54)*y(61));
dy([36,60,61]) = dy([36,60,61]) + [-rate;-rate;rate];
% IXa + AT -> IXa:AT
rate=k(55)*y(19)*y(58);
dy([19,58,63]) = dy([19,58,63]) + [-rate;-rate;rate];
% Xa + AT -> Xa:AT
rate=k(56)*y(23)*y(58);
dy([23,58,64]) = dy([23,58,64]) + [-rate;-rate;rate];
% IIa + AT -> IIa:AT
rate=k(57)*y(47)*y(58);
dy([47,58,59]) = dy([47,58,59]) + [-rate;-rate;rate];
% V + mIIa <-> V:mIIa
rate=(k(58)*y(6)*y(4)-k(59)*y(66));
dy([4,6,66]) = dy([4,6,66]) + [-rate;-rate;rate];
% V:mIIa -> Va + mIIa
rate=k(60)*y(66);
dy([66,4,8]) = dy([66,4,8]) + [-rate;rate;rate];
% VIII + mIIa <-> VIII:mIIa
rate=(k(61)*y(14)*y(4)-k(62)*y(67));
dy([4,14,67]) = dy([4,14,67]) + [-rate;-rate;rate];
% VIII:mIIa -> VIIIa + mIIa
rate=k(63)*y(67);
dy([67,4,16]) = dy([67,4,16]) + [-rate;rate;rate];
% IIa + TM <-> IIa:TM
rate=(k(64)*y(47)*y(68)-k(65)*y(69));
dy([47,68,69]) = dy([47,68,69]) + [-rate;-rate;rate];
% IIa:TM + PC <-> IIa:TM:PC
rate=(k(66)*y(69)*y(34)-k(67)*y(70));
dy([34,69,70]) = dy([34,69,70]) + [-rate;-rate;rate];
% IIa:TM:PC -> IIa:TM + PCa
rate=k(68)*y(70);
dy([70,69,26]) = dy([70,69,26]) + [-rate;rate;rate];
% mIIa + AT -> mIIa:AT
rate=k(69)*y(3)*y(58);
dy([3,58,71]) = dy([3,58,71]) + [-rate;-rate;rate];
% PCa + PS <-> PCa:PS
rate=(k(70)*y(26)*y(28)-k(71)*y(55));
dy([26,28,55]) = dy([26,28,55]) + [-rate;-rate;rate];
% XIa + IX <-> XIa:IX
rate=(k(72)*y(54)*y(18)-k(73)*y(72));
dy([18,54,72]) = dy([18,54,72]) + [-rate;-rate;rate];
% XIa:IX -> XIa + IXa
rate=k(74)*y(72);
dy([72,20,54]) = dy([72,20,54]) + [-rate;rate;rate];

dy(73) = dy(73) - nva*lipidChange;


end
