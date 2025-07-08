% TFVIIa
% TFVII
% TF
% VIIa
% VII
% IXa
% IX
% Xa
% X
% IIa
% II
% Ia
% I
% VIIIa
% VIII
% Va
% V
% XIa
% XI
% AT
% TFPI
% XaTFPI
% APC
% PC
function [t,thr,sol]=main(k,c0,maxt)
% k - Parameters
% c0 - Initial concentrations
options = odeset('AbsTol',1e-16);
sol = ode23tb(@(t,y)f(y,k),[0,maxt/60],c0/10^-9,options);
t = sol.x*60;
y = sol.y*10^-9;
thr = y(10,:);
end


function dy = f(y,k)
p = k(73);
i = k(66:72);
h = k(42:65);
n = k(38:41);
K = k(20:37);
km = k(17:19);
k = k(1:16);

dy = zeros(size(y));
% TFVIIaF = y(1)/(1+y(7)/K(1) + y(9)/K(3));
% XaTFVIIa = k(6)/(K(3)*km(19))*y(9)*TFVIIaF;
% IXaBF = y(6)*p*n(20)/(K(20)+y(6));
%Concentrations
VIIaTF = y(1);
VIITF = y(2);
TF = y(3);
VIIa = y(4);
VII = y(5);
IXa = y(6);
IX = y(7);
Xa = y(8);
X = y(9);
IIa = y(10);
II = y(11);
Ia = y(12);
I = y(13);
VIIIa = y(14);
VIII = y(15);
Va = y(16);
V = y(17);
XIa = y(18);
XI = y(19);
AT = y(20);
TFPI = y(21);
XaTFPI = y(22);
PCa = y(23);
PC = y(24);
%Indexes
iVIIaTF = 1;
iVIITF = 2;
iTF = 3;
iVIIa = 4;
iVII = 5;
iIXa = 6;
iIX = 7;
iXa = 8;
iX = 9;
iIIa = 10;
iII = 11;
iIa = 12;
iI = 13;
iVIIIa = 14;
iVIII = 15;
iVa = 16;
iV = 17;
iXIa = 18;
iXI = 19;
iAT = 20;
iTFPI = 21;
iXaTFPI = 22;
iPCa = 23;
iPC = 24;
     IXaBF = IXa*p*n(1)/(K(11)+IXa);
     VIIaTFF = VIIaTF/(1+IX/K(1)+X/K(3));
     XaVIIaTF = k(6)/(K(3)*km(3))*X*VIIaTFF;
     VaB = Va*p*n(4)/(K(18)+Va);
     XaVaB = Xa*VaB/(K(14)*(1+i(7)/K(15)+Xa/K(14))+VaB);
     VaBF = VaB-XaVaB;
     XaF = Xa-XaVaB;
     IIB = II*p*n(3)/(K(17)*(1+X/K(16)+II/K(17)));
     IIaF = IIa/(1+(Ia+I)/K(8));
     XB = X*p*n(3)/(K(16)*(1+X/K(16)+II/K(17)));
     VIIIaBF = VIIIa*p*n(2)/((K(12)+VIIIa)*(1+XB/(p*K(7))*(1+i(7)/K(13))));
     

     dy(1) = k(1)*VIIa*TF - km(1)*VIIaTFF + k(2)*VIITF*IIaF + k(3)*VIITF*XaF - h(1)*VIIaTFF*XaTFPI - h(2)*XaVIIaTF*TFPI;
     dy(2) = k(1)*VII*TF - km(1)*VIITF - k(2)*VIITF*IIaF - k(3)*VIITF*XaF;
     dy(3) =  - (k(1)*VIIa*TF - km(1)*VIIaTFF) - (k(1)*VII*TF - km(1)*VIITF);
     dy(4) =  - (k(1)*VIIa*TF - km(1)*VIIaTFF) + k(2)*VII*IIaF;
     dy(5) =  - (k(1)*VII*TF - km(1)*VIITF) - k(2)*VII*IIaF;
     dy(6) = k(4)/K(1)*IX*VIIaTFF + k(5)*IX*XIa/(K(2) + IX) - h(3)*AT*IXa;
     dy(7) =  - (k(4)/K(1))*IX*VIIaTFF - k(5)*IX*XIa/(K(2) + IX);
     dy(8) = k(6)/K(3)*X*VIIaTFF + k(7)*IXaBF*XB/(p*K(4)) + k(8)*IXaBF*VIIIaBF*XB/(p^2*K(6)*K(5)) - (k(9)*XaF*TFPI - km(2)*XaTFPI) - (h(4)*AT + h(5)*i(1) + h(6)*i(2) + h(7)*i(5))*XaF - h(8)*AT*XaVaB;
     dy(9) =  - (k(6)/K(3))*X*VIIaTFF - k(7)*IXaBF*XB/(p*K(4)) - k(8)*IXaBF*VIIIaBF*XB/(p^2*K(6)*K(5));
     dy(10) = k(10)*p*XaF*II + k(11)*XaVaB*IIB/p - (h(9)*AT + h(10)*i(1) + h(11)*i(2) + h(12)*i(5) + h(13)*i(4))*IIaF;
     dy(11) =  - k(10)*p*XaF*II - k(11)*XaVaB*IIB/p;
     dy(12) = k(12)/K(8)*I*IIaF;
     dy(13) =  - (k(12)/K(8))*I*IIaF;
     dy(14) = k(13)*VIII*IIaF/(K(9) + IIaF) - h(14)*VIIIa;
     dy(15) =  - k(13)*VIII*IIaF/(K(9) + IIaF);
     dy(16) = k(14)*V*IIaF/(K(10) + IIaF) - h(15)*PCa*VaBF;
     dy(17) =  - k(14)*V*IIaF/(K(10) + IIaF);
     dy(18) = k(15)*p*XI*IIaF - (h(16)*AT + h(17)*i(3) + h(18)*i(2) + h(19)*i(5) + h(20)*i(6))*XIa;
     dy(19) =  - k(15)*p*XI*IIaF;
     dy(20) =  - (h(3)*IXa + h(4)*XaF + h(8)*XaVaB + h(9)*IIaF + h(16)*XIa)*AT;
     dy(21) =  - (k(9)*XaF*TFPI - km(2)*XaTFPI) - h(2)*XaVIIaTF*TFPI;
     dy(22) = k(9)*XaF*TFPI - km(2)*XaTFPI - h(1)*VIIaTFF*XaTFPI;
     dy(23) = k(16)*PC*IIaF - (h(21)*i(1) + h(22)*i(3) + h(23)*i(2) + h(24)*i(5))*PCa;
     dy(24) =  - k(16)*PC*IIaF;

     
%      dy = zeros(24,1);
%      %Mass Action Law Reactions
%      %TF + VII <-> TF:VII
%      rate = k(1)*TF*VII-km(1)*VIITF;
%      dy([iVII,iTF,iVIITF]) = dy([iVII,iTF,iVIITF]) + [-rate;-rate;rate];
%      %TF + VIIa <-> TF:VIIaF
%      rate = k(1)*TF*VIIa-km(1)*VIIaTFF;
%      dy([iVIIa,iTF,iVIIaTF]) = dy([iVIIa,iTF,iVIIaTF]) + [-rate;-rate;rate];
%      %TF:VII + IIaF -> TF:VIIa + IIaF
%      rate = k(2)*VIITF*IIaF;
%      dy([iVIITF,iVIIaTF]) = dy([iVIITF,iVIIaTF]) + [-rate;rate];
%      %TF:VII + XaF -> TF:VIIa + XaF
%      rate = k(3)*VIITF*XaF;
%      dy([iVIITF,iVIIaTF]) = dy([iVIITF,iVIIaTF]) + [-rate;rate];
%      %VII + IIaF -> VIIa + IIaF
%      rate = k(2)*VII*IIaF;
%      dy([iVII,iVIIa]) = dy([iVII,iVIIa]) + [-rate;rate];
%      %TF:VIIaF + Xa:TFPI -> TF:VIIa:Xa:TFPI
%      rate = h(1)*VIIaTFF*XaTFPI;
%      dy([iVIIaTF,iXaTFPI]) = dy([iVIIaTF,iXaTFPI]) + [-rate;-rate];
%      %TF:VIIa:Xa + TFPI -> TF:VIIa:Xa:TFPI
%      rate = h(2)*XaVIIaTF*TFPI;
%      dy([iVIIaTF,iTFPI]) = dy([iVIIaTF,iTFPI]) + [-rate;-rate];
%      %IX + TF:VIIaF -> IXa + TF:VIIa
%      rate = k(4)/K(1)*IX*VIIaTFF;
%      dy([iIX,iIXa]) = dy([iIX,iIXa]) + [-rate;rate];
%      %IXa + AT -> IXa:AT
%      rate = h(3)*IXa*AT;
%      dy([iIXa,iAT]) = dy([iIXa,iAT]) + [-rate;-rate];
%      %X + TF:VIIaF -> Xa + TF:VIIa
%      rate = k(6)/K(3)*X*VIIaTFF;
%      dy([iX,iXa]) = dy([iX,iXa]) + [-rate;rate];
%      %XaF + TFPI <-> Xa:TFPI
%      rate = k(9)*XaF*TFPI-km(2)*XaTFPI;
%      dy([iXa,iTFPI,iXaTFPI]) = dy([iXa,iTFPI,iXaTFPI]) + [-rate;-rate;rate];
%      %XaF + AT -> Xa:AT
%      rate = h(4)*XaF*AT;
%      dy([iXa,iAT]) = dy([iXa,iAT]) + [-rate;-rate];
%      %XaF + a2M -> Xai + a2M
%      rate = h(5)*XaF*i(1);
%      dy(iXa) = dy(iXa) -rate;
%      %XaF + a1AT -> Xai + a1AT
%      rate = h(6)*XaF*i(2);
%      dy(iXa) = dy(iXa) -rate;
%      %XaF + PCI -> Xai + PCI
%      rate = h(7)*XaF*i(5);
%      dy(iXa) = dy(iXa) -rate;
%      %XaVaB + AT -> Xa:AT + Va
%      rate = h(8)*XaVaB*AT;
%      dy([iXa,iAT]) = dy([iXa,iAT]) + [-rate;-rate];
%      %II + XaF -> IIa + Xa
%      rate = k(10)*p*II*XaF;
%      dy([iII,iIIa]) = dy([iII,iIIa]) + [-rate;rate];
%      %IIB + XaVaB -> IIa + XaVa
%      rate = k(11)/p*IIB*XaVaB;
%      dy([iII,iIIa]) = dy([iII,iIIa]) + [-rate;rate];
%      %IIaF + AT -> IIa:AT
%      rate = h(9)*IIaF*AT;
%      dy([iIIa,iAT]) = dy([iIIa,iAT]) + [-rate;-rate];
%      %IIaF + a2M -> IIai + a2M
%      rate = h(10)*IIaF*i(1);
%      dy(iIIa) = dy(iIIa) -rate;
%      %IIaF + a1AT -> IIai + a1AT
%      rate = h(11)*IIaF*i(2);
%      dy(iIIa) = dy(iIIa) -rate;
%      %IIaF + PCI -> IIai + PCI
%      rate = h(12)*IIaF*i(5);
%      dy(iIIa) = dy(iIIa) -rate;
%      %IIaF + hep -> IIai + hep
%      rate = h(13)*IIaF*i(4);
%      dy(iIIa) = dy(iIIa) -rate;
%      %I + IIaF -> Ia + IIa
%      rate = k(12)/K(8)*I*IIaF;
%      dy([iI,iIa]) = dy([iI,iIa]) + [-rate;rate];
%      %VIIIa -> VIIIai
%      rate = h(14)*VIIIa;
%      dy(iVIIIa) = dy(iVIIIa) -rate;
%      %VaBF + PCa -> Vai + PCa
%      rate = h(15)*VaBF*PCa;
%      dy(iVa) = dy(iVa) -rate;
%      %XI + IIaF -> XIa + IIa
%      rate = k(15)*p*XI*IIaF;
%      dy([iXI,iXIa]) = dy([iXI,iXIa]) + [-rate;rate];
%      %XIa + AT -> XIa:AT
%      rate = h(16)*XIa*AT;
%      dy([iXIa,iAT]) = dy([iXIa,iAT]) + [-rate;-rate];
%      %XIa + a2AP -> XIa:a2AP
%      rate = h(17)*XIa*i(3);
%      dy(iXIa) = dy(iXIa) -rate;
%      %XIa + a1AT -> XIa:a1AT
%      rate = h(18)*XIa*i(2);
%      dy(iXIa) = dy(iXIa) -rate;
%      %XIa + PCI -> XIai + PCI
%      rate = h(19)*XIa*i(5);
%      dy(iXIa) = dy(iXIa) -rate;
%      %XIa + C1inh -> XIai + C1inh
%      rate = h(20)*XIa*i(6);
%      dy(iXIa) = dy(iXIa) -rate;
%      %PC + IIaF -> PCa + IIa
%      rate = k(16)*PC*IIaF;
%      dy([iPC,iPCa]) = dy([iPC,iPCa]) + [-rate;rate];
%      %PCa + a2M -> PCai + a2M
%      rate = h(21)*PCa*i(1);
%      dy(iPCa) = dy(iPCa) -rate;
%      %PCa + a2AP -> PCai + a2AP
%      rate = h(22)*PCa*i(3);
%      dy(iPCa) = dy(iPCa) -rate;
%      %PCa + a1AT -> PCai + a1AT
%      rate = h(23)*PCa*i(2);
%      dy(iPCa) = dy(iPCa) -rate;
%      %PCa + PCI -> PCai + PCI
%      rate = h(24)*PCa*i(5);
%      dy(iPCa) = dy(iPCa) -rate;
%      
%      
%      %Michaelis Menten Reactions
%      %IX by XIa
%      rate = k(5)*IX*XIa/(K(2)+IX);
%      dy([iIX,iIXa]) = dy([iIX,iIXa]) + [-rate;rate];
%      
%      %XB by IXaBF
%      rate = k(7)*IXaBF*XB/(p*K(4));
%      dy([iX,iXa]) = dy([iX,iXa]) + [-rate;rate];
%           
%      %XB by IXaBF*VIIIaBF
%      rate = k(8)*IXaBF*VIIIaBF*XB/(p^2*K(5)*K(6));
%      dy([iX,iXa]) = dy([iX,iXa]) + [-rate;rate];
%           
%      %VIII by IIaF
%      rate = k(13)*IIaF*VIII/(K(9)+IIaF);
%      dy([iVIII,iVIIIa]) = dy([iVIII,iVIIIa]) + [-rate;rate];
%      
%      %V by IIaF
%      rate = k(14)*IIaF*V/(K(10)+IIaF);
%      dy([iV,iVa]) = dy([iV,iVa]) + [-rate;rate];

end


