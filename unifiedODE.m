% 1 TF
% 2 VII
% 3 VIIa
% 4 TF:VII
% 5 TF:VIIa
% 6 X
% 7 Xa
% 8 II
% 9 IIa
% 10 V
% 11 Va
% 12 VIII
% 13 VIIIa
% 14 XI
% 15 XIa
% 16 IX
% 17 IXa
% 18 IXa:VIIIa
% 19 Xa:Va
% 20 TFPI
% 21 Xa:TFPI
% 22 TF:VIIa:Xa:TFPI
% 23 ATIII
% 24 Xa:ATIII
% 25 IXa:ATIII
% 26 IIa:ATIII
% 27 TF:VIIa:ATIII
% 28 XIa:AT
% 29 VII:Xa
% 30 TFVII:Xa
% 31 VII:IIa
% 32 TFVII:IIa
% 33 VII:IXa
% 34 TFVII:IXa
% 35 VII:TFVIIa
% 36 X:TFVIIa
% 37 X:VIIa
% 38 X:IXaVIIIa
% 39 X:IXa
% 40 V:IIa
% 41 V:Xa
% 42 II:Xa
% 43 II:XaVa
% 44 XI:IIa
% 45 XI:XIa
% 46 IX:TFVIIa
% 47 IX:VIIa
% 48 IX:XIa
% 49 VIII:IIa
% 50 VIII:Xa
% 51 Substrate
% 52 Active Substrate
% 53 IIa:Substrate
% 54 a1AT
% 55 IIa:a1AT
% 56 Xa:a1AT
% 57 XIa:a1AT
% 58 a2AP
% 59 XIa:a2AP
% 60 a2M
% 61 IIa:a2M
% 62 C1-inh
% 63 XIa:C1-inh
% 64 PAI1
% 65 XIa:PAI1
% 66 VIIIa1L/VIIIa2
% 67 Xa:Va:AT
% 68 IXa:VIIIa:AT
% 69 Xa:Va:a1AT
% 70 Fbg
% 71 Fbg:IIa
% 72 Fbn1
% 73 Fbn1:IIa
% 74 Fbn2
% 75 Fbn2:IIa
% 76 Fbn1:IIa:AT
% 77 Fbn2:IIa:AT
% 78 FPA
% 79 FPB

function dy = unifiedODE(y,k)
%#codegen
dy = zeros(size(y));
    % TF + VII <-> TF:VII
    rate = k(1)*y(1)*y(2)-k(2)*y(4);
    dy([1,2,4]) = dy([1,2,4]) + [-rate;-rate;rate];
    
    % TF + VIIa <-> TF:VIIa
    rate = k(3)*y(1)*y(3)-k(4)*y(5);
    dy([1,3,5]) = dy([1,3,5]) + [-rate;-rate;rate];
    
    % VII + Xa <-> VII:Xa
    rate = k(5)*y(2)*y(7)-k(6)*y(29);
    dy([2,7,29]) = dy([2,7,29]) + [-rate;-rate;rate];
    % VII:Xa -> VIIa + Xa
    rate = k(7)*y(29);
    dy([29,3,7]) = dy([29,3,7]) + [-rate;rate;rate];
    
    % TFVII + Xa <-> TFVII:Xa
    rate = k(8)*y(4)*y(7)-k(9)*y(30);
    dy([4,7,30]) = dy([4,7,30]) + [-rate;-rate;rate];
    % TFVII:Xa -> TFVIIa + Xa
    rate = k(10)*y(30);
    dy([30,5,7]) = dy([30,5,7]) + [-rate;rate;rate];
    
    % VII + IIa <-> VII:IIa
    rate = k(11)*y(2)*y(9)-k(12)*y(31);
    dy([2,9,31]) = dy([2,9,31]) + [-rate;-rate;rate];
    % VII:IIa -> VIIa + IIa
    rate = k(13)*y(31);
    dy([31,3,9]) = dy([31,3,9]) + [-rate;rate;rate];
    
    % TFVII + IIa <-> TFVII:IIa
    rate = k(14)*y(4)*y(9)-k(15)*y(32);
    dy([4,9,32]) = dy([4,9,32]) + [-rate;-rate;rate];
    % TFVII:IIa -> TFVIIa + IIa
    rate = k(16)*y(32);
    dy([32,5,9]) = dy([32,5,9]) + [-rate;rate;rate];
    
    % VII + IXa <-> VII:IXa
    rate = k(17)*y(2)*y(17)-k(18)*y(33);
    dy([2,17,33]) = dy([2,17,33]) + [-rate;-rate;rate];
    % VII:IXa -> VIIa + IXa
    rate = k(19)*y(33);
    dy([33,3,17]) = dy([33,3,17]) + [-rate;rate;rate];
    
    % TFVII + IXa <-> TFVII:IXa
    rate = k(20)*y(4)*y(17)-k(21)*y(34);
    dy([4,17,34]) = dy([4,17,34]) + [-rate;-rate;rate];
    % TFVII:IXa -> TFVIIa + IXa
    rate = k(22)*y(34);
    dy([34,5,17]) = dy([34,5,17]) + [-rate;rate;rate];
    
    % VII + TFVIIa <-> VII:TFVIIa
    rate = k(23)*y(2)*y(5)-k(24)*y(35);
    dy([2,5,35]) = dy([2,5,35]) + [-rate;-rate;rate];
    % VII:TFVIIa -> VIIa + TFVIIa
    rate = k(25)*y(35);
    dy([35,3,5]) = dy([35,3,5]) + [-rate;rate;rate];
    
    % TFVIIa + AT -> TFVIIaAT
    rate = k(26)*y(5)*y(23);
    dy([5,23,27]) = dy([5,23,27]) + [-rate;-rate;rate];
    
    % X + TFVIIa <-> X:TFVIIa
    rate = k(27)*y(6)*y(5)-k(28)*y(36);
    dy([6,5,36]) = dy([6,5,36]) + [-rate;-rate;rate];
    % X:TFVIIa -> Xa + TFVIIa
    rate = k(29)*y(36);
    dy([36,7,5]) = dy([36,7,5]) + [-rate;rate;rate];
    
    % X + VIIa <-> X:VIIa
    rate = k(30)*y(6)*y(3)-k(31)*y(37);
    dy([6,3,37]) = dy([6,3,37]) + [-rate;-rate;rate];
    % X:VIIa -> Xa + VIIa
    rate = k(32)*y(37);
    dy([37,7,3]) = dy([37,7,3]) + [-rate;rate;rate];
    
    % X + IXaVIIIa <-> X:IXaVIIIa
    rate = k(33)*y(6)*y(18)-k(34)*y(38);
    dy([6,18,38]) = dy([6,18,38]) + [-rate;-rate;rate];
    % X:IXaVIIIa -> Xa + IXaVIIIa
    rate = k(35)*y(38);
    dy([38,7,18]) = dy([38,7,18]) + [-rate;rate;rate];
    
    % X + IXa <-> X:IXa
    rate = k(36)*y(6)*y(17)-k(37)*y(39);
    dy([6,17,39]) = dy([6,17,39]) + [-rate;-rate;rate];
    % X:IXa -> Xa + IXa
    rate = k(38)*y(39);
    dy([39,7,17]) = dy([39,7,17]) + [-rate;rate;rate];
    
    % V + IIa <-> V:IIa
    rate = k(39)*y(10)*y(9)-k(40)*y(40);
    dy([10,9,40]) = dy([10,9,40]) + [-rate;-rate;rate];
    % V:IIa -> Va + IIa
    rate = k(41)*y(40);
    dy([40,11,9]) = dy([40,11,9]) + [-rate;rate;rate];
    
    % V + Xa <-> V:Xa
    rate = k(42)*y(10)*y(7)-k(43)*y(41);
    dy([10,7,41]) = dy([10,7,41]) + [-rate;-rate;rate];
    % V:Xa -> Va + Xa
    rate = k(44)*y(41);
    dy([41,11,7]) = dy([41,11,7]) + [-rate;rate;rate];
    
    % Xa + Va <-> XaVa
    rate = k(45)*y(7)*y(11)-k(46)*y(19);
    dy([7,11,19]) = dy([7,11,19]) + [-rate;-rate;rate];
    
    % Xa + AT -> XaAT
    rate = k(47)*y(7)*y(23);
    dy([7,23,24]) = dy([7,23,24]) + [-rate;-rate;rate];
    
    % XaVa + AT -> XaVaAT
    rate = k(48)*y(19)*y(23);
    dy([19,23,67]) = dy([19,23,67]) + [-rate;-rate;rate];
    
    % II + Xa <-> II:Xa
    rate = k(49)*y(8)*y(7)-k(50)*y(42);
    dy([8,7,42]) = dy([8,7,42]) + [-rate;-rate;rate];
    % II:Xa -> IIa + Xa
    rate = k(51)*y(42);
    dy([42,9,7]) = dy([42,9,7]) + [-rate;rate;rate];
    
    % II + XaVa <-> II:XaVa
    rate = k(52)*y(8)*y(19)-k(53)*y(43);
    dy([8,19,43]) = dy([8,19,43]) + [-rate;-rate;rate];
    % II:XaVa -> IIa + XaVa
    rate = k(54)*y(43);
    dy([43,9,19]) = dy([43,9,19]) + [-rate;rate;rate];
    
    % IIa + AT -> IIaAT
    rate = k(55)*y(9)*y(23);
    dy([9,23,26]) = dy([9,23,26]) + [-rate;-rate;rate];
    
    % XI + IIa <-> XI:IIa
    rate = k(56)*y(14)*y(9)-k(57)*y(44);
    dy([14,9,44]) = dy([14,9,44]) + [-rate;-rate;rate];
    % XI:IIa -> XIa + IIa
    rate = k(58)*y(44);
    dy([44,15,9]) = dy([44,15,9]) + [-rate;rate;rate];
    
    % XI + XIa <-> XI:XIa
    rate = k(59)*y(14)*y(15)-k(60)*y(45);
    dy([14,15,45]) = dy([14,15,45]) + [-rate;-rate;rate];
    % XI:XIa -> XIa + XIa
    rate = k(61)*y(45);
    dy([45,15]) = dy([45,15]) + [-rate;2*rate];
    
    % XIa + AT -> XIaAT
    rate = k(62)*y(15)*y(23);
    dy([15,23,28]) = dy([15,23,28]) + [-rate;-rate;rate];
    
    % IX + TFVIIa <-> IX:TFVIIa
    rate = k(63)*y(16)*y(5)-k(64)*y(46);
    dy([16,5,46]) = dy([16,5,46]) + [-rate;-rate;rate];
    % IX:TFVIIa -> IXa + TFVIIa
    rate = k(65)*y(46);
    dy([46,17,5]) = dy([46,17,5]) + [-rate;rate;rate];
    
    % IX + VIIa <-> IX:VIIa
    rate = k(66)*y(16)*y(3)-k(67)*y(47);
    dy([16,3,47]) = dy([16,3,47]) + [-rate;-rate;rate];
    % IX:VIIa -> IXa + VIIa
    rate = k(68)*y(47);
    dy([47,17,3]) = dy([47,17,3]) + [-rate;rate;rate];
    
    % IX + XIa <-> IX:XIa
    rate = k(69)*y(16)*y(15)-k(70)*y(48);
    dy([16,15,48]) = dy([16,15,48]) + [-rate;-rate;rate];
    % IX:XIa -> IXa + XIa
    rate = k(71)*y(48);
    dy([48,17,15]) = dy([48,17,15]) + [-rate;rate;rate];
    
    % VIII + IIa <-> VIII:IIa
    rate = k(72)*y(12)*y(9)-k(73)*y(49);
    dy([12,9,49]) = dy([12,9,49]) + [-rate;-rate;rate];
    % VIII:IIa -> VIIIa + IIa
    rate = k(74)*y(49);
    dy([49,13,9]) = dy([49,13,9]) + [-rate;rate;rate];
    
    % VIII + Xa <-> VIII:Xa
    rate = k(75)*y(12)*y(7)-k(76)*y(50);
    dy([12,7,50]) = dy([12,7,50]) + [-rate;-rate;rate];
    % VIII:Xa -> VIIIa + Xa
    rate = k(77)*y(50);
    dy([50,13,7]) = dy([50,13,7]) + [-rate;rate;rate];
    
    % VIIIa <-> VIIIa1L + VIIIa2
    rate = k(78)*y(13) - k(79)*y(66)*y(66);
    dy([13,66]) = dy([13,66]) + [-rate;rate];
    % IXa:VIIIa -> IXa + VIIIa1L + VIIIa2
    rate = k(80)*y(18);
    dy([18,17,66]) = dy([18,17,66]) + [-rate;rate;rate];
    
    % IXa + VIIIa <-> IXaVIIIa
    rate = k(81)*y(17)*y(13)-k(82)*y(18);
    dy([17,13,18]) = dy([17,13,18]) + [-rate;-rate;rate];
    
    % IXa + AT -> IXaAT
    rate = k(83)*y(17)*y(23);
    dy([17,23,25]) = dy([17,23,25]) + [-rate;-rate;rate];
    
    % IXaVIIIa + AT -> IXaVIIIaAT
    rate = k(84)*y(18)*y(23);
    dy([18,23,68]) = dy([18,23,68]) + [-rate;-rate;rate];
    
    % Xa + TFPI <-> Xa:TFPI
    rate = k(85)*y(7)*y(20)-k(86)*y(21);
    dy([7,20,21]) = dy([7,20,21]) + [-rate;-rate;rate];
    
    % TFVIIa + XaTFPI <-> TFVIIaXaTFPI
    rate = k(87)*y(21)*y(5)-k(88)*y(22);
    dy([21,5,22]) = dy([21,5,22]) + [-rate;-rate;rate];
    
    % IIa + a1AT -> IIa:a1AT
    rate = k(89)*y(9)*y(54);
    dy([9,54,55]) = dy([9,54,55]) + [-rate;-rate;rate];
    % Xa + a1AT -> Xa:a1AT
    rate = k(90)*y(7)*y(54);
    dy([7,54,56]) = dy([7,54,56]) + [-rate;-rate;rate];
    % Xa:Va + a1AT -> Xa:Va:a1AT
    rate = k(91)*y(19)*y(54);
    dy([19,54,69]) = dy([19,54,69]) + [-rate;-rate;rate];
    % XIa + a1AT -> XIa:a1AT
    rate = k(92)*y(15)*y(54);
    dy([15,54,57]) = dy([15,54,57]) + [-rate;-rate;rate];
    % XIa + a2AP -> XIa:a2AP
    rate = k(93)*y(15)*y(58);
    dy([15,58,59]) = dy([15,58,59]) + [-rate;-rate;rate];
    % IIa + a2M -> IIa:a2M
    rate = k(94)*y(9)*y(60);
    dy([9,60,61]) = dy([9,60,61]) + [-rate;-rate;rate];
    % XIa + C1inh -> XIa:C1inh
    rate = k(95)*y(15)*y(62);
    dy([15,62,63]) = dy([15,62,63]) + [-rate;-rate;rate];
    % XIa + PAI1 -> XIa:PAI1
    rate = k(96)*y(15)*y(64);
    dy([15,64,65]) = dy([15,64,65]) + [-rate;-rate;rate];
    
    % Fbg + IIa <-> Fbg:IIa
    rate = k(97)*y(70)*y(9) - k(98)*y(71);
    dy([70,9,71]) = dy([70,9,71]) + [-rate;-rate;rate];
    % Fbg:IIa -> Fbn1 + IIa + FPA
    rate = k(99)*y(71);
    dy([71,72,9,78]) = dy([71,72,9,78]) + [-rate;rate;rate;rate];
    
    % Fbn1 + IIa <-> Fbn1:IIa
    rate = k(100)*y(72)*y(9) - k(101)*y(73);
    dy([72,9,73]) = dy([72,9,73]) + [-rate;-rate;rate];
    
    % Fbn1:IIa -> Fbn2 + IIa + FPB
    rate = k(102)*y(73);
    dy([73,74,9,79]) = dy([73,74,9,79]) + [-rate;rate;rate;rate];
    
    % Fbn2 + IIa <-> Fbn2:IIa
    rate = k(103)*y(74)*y(9) - k(104)*y(75);
    dy([74,9,75]) = dy([74,9,75]) + [-rate;-rate;rate];
    
    % Fbn1:IIa + AT -> Fbn1:IIa:AT
    rate = k(105)*y(73)*y(23);
    dy([73,23,76]) = dy([73,23,76]) + [-rate;-rate;rate];
    
    % Fbn2:IIa + AT -> Fbn2:IIa:AT
    rate = k(106)*y(75)*y(23);
    dy([75,23,77]) = dy([75,23,77]) + [-rate;-rate;rate];
    
    % IIa + Substrate <-> IIa:Substrate
    rate = k(107)*y(9)*y(51)-k(108)*y(53);
    dy([9,51,53]) = dy([9,51,53]) + [-rate;-rate;rate];
    % IIa:Substrate -> IIa + ActiveSubstrate
    rate = k(109)*y(53);
    dy([53,9,52]) = dy([53,9,52]) + [-rate;rate;rate];
end
