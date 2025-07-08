load('parameters.mat', 'kParam')

c0 = setIC([15e-12, 1.4e-6, 2e-8, 1e-8, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9, 8.73e-6]);

[t,thr,OD,~]=unifiedODESolver(kParam,c0,1200,1e-14);

figure
plot(t/60,OD/10^-6,'r')
xlabel('Time (mins)')
ylabel('Active Substrate Concentration (\muM)')
title('Unified Model: Thrombin Generation Curve')
a = gca();
a.Box = "off";
a.XTick = 0:5:20;
y = ylim;
ylim([y(1),250])
xlim([0,20])
a.YTick = 0:50:250;
x_width = 3.42; y_width=2.57;
f = gcf();
f.Units = "inches";
set(gcf, 'Position', [0 0 x_width y_width]);
a.YGrid = 'on';
a.TitleFontWeight = "normal";
a.Box = 'off';
