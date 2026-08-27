%Setup
load('priors.mat')
load('parameters.mat')
iterCount = 2000;
timepoints = 0:0.1:20;

priorSTD(priorSTD==2.5) = 0.5;
priorSTD(75:84) = 0.5;

c0 = setIC([15e-12, 1.4e-6, 2e-8, 1e-8, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9, 8.73e-6]);
ODArray = zeros(length(timepoints),iterCount);
thrArray = ODArray;
for i = 1:iterCount
    k = sample(priorMean,priorSTD);
     [tmpOD,tmpThr] = simulate(k,c0,timepoints);
     ODArray(:,i) = tmpOD;
end
save('densityPlotOverPriorData.mat','ODArray','thrArray')
densityPlot(true)

%% Density over posterior
load('posterior.mat')
load('parameters.mat')
c0 = setIC([15e-12, 1.4e-6, 2e-8, 1e-8, 7e-10, 9e-8, 1.6e-7, 3e-8, 3.4e-6, 2.5e-9, 8.73e-6]);
ratesAndCounts = cell(200,2);
k=1;
for i = 1:2000
    alreadyExists = any(cellfun(@(x) isequal(x, nextPop(i,:)), ratesAndCounts));
    if ~alreadyExists
        ratesAndCounts{k,1} = nextPop(i,:);
        count = 0;
        for j = 1:2000
            if all(isequal(nextPop(i,:),nextPop(j,:)))
                count = count + 1;
            end
        end
        ratesAndCounts{k,2} = count;
        k = k + 1;
    end
end

ki = 1:84;
timepoints = 0:0.1:20;
ODArray = zeros(length(timepoints),2000);
k = 1;
for i = 1:length(ratesAndCounts)
    rates = ratesAndCounts{i,1};
    kPrior(ki) = rates;
    [OD,~] = simulate(prior2param(kPrior),c0,timepoints);
    for j = 1:ratesAndCounts{i,2}
        ODArray(:,k) = OD;
        k = k + 1;
    end
end

save('densityPlotOverPosteriorData.mat')
densityPlot(false)

%% Functions

function [OD,thr] = simulate(k,c0,timepoints)
    [~,~,~,sol] = unifiedODESolver(k,c0,timepoints(end)*60,1e-14);
    OD = deval(sol,timepoints*60,52)/1e-6;
    thr = (deval(sol,timepoints*60,9)+deval(sol,timepoints*60,53))/1e-9;
end

function k = sample(mean,std)
    kPrior = 10.^(normrnd(mean,std));
    k = prior2param(kPrior);
end

function densityPlot(prior)
    if prior
        load('densityPlotOverPriorData.mat')
    else
        load('densityPlotOverPosteriorData.mat')
    end
    %Percentiles from curves
    data = ODArray;
    percentilePoints = [5,25,50,75,95];
    timepoints = 0:0.1:20;
    percentiles = zeros(length(timepoints),length(percentilePoints));
    for i = 1:length(timepoints)
        percentiles(i,:) = prctile(data(i,:),percentilePoints);
    end
    
    %Plot curves
    figure
    hold on
    numOfBoxes = (length(percentilePoints)-1)/2;
    for i = 1:numOfBoxes
    upperBoundary = percentiles(:,i)';
    lowerBoundary = percentiles(:,length(percentilePoints)-i+1)';
    patch([timepoints fliplr(timepoints)], [upperBoundary  fliplr(lowerBoundary)], [0.85, 1, 0.85]/i);
    end
    
    plot(timepoints,percentiles(:,3),'k','LineWidth',1)
    legend({'90% Interval','50% Interval','Median'},'Location','northwest')
    if prior
        title('Density Over Priors')
    else
        title('Density Over Posteriors')
    end
    xlabel('Time (mins)')
    ylabel('Active Substrate Conc. (\muM)')
    
    x_width = 3.42; y_width=2.57;
    a = gca();
    f = gcf();
    f.Units = "inches";
    set(gcf, 'Position', [0 0 x_width y_width]);
    a.YGrid = 'on';
    a.TitleFontWeight = "normal";
    a.Box = 'off';
    if prior
        exportgraphics(f,'prior.png','Resolution',300)
    else
        exportgraphics(f,'posterior.png','Resolution',300)
    end
    if prior
        % Extra figure with just 3 curves
        figure
        hold on
        for i=1:5
            plot(timepoints, data(:,i), 'k', 'LineWidth',1)
        end
        title('Prior sampled OD curves')
        xlabel('Time (mins)')
        ylabel('Active Substrate Conc. (\muM)')
        a = gca();
        f = gcf();
        f.Units = "inches";
        set(gcf, 'Position', [0 0 x_width y_width]);
        a.YGrid = 'on';
        a.TitleFontWeight = "normal";
        a.Box = 'off';
        
        exportgraphics(f,'prior-samples.png','Resolution',300)
    end
end
