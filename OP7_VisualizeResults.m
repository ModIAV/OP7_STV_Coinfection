%%% ----------------------------------------------------------- %%%
%%% Visualization of simulation results of a mathematical model %%%
%%% of OP7 and STV coinfection                                  %%%
%%% ----------------------------------------------------------- %%%
%
%   authors: Daniel Ruediger
%   last revised: 2023/06/13

function OP7_VisualizeResults(AllTime,AllValues,AllProperties,resultSTV,resultCOI)

global p d

%% Visualize simulation results - RNAs up tp 12 hpi
close all;
colors = [0   0.3  0.6; % S5
          0   0    0;   % S7_WT
          0.9 0.42 0;   % S7_OP7
          0   0.5  0];  % S8
LW = 1.5; % LineWidth
MS = 5;   % MarkerSize
symbols   = {'d', 'x', 'o', 's'};
lineStyle = {'-', '-.', '-', '--'};

ylims     = [1e-1 7e6;
             2e-3 2e4;
             3e0 3e6];         

InfStyle = {'STV', 'COI'};
InfName  = {'RKI', 'OP7'};
RNAs     = {'mRNA', 'cRNA', 'vRNA'};
 
Segs     = {'S5', 'S7_WT', 'S7_OP7', 'S8'};
ModelRNA = {'Rm5', 'Rm7', 'Rm9', 'Rm8'; ...
            'RcSeg5', 'RcSeg7', 'RcSeg9', 'RcSeg8'; ...
            'RvSeg5', 'RvSeg7', 'RvSeg9', 'RvSeg8'};
    
h.Fig = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [1 2 14 20]);

for i = 1 : 2 % RKI / OP7
    for j = 1 : 3 % RNA species

        % RNA levels in co-infected cell
        subplot(3,2,(j-1)*2+i)
        hold on;

        for k = 1 : length(Segs)

            % include offset
            tmp = AllValues.(InfStyle{i})(:,strcmp(ModelRNA{j,k},AllProperties.(InfStyle{i}))); 
            tmp = tmp + 10.^d.(RNAs{j}).(Segs{k}).(InfName{i})(1,1);

            h.(['p', num2str(k)]) = plot(AllTime.(InfStyle{i}), tmp, ...
                 lineStyle{k}, 'Color', colors(k,:), 'LineWidth', LW);

            if ( p.LogData )
                errorbar(d.time.RNA,10.^d.(RNAs{j}).(Segs{k}).(InfName{i})(:,1), ...
                         -(10.^(d.(RNAs{j}).(Segs{k}).(InfName{i})(:,1) - d.(RNAs{j}).(Segs{k}).(InfName{i})(:,2)) - 10.^d.(RNAs{j}).(Segs{k}).(InfName{i})(:,1)), ...
                           10.^(d.(RNAs{j}).(Segs{k}).(InfName{i})(:,1) + d.(RNAs{j}).(Segs{k}).(InfName{i})(:,2)) - 10.^d.(RNAs{j}).(Segs{k}).(InfName{i})(:,1), ...
                           symbols{k} ,'Color',colors(k,:),'LineWidth', LW)             
            else
                % avoid negative errorbars in log space
                tmp = d.(RNAs{j}).(Segs{k}).(InfName{i})(:,1) < d.(RNAs{j}).(Segs{k}).(InfName{i})(:,2);
                NegError = 10.^d.(RNAs{j}).(Segs{k}).(InfName{i})(:,2);
                NegError(tmp) = 0.99999*10.^d.(RNAs{j}).(Segs{k}).(InfName{i})(tmp,1);

                errorbar(d.time.RNA,10.^d.(RNAs{j}).(Segs{k}).(InfName{i})(:,1), ...
                                    NegError, ...
                                    10.^d.(RNAs{j}).(Segs{k}).(InfName{i})(:,2), ...
                         symbols{k} ,'Color',colors(k,:),'LineWidth', LW, 'MarkerSize', MS)                        
            end
               
            if ( i == 1 ) 
                try
                   plot(d.time.Total(10:end),10.^d.RKI.ESV(10:end,strcmp(ModelRNA{j,k},d.ExpStateNames)), ...
                        symbols{k},'Color',colors(k,:),'LineWidth', LW, 'MarkerSize', MS);
                catch
                end
            end

        end

        set(gca, 'FontSize', 9, 'FontName', 'Arial',...
                 'XLim', [-0.5 12.5],  'XTick', 0:6:p.maxtime,...
                 'YLim', ylims(j,:), 'YTick', 10.^(-6:2:10), 'YScale', 'log'); %,...

        xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
        ylabel(sprintf('log_{10}(%s/cell)',RNAs{j}), 'FontName', 'Arial', 'FontSize', 10);

        if ( j == 1 )
            title(InfName{i}, 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');

            if ( i == 1 )     
                legend([h.p1, h.p2, h.p3, h.p4], 'S5', 'S7_W_T', 'S7_O_P_7', 'S8', 'Location', 'SE')
            end
        end

        box on;

    end
end

if ( p.SaveFigures )
    export_fig(sprintf('figures/%s_Fig1%i',p.PNGFileName), '-r300', '-png', '-painters', '-nocrop');
end

%% Visualize simulation results - virus titers

h.Fig = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [16 17 16 6]);

% TCID50 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,1)
hold on;

k = 1;
TCrange = 2:length(d.time.Titers);
if ( p.LogData )                
    errorbar(d.time.Titers(TCrange),10.^d.TCID.dyn.(InfName{k})(TCrange,1), ...
             -(10.^(d.TCID.dyn.(InfName{k})(TCrange,1) - d.TCID.dyn.(InfName{k})(TCrange,2)) - 10.^d.TCID.dyn.(InfName{k})(TCrange,1)), ...
               10.^(d.TCID.dyn.(InfName{k})(TCrange,1) + d.TCID.dyn.(InfName{k})(TCrange,2)) - 10.^d.TCID.dyn.(InfName{k})(TCrange,1), ...
               symbols{5*(k-1)+1} ,'Color',colors(k+1,:),'LineWidth', LW)    
else
    % avoid negative errorbars in log space
    tmp = d.TCID.dyn.RKI(:,1) < d.TCID.dyn.RKI(:,2);
    NegError = 10.^d.TCID.dyn.RKI(:,2);
    NegError(tmp) = 0.99999*10.^d.TCID.dyn.RKI(tmp,1);

    errorbar(d.time.Titers,10.^d.TCID.dyn.RKI(:,1), ...
                           NegError, ...
                           10.^d.TCID.dyn.RKI(:,2), ...
                           symbols{1} ,'Color', 'k','LineWidth', LW, 'MarkerSize', MS)            
end
            
k = 2;
TCrange = setdiff(2:length(d.time.Titers), 6);                   
if ( p.LogData )                
    errorbar(d.time.Titers(TCrange),10.^d.TCID.dyn.(InfName{k})(TCrange,1), ...
             -(10.^(d.TCID.dyn.(InfName{k})(TCrange,1) - d.TCID.dyn.(InfName{k})(TCrange,2)) - 10.^d.TCID.dyn.(InfName{k})(TCrange,1)), ...
               10.^(d.TCID.dyn.(InfName{k})(TCrange,1) + d.TCID.dyn.(InfName{k})(TCrange,2)) - 10.^d.TCID.dyn.(InfName{k})(TCrange,1), ...
               'x' ,'Color',colors(k+1,:),'LineWidth', LW)    
else
    % avoid negative errorbars in log space
    tmp = d.TCID.dyn.OP7(:,1) < d.TCID.dyn.OP7(:,2);
    NegError = 10.^d.TCID.dyn.OP7(:,2);
    NegError(tmp) = 0.99999*10.^d.TCID.dyn.OP7(tmp,1);

    errorbar(d.time.Titers,10.^d.TCID.dyn.OP7(:,1), ...
                           NegError, ...
                           10.^d.TCID.dyn.OP7(:,2), ...
                           'x' ,'Color', colors(k,:),'LineWidth', LW, 'MarkerSize', MS)                        
end
                   
% include Fpar and decay over time for TCID   
tmp = d.NumInfectedCells*AllValues.STV(:,strcmp('TCID',AllProperties.STV));
tmp = tmp .* p.FPar_STV .* exp(-p.kDegV*AllTime.STV');
plot(AllTime.STV, tmp, 'Color', 'k', 'LineWidth', LW);   

tmp = d.NumInfectedCells*AllValues.COI(:,strcmp('TCID',AllProperties.COI));
tmp = tmp .* p.FPar_COI .* exp(-p.kDegV*AllTime.COI');
plot(AllTime.COI, tmp, 'Color', colors(k+1,:), 'LineWidth', LW);   
    
set(gca, 'FontSize', 9, 'FontName', 'Arial',...
         'XLim', [0 32],  'XTick', 0:6:32,...
         'YLim', [1e1 1e9], 'YTick', 10.^(2:2:12), 'YScale', 'log');

xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel(sprintf('log_{10}(infectious\nvirions/mL)'), 'FontName', 'Arial', 'FontSize', 10);     
title('TCID_5_0', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');     

box on;

% HA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(1,2,2)
hold on;

k = 1;
TCrange = 2:length(d.time.Titers);
if ( p.LogData )
    tmp = 10.^(d.HA.dyn.(InfName{k})(TCrange,1) + d.HA.dyn.(InfName{k})(TCrange,2)) - 10.^d.HA.dyn.(InfName{k})(TCrange,1);

    a = find(d.HA.dyn.(InfName{k})(TCrange,2) > d.HA.dyn.(InfName{k})(TCrange,1));
    tmp(a) = 10.^d.HA.dyn.(InfName{k})(TCrange(a),2);

    errorbar(d.time.Titers(TCrange),10.^d.HA.dyn.(InfName{k})(TCrange,1), ...
             -(10.^(d.HA.dyn.(InfName{k})(TCrange,1) - d.HA.dyn.(InfName{k})(TCrange,2)) - 10.^d.HA.dyn.(InfName{k})(TCrange,1)), ...
               tmp, ...
               symbols{5*(k-1)+1} ,'Color',colors(k+1,:),'LineWidth', LW)  
else
    % avoid negative errorbars in log space
    tmp = d.HA.dyn.RKI(:,1) < d.HA.dyn.RKI(:,2);
    NegError = 10.^d.HA.dyn.RKI(:,2);
    NegError(tmp) = 0.99999*10.^d.HA.dyn.RKI(tmp,1);

    errorbar(d.time.Titers,10.^d.HA.dyn.RKI(:,1), ...
                           NegError, ...
                           10.^d.HA.dyn.RKI(:,2), ...
                           symbols{1} ,'Color', 'k','LineWidth', LW, 'MarkerSize', MS)                 
end

k = 2;
TCrange = setdiff(2:length(d.time.Titers), 6);
if ( p.LogData )
    tmp = 10.^(d.HA.dyn.(InfName{k})(TCrange,1) + d.HA.dyn.(InfName{k})(TCrange,2)) - 10.^d.HA.dyn.(InfName{k})(TCrange,1);

    a = find(d.HA.dyn.(InfName{k})(TCrange,2) > d.HA.dyn.(InfName{k})(TCrange,1));
    tmp(a) = 10.^d.HA.dyn.(InfName{k})(TCrange(a),2);

    errorbar(d.time.Titers(TCrange),10.^d.HA.dyn.(InfName{k})(TCrange,1), ...
             -(10.^(d.HA.dyn.(InfName{k})(TCrange,1) - d.HA.dyn.(InfName{k})(TCrange,2)) - 10.^d.HA.dyn.(InfName{k})(TCrange,1)), ...
               tmp, ...
               'x' ,'Color',colors(k+1,:),'LineWidth', LW)  
else         
    % avoid negative errorbars in log space
    tmp = d.HA.dyn.OP7(:,1) < d.HA.dyn.OP7(:,2);
    NegError = 10.^d.HA.dyn.OP7(:,2);
    NegError(tmp) = 0.99999*10.^d.HA.dyn.OP7(tmp,1);

    errorbar(d.time.Titers,10.^d.HA.dyn.OP7(:,1), ...
                           NegError, ...
                           10.^d.HA.dyn.OP7(:,2), ...
                           symbols{2} ,'Color', 'r','LineWidth', LW, 'MarkerSize', MS)                                    
end
                   
tmp = d.NumInfectedCells*AllValues.STV(:,strcmp('HA',AllProperties.STV));
plot(AllTime.STV, tmp, 'Color', 'k', 'LineWidth', LW);              

tmp = d.NumInfectedCells*AllValues.COI(:,strcmp('HA',AllProperties.COI));           
plot(AllTime.COI, tmp, 'Color', colors(k+1,:), 'LineWidth', LW);  
     

set(gca, 'FontSize', 9, 'FontName', 'Arial',...
         'XLim', [0 32],  'XTick', 0:6:32,...
         'YLim', [1e7 1e11], 'YTick', 10.^(0:2:16), 'YScale', 'log')

xlabel('time post infection (h)', 'FontName', 'Arial', 'FontSize', 10);
ylabel(sprintf('log_{10}(total\nvirions/mL)'), 'FontName', 'Arial', 'FontSize', 10);     
title('HA', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');     

box on;

if ( p.SaveFigures )
    export_fig(sprintf('figures/%s_Fig2',p.PNGFileName), '-r300', '-png', '-painters', '-nocrop');
end

%% Visualize simulation results -  release composition - normalized new
% close all;

XPos = [0.15, 0.75, 0.55]; % XStart, XWidth, XJump
YPos = [0.1, 0.35, 0.5]; 

YLims = [1e-4 1.5e1];

% extract results from experiments
ParHA_SeedRKI = 2*10^7*10^d.HA.Seed.RKI(1,1);
ParHA_SeedOP7 = 2*10^7*10^d.HA.Seed.OP7(1,1);
ParHA_12hpiRKI = 2*10^7*10^d.HA.hpi12.RKI(1,1);
ParHA_12hpiOP7 = 2*10^7*10^d.HA.hpi12.OP7(1,1);

NormSeed = [d.S5.Seed.RKI(1,1), d.S7_WT.Seed.RKI(1,1), ...
            NaN, d.S8.Seed.RKI(1,1)] ./ ParHA_SeedRKI;
Norm12hpi = [d.S5.hpi12.RKI(1,1), d.S7_WT.hpi12.RKI(1,1), ...
            NaN, d.S8.hpi12.RKI(1,1)] ./ ParHA_12hpiRKI;  
NormSeed(3) = NormSeed(2);        
Norm12hpi(3) = Norm12hpi(2);       

% normalize to 8 segments per particle
SegPerParSeed = p.ParticleRatioIn_COI;

% calculate simulation results
t_12hpi = find(resultCOI.time>12,1,'first');

Segs = [5 7 8 9];
SegPerPar12hpi = zeros(1,length(Segs));
for i = 1 : length(Segs)
        
    % Vrel
    try
        if ( sum(Segs(i) == [5 7 8]) )
            SegPerPar12hpi(i) = SegPerPar12hpi(i) + AllValues.COI(t_12hpi, strcmp('Vrel',AllProperties.COI));
        end
    catch
    end
    
    % Drel
    try
        if ( sum(Segs(i) == [5 8 9]) )
            SegPerPar12hpi(i) = SegPerPar12hpi(i) + AllValues.COI(t_12hpi, strcmp('Drel',AllProperties.COI));
        end    
    catch
    end
    
    % OP7UP + OP7 RNG
    if ( Segs(i) == 9 )
        SegPerPar12hpi(i) = SegPerPar12hpi(i) + AllValues.COI(t_12hpi, strcmp('OP7rel',AllProperties.COI));
    end    
    
    % SXrel
    try
        if ( sum(Segs(i) == [5 8]) )
            SegPerPar12hpi(i) = SegPerPar12hpi(i) + AllValues.COI(t_12hpi, strcmp('SXrel',AllProperties.COI));
        end    
    catch
    end    
    
    % S7rel
    try
        if ( Segs(i) == 7 )
            SegPerPar12hpi(i) = SegPerPar12hpi(i) + AllValues.COI(t_12hpi, strcmp('S7rel',AllProperties.COI));
        end    
    catch
    end       
    
end   
SegPerPar12hpi = SegPerPar12hpi/AllValues.COI(t_12hpi, strcmp('HA',AllProperties.COI));

h.Fig = figure('color', 'w',...
               'paperpositionmode', 'auto',...
               'paperunits', 'centimeters',...
               'paperposition', [0 0 16 10],...
               'units', 'centimeters',...
               'position', [16 2 10 12]);

%%% OP7 segment seed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes; hold on;

plot([-1 11], [1 1], 'k--', 'LineWidth', 1);

bar(1, 10^d.S5.LogSeed(1), 'FaceColor', colors(1,:))
bar(2, 10^d.S7_WT.LogSeed(1), 'FaceColor', colors(2,:))
bar(3, 10^d.S7_OP7.LogSeed(1), 'FaceColor', colors(3,:))
bar(4, 10^d.S8.LogSeed(1), 'FaceColor', colors(4,:))

plot([5 5], [eps 1e20], 'k-', 'LineWidth', 1.5);

bar(6, SegPerParSeed (2), 'FaceColor', colors(1,:))
bar(7, SegPerParSeed (1), 'FaceColor', colors(2,:))
bar(8, SegPerParSeed (3), 'FaceColor', colors(3,:))
bar(9, SegPerParSeed (2), 'FaceColor', colors(4,:))

set(gca, 'FontSize', 9, 'FontName', 'Arial',...
         'XLim', [0 10],  'XTick', [1:4, 6:9], 'XTickLabel', {'5', 'WT7', 'OP7', '8'}, ...
         'YLim', YLims, 'YTick', 10.^(-4:1:20),...
         'YScale', 'log',...
         'Position', [XPos(1)+0.05, YPos(1) + YPos(3), XPos(2), YPos(2)-0.02]); %,...
     
xlabel('segment', 'FontName', 'Arial', 'FontSize', 10);
ylabel(sprintf('seed: log_{10}(average\nvRNAs/virion)'), 'FontName', 'Arial', 'FontSize', 10);

X = get(gca,'XLim'); Y = get(gca,'YLim');
text(0.2*X(2),2*Y(2),'data', 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');
text(0.7*X(2),2*Y(2),'sim', 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');

box on;

%%% OP7 segment hpi 12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes; hold on;

plot([-1 11], [1 1], 'k--', 'LineWidth', 1);

bar(1, 10^d.S5.LogHpi12(1), 'FaceColor', colors(1,:))
bar(2, 10^d.S7_WT.LogHpi12(1), 'FaceColor', colors(2,:))
bar(3, 10^d.S7_OP7.LogHpi12(1), 'FaceColor', colors(3,:))
bar(4, 10^d.S8.LogHpi12(1), 'FaceColor', colors(4,:))

plot([5 5], [eps 1e20], 'k-', 'LineWidth', 1.5);

bar(6, SegPerPar12hpi(1), 'FaceColor', colors(1,:))
bar(7, SegPerPar12hpi(2), 'FaceColor', colors(2,:))
bar(8, SegPerPar12hpi(4), 'FaceColor', colors(3,:))
bar(9, SegPerPar12hpi(3), 'FaceColor', colors(4,:))

set(gca, 'FontSize', 9, 'FontName', 'Arial',...
         'XLim', [0 10],  'XTick', [1:4, 6:9], 'XTickLabel', {'5', 'WT7', 'OP7', '8'}, ...
         'YLim', YLims, 'YTick', 10.^(-4:1:20),...
         'YScale', 'log',...
         'Position', [XPos(1)+0.05, YPos(1), XPos(2), YPos(2)-0.02]); %,...
     
xlabel('segment', 'FontName', 'Arial', 'FontSize', 10);
ylabel(sprintf('12 hpi: log_{10}(average\nvRNAs/virion)'), 'FontName', 'Arial', 'FontSize', 10);

X = get(gca,'XLim'); Y = get(gca,'YLim');
text(0.2*X(2),2*Y(2),'data', 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');
text(0.7*X(2),2*Y(2),'sim', 'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');

box on;

if ( p.SaveFigures )
    export_fig(sprintf('figures/%s_Fig3',p.PNGFileName), '-r300', '-png', '-painters', '-nocrop');
end






















































