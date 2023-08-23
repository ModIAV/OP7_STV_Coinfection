%%% --------------------------------------------------------------- %%%
%%% Evaluation, Optimization and Simulation of a mathematical model %%%
%%% of OP7 and STV coinfection                                      %%%
%%% --------------------------------------------------------------- %%%
%
%   authors: Daniel Ruediger
%   last revised: 2023/06/13

%% Housekeeping commands
close all; clear all; clc;

global p d

% define compiler (example, if not already integrated)
% setenv('MW_MINGW64_LOC','C:\...\mingw64');

% check for IQM toolbox (formerly known as SB Toolbox 2)
if ( ~exist('IQMsimulate','file') )
    error(['Please install the IQM toolbox as it is required for model simulation!\n', ...
                   'Check: https://iqmtools.intiquan.com/']);
end 

% check for fSSm algorithm
ChangeSearchAlgorithm = 0;
if ( ~exist('fssm_kernel','file') )
    ChangeSearchAlgorithm = 1;
	fprintf('fSSm algortihm not found. Switching to "fmincon".\n\n') 
end

tmp = dir('MexModelFile_Opt.mexa64'); % delete MEX files to prevent occasional crashes
for i = 1:length(tmp);  delete(tmp(i).name);  end

%% Options
p.ModelName          = 'OP7_STV_ModelFile.txt'; % name of model file
p.OptimizeParameters = 0; % enable parameter optimization
    %%% Fitting strategy
    p.FittingStrategy = 1; % 1 - only fit RKI data
                           % 2 - fit OP7 data w/ estimated baseline parameters

p.LogData    = 1; % use log data                              
p.maxtime    = 30; % model simualtion time
p.MOI        = 10; % applied MOI

p.ParticleRatioIn_COI_override = []; % test other COI input scenario ([S7-STV, S1-6/S8, S7-OP7])
                                     % [] = use initial values from data

% plotting/save options
p.SaveFigures = 0;
p.PNGFileName = '';
               
%% Model structure options
% degradation of infectious virions
p.kDegV = 0.1; 

%%% Advantages %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.ModOpt.Combined_cvRNA_Adv = 1; % advantage in vRNA transcription
    p.OptChange.Fadv_cvRNA = 0.412724481920294; % 'ratio' - calculate Fadv from vRNA length ratio // #number# - apply value for Fadv

p.ModOpt.mRNA_Adv = 1; % reduction of mRNA transcription
    p.OptChange.Fadv_mRNA = -0.9312; % #number# - apply value for mRNA reduction 
    
p.ModOpt.Lv9  = 'mut';   % 'mut' - same length as normal S7 // 'del' - use given length (nt)
    p.OptChange.Lv9  = 1027;

p.ModOpt.Lm9  = 'mut';   % 'mut' - same length as normal S7 // 'del' - use given length (nt)
    p.OptChange.Lm9  = 1005;
    
%%% Different M1 + M2 function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
p.ModOpt.OP7_M1_DiffBind = 1; % OP7-derived M1 binds with different rate to vRNA in the nucleus than normal M1
    p.OptChange.Fbind = 0; % define value change for binding rate: 1 = normal, 0 = no binding, # = use factor
    
% not used in final model    
p.ModOpt.OP7_M1_noRel = 0; % OP7-derived M1 does not contribute to virion release

p.ModOpt.OP7_M2_noRel = 0; % OP7-derived M2 does not contribute to virion release    
    
p.ModOpt.OP7_M1_DiffNucExp = 0; % OP7-derived M1 leads to different nuclear export rate of vRNP than normal M1
    p.OptChange.Fexp = 0; % define value change for nuclear export rate: 1 = normal, 0 = no binding, # = use factor    0
    
p.ModOpt.OP7_M1_DiffM19Syn = 0; % OP7-derived M1 is synthesized at a different rate than normal M1
    p.OptChange.Fsyn = 0; % define value change for synthesis rate: 1 = normal, 0 = no binding, # = use factor        
    
p.ModOpt.OP7_M1_DiffVpDeg = 0; % vRNPs covered by OP7-derived M1 are degraded at a different rate than normal M1
    p.OptChange.FdegVp = 0; % define value change for degradation rate: 1 = normal, 0 = no binding, # = use factor            

%% Other options
p.SimTime = p.maxtime*1.05;        %h  maximum simulation time
   
p.ChangedParaNames  = {'kSynV', 'kSynC',  'kSynM', 'kSynP', 'kRel', 'KvRel', 'kBindM1', 'kBindNp'};   % names  of parameters that ought to be changed in the model
p.ChangedParaValues = [ 1.01e-1, 1.98e-4,  2.68e5,  3.81e4,  5.05e2, 2.03e0,  5.75e-6,   1.37e-6];   % values of parameters that ought to be changed in the model            
p.FPar_STV = 0.199;     

if ( p.OptimizeParameters )
    fprintf('Use final parameters with random noise for estimation.\n\n');
    
    if ( p.FittingStrategy == 1 )
        p.ChangedParaValues = p.ChangedParaValues + randn(size(p.ChangedParaValues))*0.1.*p.ChangedParaValues;
        p.FPar_STV          = randn*0.1*p.FPar_STV;
    else
        p.OptChange.Fadv_cvRNA = p.OptChange.Fadv_cvRNA + randn*0.1*p.OptChange.Fadv_cvRNA;
        p.OptChange.Fadv_mRNA  = p.OptChange.Fadv_mRNA  + randn*0.1*p.OptChange.Fadv_mRNA;
    end
end

p.FPar_COI = p.FPar_STV;

p.OdeOptions.abstol = 1e-30;    % absolute tolerance of the ODE solver
p.OdeOptions.reltol = 1e-10;    % relative tolerance of the ODE solver

%% Optimization options
p.MaxTime            = 1*60*60; %maximum run time of optimization solver in seconds
p.LogPara            = 0; % use logarithmic parameters for estimation (converted internally)
                       
if ( p.FittingStrategy == 1 )
    p.PaOpt        = {'kSynV', ...
                      'kSynC', ...
                      'kSynM', ... 
                      'kSynP', ...
                      'kRel', ...
                      'KvRel', ... 
                      'kBindM1', ...
                      'kBindNp', ...
                      }; %parameters subject to optimization in the intracellular model          
    p.PaOpt_Titers = {'FPar_STV'}; %parameters subject to optimization applied to virus titers  
    p.IcOpt        = {};    %initial conditions subject to optimization in the intracellular model, set Vex as the first variable

    p.States2FitSTV = {'RvSeg5', 'RvSeg7', 'RvSeg8', ... % fit time course of states
                       'Rm5',    'Rm7',    'Rm8', ...
                       'RcSeg5', 'RcSeg7', 'RcSeg8',...
                       'TCID', 'HA', ... 
                       'PB1tot', 'NPtot', 'M1_WTtot', ...
                       }; 
                   
    p.NucPercentageSTV   = 0;      
    p.Composition2FitSTV = 0;     % fit virus titers determined at start or end

    p.States2FitCOI = {};
    p.Composition2FitCOI = 0;
    p.NucPercentageCOI   = 0; 
    
    p.Composition2FitPassages = 0;    % fit virus titers determined at 12 hpi    

elseif ( p.FittingStrategy == 2 )
    p.PaOpt = {'Fadv_cvRNA', ...
               'Fadv_mRNA', ...
               }; %parameters subject to optimization in the intracellular model          
    p.PaOpt_Titers = {... 'FPar_COI'
                      }; %parameters subject to optimization applied to virus titers           
    p.IcOpt = {};    %initial conditions subject to optimization in the intracellular model, set Vex as the first variable

    p.States2FitSTV      = {};
    p.Composition2FitSTV = 0;     % fit virus titers determined at start or end
    p.NucPercentageSTV   = 0; 

    p.States2FitCOI = {'RvSeg5', 'RvSeg7', 'RvSeg9', 'RvSeg8', ... % fit time course of states
                       'Rm5',    'Rm7',    'Rm9',    'Rm8', ...
                       'RcSeg5', 'RcSeg7', 'RcSeg9', 'RcSeg8',...
                       'HA', 'TCID', ...
                       'PB1tot', 'NPtot', 'M1_WTtot', 'M1_OP7tot', ...
                       }; 
                   
    p.Composition2FitCOI = 0;    % fit virus titers determined at 12 hpi
    p.NucPercentageCOI   = 0;     
    
    p.Composition2FitPassages = 0;    % fit virus titers determined at 12 hpi    
end

p.RelStateWeightSTV = ones(size(p.States2FitSTV));     %relative weight of the deviation of each fitted state in the objective function
p.RelStateWeightSTV2 = ones(size(p.Composition2FitSTV));     %relative weight of the deviation of each fitted state in the objective function
    
p.RelStateWeightCOI = ones(size(p.States2FitCOI));     %relative weight of the deviation of each fitted state in the objective function
p.RelStateWeightCOI2 = ones(size(p.Composition2FitCOI));     %relative weight of the deviation of each fitted state in the objective function
    
p.nPaOpt = length(p.PaOpt);                          %number of parameters subject to optimization
p.nOpt   = length(p.PaOpt)+length(p.IcOpt);       %number of parameters and initial conditions subject to optimization 

p.PaDev = 100;                 %upper and lower bounds of parameters deviate by */ PaDev from initial parameter guess
p.IcDev = 10;                 %upper and lower bounds of IC deviate by */ IcDev from IC

% Optimization minor options
if ( ChangeSearchAlgorithm )
    p.Solver  = 'fmincon';             %optimization solver: 'fminsearch', 'fmincon', 'fSSm'
else
    p.Solver  = 'fSSm';             
end
p.MaxIter = 1e10; % e8                %maximum iterations
p.MaxEval = 1e9; % e7                %maximum function evaluations
p.SearchRadiusConfidence = []; 0.2; % parameter for CMA-ES algorithm, defines search range (confidence radius), not required when upper + lower bounds are defined

p.SmallRefset    = 1; % reduce Refset size for faster (but less robust) estimation   

p.ObjFunCal = 'log';      %'normal': model and experiment are compared quantitatively 
                                       %'Log': model states and experimental values are logarithmized prior to comparison
                                       %'log-extra_0_rule': simulation number that go to 0 are handled differently, calculate percentage of ExpVal to max(ExpVal)
                                       %'NormToMax': model and experimental values are normalized to maximum prior to comparison
                                       %'NormToState': model and experimental values are normalized to a specific state given by p.NormState prior to comparison
                                       %'NormToInitial': model and experimental values are normalized to inital values of each state prior to comparison
                                       %'NormToInitialState: model and experimental values are normalized to initial value of a specific state given by p.NormState prior to comparison                              
p.NormState = '';            %state variables and experiment values are normalized to the inital value of this state prior to their comparison when p.ObjFunCal = 'NormToInitialState'

%% Load experimental data
load('ExpData_d.mat') 

%% Load, adjust and compile model

% check which paras exist in the basic model
Model = IQMmodel(fullfile(strcat(pwd, filesep), p.ModelName));
[p.ParaNames, ~] = IQMparameters(Model);

tmp = boolean(sum(strcmp(repmat(p.ChangedParaNames,length(p.ParaNames),1),...
                   repmat(p.ParaNames,1,length(p.ChangedParaNames))),1));
               
p.CPN_pre  = p.ChangedParaNames(tmp);
p.CPV_pre  = p.ChangedParaValues(tmp);
p.CPN_post = p.ChangedParaNames(~tmp);
p.CPV_post = p.ChangedParaValues(~tmp);

% update parameters + load model
Model = IQMparameters(Model, p.CPN_pre, p.CPV_pre);
p.VariableNames             =  IQMvariables(Model);
[p.StateNames, ~, p.Ic]     =     IQMstates(Model);
[p.ParaNames, p.ParaValues] = IQMparameters(Model);

MS = IQMstruct(Model);

%%% Lv9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( strcmp('mut',p.ModOpt.Lv9) )
    MS.parameters(strcmp('Lv9', p.ParaNames)).value = MS.parameters(strcmp('Lv7', p.ParaNames)).value;
elseif ( strcmp('del',p.ModOpt.Lv9) )
    MS.parameters(strcmp('Lv9', p.ParaNames)).value = p.OptChange.Lv9;
else
    error('Check options!');
end
%%% advantage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( p.ModOpt.Combined_cvRNA_Adv )
    % reaction
    tmp = find(strcmp('rSynRc9', p.VariableNames));
    MS.variables(tmp).formula = '(Fadv_cvRNA+1)*kSynC*VpNuc9*P_Rdrp';    
    
    tmp = find(strcmp('rSynRv9', p.VariableNames));
    MS.variables(tmp).formula = '(Fadv_cvRNA+1)*kSynV*Cp9*P_Rdrp';     
    
    % parameter
    MS.parameters(end+1).name  = 'Fadv_cvRNA';
    if ( strcmp('ratio',p.OptChange.Fadv_cvRNA) )
        MS.parameters(end).value   = 0.5 ...
                                     * (MS.parameters(strcmp('Lv7', p.ParaNames)).value ...
                                     / MS.parameters(strcmp('Lv9', p.ParaNames)).value ...
                                     - 1);
    else
        MS.parameters(end).value = p.OptChange.Fadv_cvRNA;
    end  
end

if ( p.ModOpt.mRNA_Adv )
    % reaction
    tmp = find(strcmp('rSynRm9', p.VariableNames));
    MS.variables(tmp).formula = '(Fadv_mRNA+1)*kSynM/Lm9*VpNuc9';
    
    % parameter
    MS.parameters(end+1).name  = 'Fadv_mRNA';
    if ( strcmp('ratio',p.OptChange.Fadv_mRNA) )
        MS.parameters(end).value   = MS.parameters(strcmp('Lv7', p.ParaNames)).value ...
                                     / MS.parameters(strcmp('Lv9', p.ParaNames)).value ...
                                     - 1;
    else
        MS.parameters(end).value = p.OptChange.Fadv_mRNA;
    end  
end

%%% OP7-derived M1 binding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( p.ModOpt.OP7_M1_DiffBind ) % binding factor of OP7-derived OP7M1 to vRNA   
    % reaction
    tmp = find(strcmp('rBindM19_1', p.VariableNames));
    for i = 1 : 9
        MS.variables(tmp+i-1).formula = sprintf('Fbind*kBindM1*P_M19*VpNuc%i',i);      
    end
    
    % parameter
    MS.parameters(end+1).name = 'Fbind';
    MS.parameters(end).value  = p.OptChange.Fbind;    
end
%%% OP7-derived M1 virion release %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( p.ModOpt.OP7_M1_noRel ) 
    % ODE
    tmp = find(strcmp('P_M19', p.StateNames));
    MS.states(tmp).ODE = 'rSynM19 - rBindM19';
    
    tmp = find(strcmp('M1_OP7tot', p.StateNames));
    MS.states(tmp).ODE = 'rSynM19';    
    
    % variables
    tmp = find(strcmp('P_M1tot', p.VariableNames));
    MS.variables(tmp).formula = 'P_M17';        
end
%%% OP7-derived M2 virion release %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( p.ModOpt.OP7_M2_noRel ) 
    % ODE
    tmp = find(strcmp('P_M29', p.StateNames));
    MS.states(tmp).ODE = 'rSynM29';    
    
    % variables
    tmp = find(strcmp('P_M2tot', p.VariableNames));
    MS.variables(tmp).formula = 'P_M27';        
end
%%% OP7-derived M1 nuclear export %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( p.ModOpt.OP7_M1_DiffNucExp ) % nuclear export factor of vRNPs bound by OP7-derived OP7M1
    % reaction
    tmp = find(strcmp('rExp9_1', p.VariableNames));
    for i = 1 : 9
        MS.variables(tmp+i-1).formula = sprintf('Fexp*kExp*P_Nep*VpNucM19_%i',i);      
    end
    
    % parameter
    MS.parameters(end+1).name = 'Fexp';
    MS.parameters(end).value  = p.OptChange.Fexp;      
end
%%% OP7-derived M1 synthesis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( p.ModOpt.OP7_M1_DiffM19Syn ) % transcription factor of OP7-derived OP7M1
    % reaction
    tmp = find(strcmp('rSynM19', p.VariableNames));
    MS.variables(tmp).formula = sprintf('Fsyn*kSynP/Drib*(1-Fspl7)*Rm9');      
    
    % parameter
    MS.parameters(end+1).name = 'Fsyn';
    MS.parameters(end).value  = p.OptChange.Fsyn;     
end
%%% degradation of vRNP covered in OP7-derived M1 %%%%%%%%%%%%%%%%%%%%%%%%%
if ( p.ModOpt.OP7_M1_DiffVpDeg ) % degradation factor of vRNPs bound by OP7-derived OP7M1
    % reaction
    tmp = find(strcmp('rDegVpNucM19_1', p.VariableNames));
    for i = 1 : 9
        MS.variables(tmp+i-1).formula = sprintf('FdegVp*kDegRnp*VpNucM19_%i',i);      
    end
    tmp = find(strcmp('rDegVpCytM19_1', p.VariableNames));
    for i = 1 : 9
        MS.variables(tmp+i-1).formula = sprintf('FdegVp*kDegRnp*VpCytM19_%i',i);      
    end    
    
    % parameter
    MS.parameters(end+1).name = 'FdegVp';
    MS.parameters(end).value  = p.OptChange.FdegVp;        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model = IQMmodel(MS);

Model = IQMparameters(Model, p.CPN_post, p.CPV_post);
p.VariableNames             =  IQMvariables(Model);
[p.StateNames, ~, p.Ic]     =     IQMstates(Model);
[p.ParaNames, p.ParaValues] = IQMparameters(Model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    IQMmakeMEXmodel(Model, 'MexModelFile_Opt');
    p.CompileFlag = 1;  % compilation successful
catch ModelCompileError
    p.CompileFlag = 0;  % compilation failed
    
    p.Model = Model;
    fprintf(['\nWarning: Model could not be compiled into a MEX-file.\n', ...
         '         Simulation will take longer than using a compiler and\n', ...
         '         may generate additional warnings during simulation.\n', ...
         '         The usage of a compiler is advised.\n\n']);
end

%% Set initial conditions
% extract input from experiments
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

p.ParticleRatioIn_STV    = [1, 1, 0];
p.ParticleRatioIn_COI    = 10.^[d.S7_WT.LogSeed(1), ...
                                mean([d.S5.LogSeed(1), d.S8.LogSeed(1)]),...
                                d.S7_OP7.LogSeed(1)];  
                            
p.InfectingParticles_STV = p.MOI;
p.InfectingParticles_COI = p.MOI;                            
        
%% Model optimization
if ( p.OptimizeParameters )
    p.Ic = zeros(size(p.Ic));
    
    Opt = Opt_Optimization;

    % convert logarithmic paras
    if ( p.LogPara )
        p.Opt.result = 10.^Opt.xbest;
    else
        p.Opt.result = Opt.xbest;
    end

    % CMAES does not privide best result at end?
    if ( p.Opt.result ~= p.safety.xbest )
        warning('Best parameter set changed to actual optimal result!');
        p.Opt.result = p.safety.xbest;
    end    

    % display best para set
    tmp = [p.PaOpt, p.PaOpt_Titers];
    fprintf('\n');
    for i = 1 : length(tmp)
        fprintf('%s\t%.3e\n',tmp{i},p.Opt.result(i));
    end
    fprintf('\n');
    
    if ( p.FittingStrategy == 1 )
        p.FPar_STV = p.Opt.result(end);
        p.FPar_COI = p.FPar_STV;
    end

    % save results
    Pnames        = p.PaOpt;
    Pnames_Titers = p.PaOpt_Titers;
    Pvalues       = p.Opt.result;
    ObjFunVal     = Opt.ObjFunVal;
    
else
    p.PaOpt      = {};
    p.nPaOpt     = 0;
    p.Opt.result = [];
end

%% simulate model
p.Ic = zeros(size(p.Ic));

% override
if ( ~isempty(p.ParticleRatioIn_COI_override) )
    p.ParticleRatioIn_COI    = p.ParticleRatioIn_COI_override;
end

% infection with standard viruses and no DIP or OP7
p.Ic(strcmp('Vex',p.StateNames))  = p.InfectingParticles_STV;
p.Ic(strcmp('RpII',p.StateNames)) = 1;
if(p.CompileFlag)
    resultSTV = IQMPsimulate('MexModelFile_Opt', p.SimTime, p.Ic, ...
                  [p.PaOpt, {'pInLowSeg', 'pInMedSeg', 'pInHighSeg'}], ...
                  [p.Opt.result(1:p.nPaOpt), p.ParticleRatioIn_STV], ...
                  p.OdeOptions);              
else        
    tmp = Model;
    tmp = IQMparameters(tmp, [p.PaOpt, {'pInLowSeg', 'pInMedSeg', 'pInHighSeg'}], ...
                             [p.Opt.result(1:p.nPaOpt), p.ParticleRatioIn_STV]);            
    
    resultSTV = IQMsimulate(tmp, 'ode23s', p.SimTime, p.Ic, p.OdeOptions);
    resultSTV.time = resultSTV.time';
end

% infection with one standard virus and one defective interfering particle
p.Ic(strcmp('Vex',p.StateNames))  = p.InfectingParticles_COI;
p.Ic(strcmp('RpII',p.StateNames)) = 1;
if(p.CompileFlag)
    resultCOI = IQMPsimulate('MexModelFile_Opt', p.SimTime, p.Ic, ...
                  [p.PaOpt, {'pInLowSeg', 'pInMedSeg', 'pInHighSeg'}], ...
                  [p.Opt.result(1:p.nPaOpt), p.ParticleRatioIn_COI], ...
                  p.OdeOptions);
else
    
    tmp = Model;
    tmp = IQMparameters(tmp, [p.PaOpt, {'pInLowSeg', 'pInMedSeg', 'pInHighSeg'}], ...
                             [p.Opt.result(1:p.nPaOpt), p.ParticleRatioIn_COI]);            
    
    resultCOI = IQMsimulate(tmp, 'ode23s', p.SimTime, p.Ic, p.OdeOptions);
    resultCOI.time = resultCOI.time';
end

AllTime.STV       = resultSTV.time;
AllProperties.STV = [resultSTV.states, resultSTV.variables];
AllValues.STV     = [resultSTV.statevalues, resultSTV.variablevalues];

AllTime.COI       = resultCOI.time;
AllProperties.COI = [resultCOI.states, resultCOI.variables];
AllValues.COI     = [resultCOI.statevalues, resultCOI.variablevalues];

AllValues.STV(AllValues.STV<0) = 0;
AllValues.COI(AllValues.COI<0) = 0;

%% Plot figures
OP7_VisualizeResults(AllTime,AllValues,AllProperties,resultSTV,resultCOI)














































