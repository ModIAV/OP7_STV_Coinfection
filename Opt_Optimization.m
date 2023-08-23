%%% --------------------------------------------- %%%
%%% Initialization for the parameter optimization %%%
%%% for a model of OP7 and STV coinfection        %%%
%%% --------------------------------------------- %%%
%
%   authors: Daniel Ruediger
%   last revised: 2023/06/13

function result = Opt_Optimization

global p
global d
global stopOptimization

stopOptimization = 0;

%%
p.TotalOptTime = 0;
p.TotalTimeCounter = 1;

result = []; 

%% Optimization solver options                    

%%% fSSm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( p.LogPara )
    p.fSSmOptions.log_var  = []; % 1:p.nOpt+p.Ex.nOpt;
else
    p.fSSmOptions.log_var  = []; % 1:p.nOpt+p.Ex.nOpt;
end
p.fSSmOptions.maxtime      = p.MaxTime;            % maximal time limit in seconds (e.g. 1 day = 60*60*24)
p.fSSmOptions.maxeval      = p.MaxEval;            % maximal number of function evaluations (e.g. 1000000)
p.fSSmOptions.local.solver = 0;

if ( p.SmallRefset )
    nvar = p.nOpt;
    p.fSSmOptions.ndiverse     = 4*nvar; % 2 = arbitrary factor, normally 10, changed to 4 here
    tmp = roots([1, -1, -3*nvar]);       % 3 = arbitrary factor, normally 10
    p.fSSmOptions.dim_refset   = ceil(tmp(tmp>0));
    if ( mod(p.fSSmOptions.dim_refset,2) )
        p.fSSmOptions.dim_refset = p.fSSmOptions.dim_refset + 1; % dim_refset should be an even number
    end
end

%%% Fminsearch and Fmincon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.FminOptions = optimset('MaxIter',     p.MaxIter,...
                         'MaxFunEvals', p.MaxEval,...
                         'Display',     'iter');
                         
%% Prepare optimization of the intracellular model
% Find position of initial conditions subject to optimization in p.Ic
if(isempty(p.IcOpt))
    p.PosIcOpt = [];
else
    p.PosIcOpt = zeros(length(p.IcOpt),1);
    for CountIc=1:length(p.IcOpt)
        if(sum(strcmp(p.IcOpt(CountIc), p.StateNames))==0)
            disp(char(strcat('Warning(Optimization): "', p.IcOpt(CountIc), '" is no initial condition in the intracellular model.')));
            return
        else
            p.PosIcOpt(CountIc) = find(strcmp(p.IcOpt(CountIc), p.StateNames));
        end
    end
end

% Set optimization bounds and initial guess for the intracellular model
p.x_L = zeros(p.nOpt,1);
p.x_U = zeros(p.nOpt,1);
p.x0  = zeros(p.nOpt,1);
for CountOptPa  = 1 : length(p.PaOpt)    
    if(sum(strcmp(p.PaOpt(CountOptPa ), p.ParaNames))==0)
        disp(strcat('Warning(Optimization): "', p.PaOpt(CountOptPa ),...
            '" is no parameter in the intracellular model.'));
        return
    end
    p.x_L(CountOptPa) = p.ParaValues(strcmp(p.PaOpt(CountOptPa), p.ParaNames))/p.PaDev;
    p.x_U(CountOptPa) = p.ParaValues(strcmp(p.PaOpt(CountOptPa), p.ParaNames))*p.PaDev;
    p.x0(CountOptPa)  = p.ParaValues(strcmp(p.PaOpt(CountOptPa), p.ParaNames));
end

if(~isempty(p.IcOpt))
    for CountOptIc=1:length(p.IcOpt)
        p.x_L(length(p.PaOpt)+CountOptIc) = 0;
        p.x_U(length(p.PaOpt)+CountOptIc) = 0;
        p.x0(length(p.PaOpt)+CountOptIc)  = 0;
    end
end

% Find position of states subject to optimization
p.PosState2FitSTV = zeros(0,2);
for CountFit=1:length(p.States2FitSTV)
    if(sum(strcmp(p.States2FitSTV(CountFit), [p.StateNames; p.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.States2FitSTV(CountFit), '" is no state or variable in the intracellular model.'));
        return
    elseif(sum(strcmp(p.States2FitSTV(CountFit), d.ExpStateNames))==0)
        
        p.States2FitSTV(CountFit)
        d.ExpStateNames        
        
        disp(strcat('Warning(Optimization): "', p.States2FitSTV(CountFit), '" has not been measured.'));
        return
    else
        p.PosState2FitSTV(CountFit, 1) = find(strcmp(p.States2FitSTV(CountFit), [p.StateNames; p.VariableNames]));
        p.PosState2FitSTV(CountFit, 2) = find(strcmp(p.States2FitSTV(CountFit), d.ExpStateNames));
    end
end

p.PosState2FitCOI = zeros(0,2);
for CountFit=1:length(p.States2FitCOI)
    if(sum(strcmp(p.States2FitCOI(CountFit), [p.StateNames; p.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.States2FitCOI(CountFit), '" is no state or variable in the intracellular model.'));
        return
    elseif(sum(strcmp(p.States2FitCOI(CountFit), d.ExpStateNames))==0)
        
        p.States2FitCOI(CountFit)
        d.ExpStateNames        
        
        disp(strcat('Warning(Optimization): "', p.States2FitCOI(CountFit), '" has not been measured.'));
        return
    else
        p.PosState2FitCOI(CountFit, 1) = find(strcmp(p.States2FitCOI(CountFit), [p.StateNames; p.VariableNames]));
        p.PosState2FitCOI(CountFit, 2) = find(strcmp(p.States2FitCOI(CountFit), d.ExpStateNames));
    end
end

% Find index for normalization to inital state value
if(strcmpi(p.ObjFunCal, 'normtoinitialstate') || strcmpi(p.ObjFunCal, 'normtostate'))
    if(sum(strcmp(p.NormState, d.ExpStateNames))==0 ||...
       sum(strcmp(p.NormState, [p.StateNames; p.VariableNames]))==0)
        disp(strcat('Warning(Optimization): "', p.NormState, '" is no state in the intracellular model and/or data set.'));
        return
    else
        p.IndexNormSim = find(strcmp(p.NormState, [p.StateNames; p.VariableNames]));
        p.IndexNormExp = find(strcmp(p.NormState, d.ExpStateNames));
    end
end

%% Prepare vector to realize offset in RNA concentrations for intracellular model
RNANames = {'RvSeg5', 'RvSeg7', 'RvSeg9', 'RvSeg8', ...
            'Rm5', 'Rm7', 'Rm9', 'Rm8', ...
            'RcSeg5', 'RcSeg7', 'RcSeg9', 'RcSeg8', ...
            'M1_OP7tot', 'M1_WTtot', 'NPtot', 'PB1tot'}; % need to be part of States2Fit 
        
% RNANames = p.States2FitSTV;
p.OffsetIndexSTV = zeros(0, 2);
if ( ~isempty(d.ExpStateNames) && ~isempty(p.States2FitSTV) )
    for CountRnas = 1:length(RNANames)
        if ( sum(strcmp(RNANames(CountRnas),p.States2FitSTV)) )
            p.OffsetIndexSTV(end+1, 1) = find(strcmp(RNANames(CountRnas), p.States2FitSTV));
            p.OffsetIndexSTV(end, 2) = find(strcmp(RNANames(CountRnas), d.ExpStateNames));
        end
    end
end 

% RNANames = p.States2FitCOI;
p.OffsetIndexCOI = zeros(0, 2);
if ( ~isempty(d.ExpStateNames) && ~isempty(p.States2FitCOI) )
    for CountRnas = 1:length(RNANames)
        if ( sum(strcmp(RNANames(CountRnas),p.States2FitCOI)) )
            p.OffsetIndexCOI(end+1, 1) = find(strcmp(RNANames(CountRnas), p.States2FitCOI));
            p.OffsetIndexCOI(end, 2) = find(strcmp(RNANames(CountRnas), d.ExpStateNames));
        end
    end
end 

%% Consider parameters applied on the "extracellular level" related to post-release calculations
for i = 1 : length(p.PaOpt_Titers)
    p.x0(end+1)  = p.(p.PaOpt_Titers{i});
    p.x_L(end+1) = p.(p.PaOpt_Titers{i})/p.PaDev;
    p.x_U(end+1) = p.(p.PaOpt_Titers{i})*p.PaDev;
end

%% Start parameter estimation

% Create safety measures in case network connection is lost or strange error occurs
p.safety.Bestf = Inf;
p.safety.xbest = ones(size([p.x0],1),1);

AllEstParas = [p.PaOpt p.IcOpt p.PaOpt_Titers];
% Set upper and lower bounds for special parameters
        p.x_L(strcmp('KvRel',AllEstParas)) = 1;
        
        p.x_U(strcmp('FPar_STV',AllEstParas)) = 1;        
        p.x_U(strcmp('FPar_COI',AllEstParas)) = 1;        
        
        if ( ~p.LogPara )
            p.x_L(strcmp('Fadv_cRNA',AllEstParas)) = 1e-3;
            p.x_L(strcmp('Fadv_mRNA',AllEstParas)) = -1;
            p.x_L(strcmp('Fadv_vRNA',AllEstParas)) = 1e-3;        
            p.x_L(strcmp('Fadv_cvRNA',AllEstParas))= -1;       
            p.x_L(strcmp('Fexp',AllEstParas))      = 0;
            p.x_L(strcmp('Fbind',AllEstParas))     = 0;
            
            p.x_U(strcmp('Fadv_cRNA',AllEstParas)) = 10;
            p.x_U(strcmp('Fadv_mRNA',AllEstParas)) = 10;
            p.x_U(strcmp('Fadv_vRNA',AllEstParas)) = 10;    
            p.x_U(strcmp('Fadv_cvRNA',AllEstParas))= 10;
            p.x_U(strcmp('Fexp',AllEstParas))      = 1;
            p.x_U(strcmp('Fbind',AllEstParas))     = 1;                   
        end
        
        if ( p.LogPara )
            p.x_L = log10(p.x_L);
            p.x_U = log10(p.x_U);
            p.x0  = log10(p.x0);
        end                 
        
switch lower(p.Solver)
    case {'fminsearch'}
        %% Fminsearch
        [xbest, result.ObjFunVal] = ...
            fminsearch(@Opt_ObjFun, [p.x0;], p.FminOptions);
        result.xbest = xbest;

    case {'fmincon'}
        [xbest, result.ObjFunVal] = ...
            fmincon(@Opt_ObjFun, [p.x0],[],[],[],[],...
            p.x_L,p.x_U,[],p.FminOptions);
        result.xbest = xbest;

    case {'fssm'}
        %% fSSm
        RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));    

        p.problem.x_0     = [p.x0];
        p.problem.x_L     = [p.x_L];
        p.problem.x_U     = [p.x_U];
        p.problem.int_var = 0;
        p.problem.bin_var = 0;
        p.problem.c_L     = [];
        p.problem.c_U     = [];
        p.problem.N       = numel(d.time.RNA);

        % initialize parameter set recording matrix
        p.record = zeros(0,size(p.problem.x_0,1));

        p.problem.f       = 'Opt_ObjFun';
        [result.fSSm] = fssm_kernel(p.problem, p.fSSmOptions);


        xbest  = result.fSSm.xbest;

        result.xbest  = xbest;
        result.ObjFunVal = result.fSSm.fbest;

    otherwise
        warning(char(['Optimization - No optimization algorithm "',p.solver,'" found!']));
        result = [];
end

