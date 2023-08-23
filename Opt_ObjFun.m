%%% -------------------------------------------------------------- %%%
%%% Evaluation of the deviation between data and model simulation  %%%
%%% for a model of OP7 and STV coinfection                         %%%
%%% -------------------------------------------------------------- %%%
%
%   authors: Daniel Ruediger
%   last revised: 2023/06/13

function ObjFunVal = Opt_ObjFun(x,varargin)

global p
global d
global stopOptimization % fssm_kernel stop variable

% clear mex

TSTART = tic;

% convert x into a row vector if CMEAS optimizer is used
if ( size(x,1) > 1 && size(x,2) > 1 )
    x = x(:,1)';
elseif ( size(x,1) > 1 )
    x = x';
else
%     error('ObjFun input has wrong format!')
end

% transform parameters if necessary
if ( p.LogPara )
    p.x = 10.^(x(1:p.nOpt));
else
    p.x = x(1:p.nOpt);
end

%% Show current paras
% tmp = [p.PaOpt, p.PaOpt_Titers];
% for i = 1 : length(tmp)
%     fprintf('%s\t%.3e\n',tmp{i},x(i));
% end
% fprintf('\n');

%% Intracellular model
p.Ic_STV = p.Ic;
p.Ic_COI = p.Ic;

p.Ic_STV(strcmp('Vex',p.StateNames)) = p.InfectingParticles_STV;
p.Ic_COI(strcmp('Vex',p.StateNames)) = p.InfectingParticles_COI;

%adopt new initial conditions
if ( ~isempty(p.PosIcOpt) )
    p.Ic_STV(p.PosIcOpt) = p.x(p.nPaOpt+1:p.nOpt);
    p.Ic_COI(p.PosIcOpt) = p.x(p.nPaOpt+1:p.nOpt);
end

%% Simulate system with current parameters
% update post-intracellular parameters
for i = 1 : length(p.PaOpt_Titers)
    p.(p.PaOpt_Titers{i}) = x(p.nOpt+i);
end
p.FPar_COI = p.FPar_STV;

% shorten sim time if only one side is optimized
if ( p.FittingStrategy == 1 )
    SimTimeSTV = p.SimTime;
    SimTimeCOI = [0 0.1];
elseif ( p.FittingStrategy == 2 )
    SimTimeSTV = [0 0.1];
    SimTimeCOI = p.SimTime;
else
    SimTimeSTV = p.SimTime;
    SimTimeCOI = p.SimTime;
end

try           
    % STV infection
    if( p.CompileFlag )
            resultSTV = IQMPsimulate('MexModelFile_Opt', SimTimeSTV, p.Ic_STV, ...
                          [p.PaOpt, {'pInLowSeg', 'pInMedSeg', 'pInHighSeg'}], ...
                          [p.x(1:p.nPaOpt), p.ParticleRatioIn_STV], ...        
                          p.OdeOptions);                      
    else
        tmp = p.Model;
        tmp = IQMparameters(tmp, [p.PaOpt, {'pInLowSeg', 'pInMedSeg', 'pInHighSeg'}], ...
                                 [p.x(1:p.nPaOpt), p.ParticleRatioIn_STV]);

        resultSTV = IQMsimulate(tmp, 'ode23s', SimTimeSTV, p.Ic_STV, p.OdeOptions);
        resultSTV.time = resultSTV.time';        
    end    

    % co-infection
    if ( p.CompileFlag )
        resultCOI = IQMPsimulate('MexModelFile_Opt', SimTimeCOI, p.Ic_COI, ...
                      [p.PaOpt, {'pInLowSeg', 'pInMedSeg', 'pInHighSeg'}], ...
                      [p.x(1:p.nPaOpt), p.ParticleRatioIn_COI], ...   
                      p.OdeOptions);                  
    else
        tmp = p.Model;
        tmp = IQMparameters(tmp, [p.PaOpt, {'pInLowSeg', 'pInMedSeg', 'pInHighSeg'}], ...
                                 [p.x(1:p.nPaOpt), p.ParticleRatioIn_COI]);            

        resultCOI = IQMsimulate(tmp, 'ode23s', SimTimeCOI, p.Ic_COI, p.OdeOptions);
        resultCOI.time = resultCOI.time';                
    end
catch  

    ObjFunVal = 1e9 + rand*1e8;
    return

end

% collect all results
AllTime.STV       = resultSTV.time;
AllProperties.STV = [resultSTV.states, resultSTV.variables];
AllValues.STV     = [resultSTV.statevalues, resultSTV.variablevalues];

AllTime.COI       = resultCOI.time;
AllProperties.COI = [resultCOI.states, resultCOI.variables];
AllValues.COI     = [resultCOI.statevalues, resultCOI.variablevalues];

AllValues.STV(AllValues.STV<0) = 0;
AllValues.COI(AllValues.COI<0) = 0;

if ( resultSTV.statevalues(end,strcmp('SafeGuardFlag',p.StateNames)) == 1 ...
     || resultCOI.statevalues(end,strcmp('SafeGuardFlag',p.StateNames)) == 1 )
    
    ObjFunVal = 2e9 + rand*1e8;
    return 
end

%% Re-organize simulation results

% find time point corresponding to experiments in simulation
tmp = repmat(AllTime.STV,length(d.time.Total),1) - repmat(d.time.Total,1,length(AllTime.STV));
[~,TimePosSTV] = find(repmat(min(abs(tmp),[],2),1,length(AllTime.STV)) == abs(tmp));

tmp = repmat(AllTime.COI,length(d.time.Total),1) - repmat(d.time.Total,1,length(AllTime.COI));
[~,TimePosCOI] = find(repmat(min(abs(tmp),[],2),1,length(AllTime.COI)) == abs(tmp));

SimValues.STV = AllValues.STV(TimePosSTV,p.PosState2FitSTV(:,1));
SimValues.COI = AllValues.COI(TimePosCOI,p.PosState2FitCOI(:,1));

%% OFFSET - intracellular level
% To account for free viral RNAs in the seed virus which might adhere to
% cells and cause a constant offset in measurements, we add the 0 hpi
% simulation value of each RNA species to the subsequent time points   
if ( ~isempty(d.ExpStateNames) && ( ~isempty(p.States2FitSTV) || ~isempty(p.States2FitCOI) ) )

        AddOffset = 10.^d.RKI.ESV(1,p.OffsetIndexSTV(:,2));  
        SimValues.STV(:,p.OffsetIndexSTV(:,1)) = SimValues.STV(:,p.OffsetIndexSTV(:,1)) ...
                                              + repmat(AddOffset, length(d.time.Total), 1);                                          

        AddOffset = 10.^d.OP7.ESV(1,p.OffsetIndexCOI(:,2));
        SimValues.COI(:,p.OffsetIndexCOI(:,1)) = SimValues.COI(:,p.OffsetIndexCOI(:,1)) ...
                                              + repmat(AddOffset, length(d.time.Total), 1);                                                   
end

% calculate titers based on total cell numbers
TiterNames = {'S5', 'S7_WT', 'S7_OP7', 'HA', 'TCID'};
for i = 1 : length(p.States2FitSTV)
    if ( sum(strcmp(p.States2FitSTV(i),TiterNames)) )
       if ( strcmp(p.States2FitSTV(i),'TCID') )           
            SimValues.STV(:,strcmp(p.States2FitSTV(i),p.States2FitSTV)) = ...
              d.NumInfectedCells * SimValues.STV(:,strcmp(p.States2FitSTV(i),p.States2FitSTV)) ...
                .* exp(-p.kDegV*AllTime.STV(TimePosSTV))' * p.FPar_STV;    
       else
            SimValues.STV(:,strcmp(p.States2FitSTV(i),p.States2FitSTV)) = ...
              d.NumInfectedCells * SimValues.STV(:,strcmp(p.States2FitSTV(i),p.States2FitSTV));           
       end
    end
end
for i = 1 : length(p.States2FitCOI)
    if ( sum(strcmp(p.States2FitCOI(i),TiterNames)) )         
       if ( strcmp(p.States2FitCOI(i),'TCID') )
            SimValues.COI(:,strcmp(p.States2FitCOI(i),p.States2FitCOI)) = ...
              d.NumInfectedCells * SimValues.COI(:,strcmp(p.States2FitCOI(i),p.States2FitCOI)) ...
                .* exp(-p.kDegV*AllTime.COI(TimePosCOI))' * p.FPar_COI;
       else
            SimValues.COI(:,strcmp(p.States2FitCOI(i),p.States2FitCOI)) = ...
              d.NumInfectedCells * SimValues.COI(:,strcmp(p.States2FitCOI(i),p.States2FitCOI));           
       end
    end
end    
    
%% RKI data - calculate objective function values
MeasuredValues.STV    = 10.^d.RKI.ESV(:,p.PosState2FitSTV(:,2));
MeasuredValuesStd.STV = 10.^d.RKI.ESVStd(:,p.PosState2FitSTV(:,2));

% Define weighting matrix for extracellular model + logarithmization
switch lower(p.ObjFunCal)
    case {'normal'}

    case {'log'}                
        SimValues.STV         = log10(SimValues.STV);
        MeasuredValues.STV    = log10(MeasuredValues.STV);
        MeasuredValuesStd.STV = log10(MeasuredValuesStd.STV);
        
        SimValues.STV(SimValues.STV==-Inf)           = 0;
        MeasuredValues.STV(MeasuredValues.STV==-Inf) = 0;                            
        MeasuredValuesStd.STV(MeasuredValuesStd.STV==-Inf) = 0;
        MeasuredValuesStd.STV(MeasuredValuesStd.STV==0)    = 1.5;
    otherwise
        error('ERRTXT');                
end

DevStatesSTV = (SimValues.STV - MeasuredValues.STV) ./ MeasuredValuesStd.STV;
DevStatesSTV(isnan(DevStatesSTV)) = 0;

ObjFunValSTV = sum(sum((DevStatesSTV(2:end,:)).^2).*p.RelStateWeightSTV);

%% OP7 data - calculate objective function values
MeasuredValues.COI    = 10.^d.OP7.ESV(:,p.PosState2FitCOI(:,2));
MeasuredValuesStd.COI = 10.^d.OP7.ESVStd(:,p.PosState2FitCOI(:,2));

% Define weighting matrix for extracellular model + logarithmization
switch lower(p.ObjFunCal)
    case {'normal'}

    case {'log'}                
        SimValues.COI         = log10(SimValues.COI);
        MeasuredValues.COI    = log10(MeasuredValues.COI);
        MeasuredValuesStd.COI = log10(MeasuredValuesStd.COI);
        
        SimValues.COI(SimValues.COI==-Inf)                 = 0;
        MeasuredValues.COI(MeasuredValues.COI==-Inf)       = 0;        
        MeasuredValuesStd.COI(MeasuredValuesStd.COI==-Inf) = 0;
        MeasuredValuesStd.COI(MeasuredValuesStd.COI==0)    = 1.5;        
    otherwise
        error('ERRTXT');                
end

DevStatesCOI = (SimValues.COI - MeasuredValues.COI) ./ MeasuredValuesStd.COI;
DevStatesCOI(isnan(DevStatesCOI)) = 0;

ObjFunValCOI = sum(sum((DevStatesCOI(2:end,:)).^2).*p.RelStateWeightCOI);

%% Combine objective function values
ObjFunValSTV = sum(ObjFunValSTV(~isnan(ObjFunValSTV))); 
ObjFunValCOI = sum(ObjFunValCOI(~isnan(ObjFunValCOI))); 

ObjFunVal = (ObjFunValSTV + ObjFunValCOI) / ...
       ( (numel(MeasuredValues.STV) - sum(sum(isnan(MeasuredValues.STV)))) ...
         + (numel(MeasuredValues.COI) - sum(sum(isnan(MeasuredValues.COI)))) );

% save best result in case of crash (does not work during parallel computation)
if ( ObjFunVal < p.safety.Bestf )
    p.safety.Bestf = ObjFunVal;
    p.safety.xbest = x;
    
%     fprintf('\nNew Bestf: %.2f\n',ObjFunVal*1e6);
end

% fprintf('\nNew ObjFun: %.6f, reached in %.1f s\n',ObjFunVal,toc(TSTART)); 

%% optimization time calculation + display
p.TotalOptTime = p.TotalOptTime + toc(TSTART);
if ( p.TotalOptTime > p.MaxTime )
   stopOptimization = 1; 
elseif ( p.TotalOptTime > p.TotalTimeCounter*3600 )
    fprintf('\n####################### Current simulation time: %i h #######################\n\n', p.TotalTimeCounter);
    p.TotalTimeCounter = p.TotalTimeCounter + 1;
end

return

%% use safety
Opt.xbest = p.safety.xbest;
Opt.ObjFunVal = p.safety.Bestf;

Opt.In.xbest = bestever.x(1:p.nOpt)';
Opt.Ex.xbest = bestever.x(p.nOpt+1:end)';
Opt.ObjFunVal = bestever.f;

