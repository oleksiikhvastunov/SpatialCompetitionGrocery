clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the type of estimation 

% year can be 2006 
year=2006;

% typeModel 0 - fixed effect is the same for all chains
% typeModel 1 - fixed effect has different intercept for chains and same slope
% typeModel 2 - fixed effect has different intercept and slope for chains
typeModel=2;

% loadStataCsv 0 - do not load files from stata and use already created mat
% file

% loadStataCsv 1 - load files from stata and create mat
% file
loadStataCsv=0;

% derivativeCheck 0 - do not check analytical=numerical derivative 
% derivativeCheck 1 - check analytical=numerical derivative
derivativeCheck=0;

% loadStartPoint 0 - do not load starting point, initial point is all zeros
% but budget parameter 
% loadStartPoint 1 - load starting point which is based on the point at
% which convergence was achived (constant is added to all parameters)
% loadStartPoint 2 - load starting point which is based on typeModel==0
% where all chain effects are the same
loadStartPoint=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read csv file and create mat file or read data directly from mat file

if (year==2006)
    % Name of the csv file produces by Stata
    FileNameStructuresStata='StataStructures2010MSAClubForPaper.csv';
    
    % Name of the Matlab file which contains matrix from csv file
    FileNameStructuresMatlab='MatlabStructures2010MSAClubForPaper.mat';
    
    % Name of the Matlab file which contains data structures for the main
    % demand estimation script
    DemandStructsName='demandStructsOP2010MSAClubForPaper.mat';
    if (loadStataCsv==1)
        data=CsvToMat(FileNameStructuresStata,FileNameStructuresMatlab);
        [ts,storeRevenue]=setupMatlabStructsClub(data,DemandStructsName,year);
    end;
    if (loadStataCsv==0)
        load(DemandStructsName);
    end;
end;

%% Load initial point if necessary 

numChainsC=max(ts.chainIDC);

storeR=storeRevenue;


storeRevenue=storeR;

switch loadStartPoint
    case 0
        % Point which is not based on any previous estimations. All
        % parameters but the one which governs badget are zero.
        params=zeros(17+2*(typeModel==0)+(1+numChainsC)*(typeModel==1)+(2*numChainsC)*(typeModel==2),1);
        params(end,1)=0.12;
    case 1 
        % Constant is added to all parameters which represent a point at
        % which convergence were achived
        % Converges within 30-40 minutes
        load('resultsMSAClubDistForPaperHanaf_2_2006.mat','view');
        params=view(:,1);
        params=params+0.01;
    case 2 
        % Results where all chain intercepts and slopes are the same (modelType==0) are
        % taken as initial guess, converges within 60-80 minutes
        load('resultsMSAClubDist_0_2006.mat','view');
        params=[view(1:14,1);ones(numChainsC,1)*view(15,1);ones(numChainsC,1)*view(16,1);view(end-2:end,1)];
end;


    
% ts - is a structure, which includes utility variables. Utility variables
% are divided into common (the ones common for models 0 to 2) and different 
% parts (different part is saved into cells that correspond to a particular model) 

% Create a full utility Variables structure and remove auxiliary, which
% were used to save space and have structures for all three types of models
% for a particular year in one file

if (typeModel==0)
    ts.utilVar=[ts.utilVarCommon,ts.utilVarDifferent{1}];
end;
if (typeModel==1)
    ts.utilVar=[ts.utilVarCommon,ts.utilVarDifferent{2}];
end;
if (typeModel==2)
    ts.utilVar=[ts.utilVarCommon,ts.utilVarDifferent{3}];
end;
ts.utilVarCommon=[];
ts.utilVarDifferent=[];


% Compare analytical and numerical derivative
if (derivativeCheck==1)
    Pass=ObjectiveDerivativeCheck(params,ts,storeRevenue);
end;

if ((derivativeCheck==0) || (Pass==1))
    options=optimset('Display','iter','TolFun',1e-6,'TolX',1e-6,'GradObj','on');

    exitflag=-999;
    maxtrials=1;
    trial=0;

    
    while ((trial<maxtrials) && (exitflag<0))
        %Create upper and lower bounds for alpha
        ub = inf*ones(size(params));
        lb = -inf*ones(size(params));
        %Set alpha, the final parameter:
        ub(end) = 1;
        lb(end) = 0;
        [xx, fval, exitflag, output] =knitromatlab(@(x) demandObjective(x,ts,storeRevenue),params,[],[],[],[],lb,ub,[],[],options,'koptions.opt');
        % If the algorithm gets stuck it starts in the point nearby
        %params=xx+randn(size(params,1),1)/100;
        params=xx;
        trial=trial+1;
        exitflag
    end;
    save(strcat('resultsMSAClubDistxxForPaperHanafCheck_',num2str(typeModel),'_',num2str(year)),'xx','-v7.3');
    
    % Computes variance-covariance matrix and standard errors
    [ varcovar, standardErrors ] = getSE_opOPalphadensity( xx, ts, storeRevenue );
    % Clear view which was loaded from initial point, if the initial point
    % has not been loaded the command just does not do anything
    clear view;
    % Estimated Parameters, Standard Errors and t-values are created
    % This type of matrix is used as an imput to functions that produce latex tables 
    view=[xx,standardErrors,xx./standardErrors];
    % Results are saved in a mat file
    save(strcat('resultsMSAClubDistForPaperHanafCheck_',num2str(typeModel),'_',num2str(year)),'view','storeRevenue','fval','typeModel','-v7.3');

end;
