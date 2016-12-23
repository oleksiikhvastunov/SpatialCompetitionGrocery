function [ varcovar, standardErrors ] = getSE_opOPalphadensity( params, ts, storeRevenue )

    % function that calculates standard errors of the estimated parameters
    % params
      
    betas=params(1:end-3);
    aalpha=params(end-2:end-1);
    alpha=params(end);
    
    u = getUtilityGen(ts, betas); 
    ts_shares = getShareGenOPdensity3(ts,u,aalpha);
    rev_hat = getRevOPalpha(ts, ts_shares, alpha);  
    
    numStores=length(storeRevenue);
    numParams=length(params);
    
     
    
    [SharesTimesUtilVarByTract] = getSharesTimesUtilVarByTractOP(ts_shares, ts);
    % derErrorByStore(i,j) - derivative of store i revenue with respect to
    % parameter j
    [derErrorByStore] = getErrorDerivativeByStoreOPalphadensity(SharesTimesUtilVarByTract, ts_shares, ts, storeRevenue,rev_hat,alpha);
    
    
    varcovar=zeros(numParams,numParams);
    outerProduct=zeros(numParams,numParams);
    for i=1:1:numStores
        
        outerProduct=outerProduct+(1/numStores)*derErrorByStore(:,i)*derErrorByStore(:,i)';
    end;
    varcovar=inv(outerProduct)/numStores;
    standardErrors=diag(varcovar).^0.5;
end

