function [ obj, grad ] = demandObjective( params, ts, storeRevenue )
    % Objective function
    % obj - value of the objective function,
    % grad - analytical gradient of the Objective function
    % params - value of the parameters
    % ts - data structure
    % storeRevenue - vector of store Revenues

     
  
    % Unpack the parameter vector into betas which correspond to utility
    % Variables, aalpha which correspond to density variables in the
    % outside option and alpha which governs food budget as a share of
    % income
    betas=params(1:end-3);
    aalpha=params(end-2:end-1);
    alpha=params(end);
    
    % Calculate utilities of the store ts.storeID which are in the choice set of ts.tractID 
    % (u is indexed by store-tract pair)
    u = getUtilityGen(ts, betas); 
    
    % Calculate shares based on utilities u, ts_shares have the same size
    % as u
    ts_shares = getShareGenOPdensity3(ts,u,aalpha);
    
    % Calculate revenue of the stores predicted by the model
    rev_hat = getRevOPalpha(ts, ts_shares, alpha); 
    obj = sum((log(rev_hat) - log(storeRevenue)).^2); 
    
    % If function has more than 1 outputs analytical derivative is calculated 
    if nargout > 1        
        [SharesTimesUtilVarByTract] = getSharesTimesUtilVarByTractOP(ts_shares, ts);
        [derUtilVar] = gradObjOPalphadensity(SharesTimesUtilVarByTract, ts_shares, ts,   storeRevenue,rev_hat, alpha);
        grad=derUtilVar;
    end 
end

