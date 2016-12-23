function [derErrorByStore] = getErrorDerivativeByStoreOPalphadensity(SharesTimesUtilVarByTract, shares, ts, storeRevenue,rev_hat,alpha)
    % This function produces matrix derErrorByStore where
    % derErrorByStore(i,j) - derivative of store i revenue with respect to variable j.
    % ContractedForUtilVar(i,j) - Derivative of store i revenue with respect constant from chain j fixed effect 
%   twoTimesRevenueDifference ... Derivative of objective function with respect to rev_hat. 


    numStores = max(ts.storeID);
    numUtilVar=size(ts.utilVar,2);

    twoTimesRevenueDifference=2*((log(rev_hat)- log(storeRevenue))./rev_hat);

    
    ContractedForUtilVar=zeros(numStores,numUtilVar);
    forUtilVar=zeros(numStores,numUtilVar);
    [fordensityvar]=getSharesTimesdensity(ts, shares, alpha);
    
    for i=1:1:numUtilVar
        ContractedForUtilVar(:,i)=accumarray(ts.storeID, (ts.inc).*ts.pop.*(shares.*ts.utilVar(:,i)-shares.*SharesTimesUtilVarByTract(:,i)));
        
        forUtilVar(:,i)=ContractedForUtilVar(:,i)'.*twoTimesRevenueDifference';
           
    end;
    
    derErrorByStore=[forUtilVar';repmat(twoTimesRevenueDifference',2,1).*fordensityvar';twoTimesRevenueDifference'.*rev_hat'/alpha];
    
end

