function [derUtilVar ] = gradObjOPalphadensity(SharesTimesUtilVarByTract, shares, ts, storeRevenue,rev_hat, alpha)
%   This function constructs the derivative of the objective function. 
%
%   ContractedForUtilVar(i,j) - Derivative of store i revenue with respect constant from chain j fixed effect 
%   twoTimesRevenueDifference - Derivative of objective function with respect to rev_hat. 


    numStores=max(ts.storeID);
    numChains=max(ts.chainID);
    numChainsC=max(ts.chainIDC);

    twoTimesRevenueDifference=2*((log(rev_hat)- log(storeRevenue))./rev_hat);

    numUtilVar=size(ts.utilVar,2);
    ContractedForUtilVar=zeros(numStores,numUtilVar);
    derUtilVarExtended=zeros(numUtilVar,1);
    
    
    for i=1:1:numUtilVar
        ContractedForUtilVar(:,i)=accumarray(ts.storeID,ts.inc.*ts.pop.*(shares.*ts.utilVar(:,i)-shares.*SharesTimesUtilVarByTract(:,i)));
        
        derUtilVarExtended(i,1)=ContractedForUtilVar(:,i)'*twoTimesRevenueDifference;
        
    end;
    
    sharestimesdensity=getSharesTimesdensity(ts, shares, alpha);
    derUtilVar=[alpha*derUtilVarExtended;sharestimesdensity'*twoTimesRevenueDifference;twoTimesRevenueDifference'*rev_hat/alpha];
    
    
end

