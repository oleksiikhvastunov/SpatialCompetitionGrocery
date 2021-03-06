function [fordensityvar]=getSharesTimesdensity(ts, shares, alpha)
    densityvar=[log(max(ts.density/1000,1)),(log(max(ts.density/1000,1))).^2];
    inopt=accumarray(ts.tractID,shares);
    outopt=1-inopt;
    outoptextended=outopt(ts.tractID);
    extended=-1*alpha*[ts.pop.*ts.inc.*shares.*outoptextended.*densityvar(:,1),ts.pop.*ts.inc.*shares.*outoptextended.*densityvar(:,2)];
    fordensityvar=[accumarray(ts.storeID,extended(:,1)),accumarray(ts.storeID,extended(:,2))];