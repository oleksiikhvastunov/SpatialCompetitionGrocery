function [ share ] = getShareGenOPdensity3(ts, u, aalpha)
    
    % Function that calculats stores' shares based on their utilities u
    % aalpha is used to calculate utilities of outside option

    
    exp_u = exp(u);
    
    % Log Desity bounded below by 1000 which is inline with Holmes.
    capdensity=log(max(ts.density/1000,1));
    
    u0=aalpha(1)*capdensity+aalpha(2)*(capdensity.^2);
    %denom will be used for the denominator (sum of all exp_u within a tract);
    %ordered by tract
    denom = accumarray(ts.tractID, exp_u);
    share = exp_u ./ (exp(u0)+denom(ts.tractID));
end

