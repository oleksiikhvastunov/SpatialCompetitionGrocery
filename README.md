# SpatialCompetitionGrocery
Estimation of the Demand Model for Spatial Competiton of Grocery Chains

(1) data=CsvToMat(FileNameStructuresStata,FileNameStructuresMatlab) - function that takes Stata produced csv file FileNameStructuresStata and transforms it into matrix data which is saved into mat file FileNameStructuresMatlab   

(2) [ts,storeRevenue]=setupMatlabStructsClub(data,DemandStructsName,year) - function that take data matrix produced by CsvToMat and transform it into structure ts and vector of stores' revenues storeRevenue which are saved to mat file DemandStructsName, year is number which is used for normalization of variables (different years have dfifferent mean income which affects normalization)

(3) mainDemandEstimationMSAClubDistForPaper - main function that contains optimization routine. In the function the parameters that govern it work can be specified, for example if you want to load data from Stata or you already have data structures in mat file. Derivative check can also be specifed. In addition, initial point for optimization routine can be loaded. 

(4) Pass=ObjectiveDerivativeCheck(params,ts,storeRevenue) - function which checks if the value of analytical derivative is close to numerical at the point params. 

(5) [obj, grad]=demandObjective(x,ts,storeRevenue) - function that calculates value of objective function obj and value of analytical derivative at point x. Objective function is a sum of squared differences between logs of store revenue in the data and produced by the model.

(6) u = getUtilityGen(ts, betas) - function that calculates utility of the store from tract t, betas are the subvector of barameters that correspond to utility variables.

(7) ts_shares = getShareGenOPdensity3(ts,u,aalpha)- function that calculates shares of tract grocery budget that go to the stores in the choice set, aalpha parameters that correspond to density and used to calculate utility of the outside option in the tract. 

(8) rev_hat = getRevOPalpha(ts, ts_shares, alpha) - function that calculates revenue of the stores, based on the shares ts_shares and alpha which is a parameter that governs grocery bugdet as a share of income.

(9) [SharesTimesUtilVarByTract] = getSharesTimesUtilVarByTractOP(ts_shares, ts) - function that calculates sum of shares time utility variables across each tract. Output is used to calculate derivative of the objective function and standard errors of the estimates.

(10) [derUtilVar] = gradObjOPalphadensity(SharesTimesUtilVarByTract, ts_shares, ts,   storeRevenue,rev_hat, alpha) - function that calculates analytical derivatives of the objective function. SharesTimesUtilVarByTract are based on the point at which analytical derivative calculated.

(11) [ varcovar, standardErrors ] = getSE_opOPalphadensity( xx, ts, storeRevenue ) - function that calculates standard errors of the estimates xx.

(12) [derErrorByStore] = getErrorDerivativeByStoreOPalphadensity(SharesTimesUtilVarByTract, ts_shares, ts, storeRevenue,rev_hat,alpha) - subroutine of getSE_opOPalphadensity which computes how revenue of every store are affected by change in each parameter. Later on the outer product is used to calculate standard errors of the estimates.
