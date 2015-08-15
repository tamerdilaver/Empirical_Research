/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Simulate HF data and compute RV_sec, RV_5min and RK for simulated data
**
**  Date:
**    	24/01/2015
**
**  Author:
**	  	Tamer Dilaver, Koen de Man & Sina Zolnoor
**
**	Supervisor:
**		L.H. Hoogerheide & S.J. Koopman
**
*/


#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>

/*
**  Function:	Compute Newey-West standard errors (code from Hoogerheide)
**
**  Input: 		mX: Matrix with prices
**				vResiduals: Vector with residuals
**
**  Output: 	vNW_std_errors: Newey-West standard errors
*/
fNeweyWestStandardError(const mX, const vResiduals) {

	decl iTmp, iT, iK, iM, vIotaT, vOneToT, mAbs_i_minus_j,
	mWeights, mOmegaHat, mInvXX, mCovOLS_estimator, vNW_std_errors;

	iT 	= rows(mX);
	iK 	= columns(mX);

	iM = ceil(4*(iT/100)^(2/9));
	vIotaT = ones(iT,1);
		
	vOneToT = zeros(iT,1);
	for (iTmp = 0;iTmp<iT; ++iTmp){
		vOneToT[iTmp][0] = iTmp+1;
	}

	mAbs_i_minus_j = fabs(vIotaT*vOneToT'-vOneToT*vIotaT');
	mWeights = (ones(iT,iT) - mAbs_i_minus_j*(1/iM)) .* (mAbs_i_minus_j .<= iM);
	mOmegaHat = (vResiduals * vResiduals') .* mWeights;

	mInvXX = invertsym(mX'*mX);
	mCovOLS_estimator = (iT/(iT-iK)) * mInvXX * (mX'*mOmegaHat*mX) * mInvXX;
	vNW_std_errors = sqrt(diagonal(mCovOLS_estimator)');
	
	return vNW_std_errors;
}

fDieboldMariano(const vForecast1, const vForecast2)
{
	decl vLossDifferential, dMean, dNeweyWestStandardError, dDieboldMarianoStatistic; 
	
	/*--- Define the loss differential (KLIC) ---*/
	vLossDifferential = vForecast1 - vForecast2;	  //log removed upon Hoogerheides recommendation
	
	/*--- Recalculate the HAC variance of the loss differential ---*/
	dMean = meanc(vLossDifferential);
	dNeweyWestStandardError = fNeweyWestStandardError(ones(sizerc(vLossDifferential),1), vLossDifferential - dMean);	  

	/*--- Retrieve the Diebold Mariano statistic ---*/
	dDieboldMarianoStatistic = dMean/dNeweyWestStandardError;
	//println("DM-statistic: ", dDieboldMarianoStatistic, "\n\n" );
	return dDieboldMarianoStatistic;
}

//fDieboldMariano(const vForecast1, const vForecast2);
//fNeweyWestStandardError(const mX, const vResiduals);

main()
{
	decl iModels, mAE, mSE;

	iModels = 4; 										//nGAS, tGAS, alGAS and astGAS
	mAE = 		loadmat("vAE_RV_GARCH.xls");	//vector one	
	mAE = mAE ~ loadmat("vAE_RV_tGAS.xls");		//vector two
	mAE = mAE ~ loadmat("vAE_RV_alGAS.xls"); 	//etc...
	mAE = mAE ~ loadmat("vAE_RV_astGAS.xls");

	mAE = mAE ~ loadmat("vAE_BV_GARCH.xls"); 	
	mAE = mAE ~ loadmat("vAE_BV_tGAS.xls");	
	mAE = mAE ~ loadmat("vAE_BV_alGAS.xls"); 	
	mAE = mAE ~ loadmat("vAE_BV_astGAS.xls");

	mAE = mAE ~ loadmat("vAE_RK_GARCH.xls"); 	
	mAE = mAE ~ loadmat("vAE_RK_tGAS.xls");	
	mAE = mAE ~ loadmat("vAE_RK_alGAS.xls"); 	
	mAE = mAE ~ loadmat("vAE_RK_astGAS.xls");

	mSE = 		loadmat("vSE_RV_GARCH.xls");	//vector one	
	mSE = mSE ~ loadmat("vSE_RV_tGAS.xls");		//vector two
	mSE = mSE ~ loadmat("vSE_RV_alGAS.xls"); 	//etc...
	mSE = mSE ~ loadmat("vSE_RV_astGAS.xls");

	mSE = mSE ~ loadmat("vSE_BV_GARCH.xls"); 	
	mSE = mSE ~ loadmat("vSE_BV_tGAS.xls");	
	mSE = mSE ~ loadmat("vSE_BV_alGAS.xls"); 	
	mSE = mSE ~ loadmat("vSE_BV_astGAS.xls");

	mSE = mSE ~ loadmat("vSE_RK_GARCH.xls"); 	
	mSE = mSE ~ loadmat("vSE_RK_tGAS.xls");	
	mSE = mSE ~ loadmat("vSE_RK_alGAS.xls"); 	
	mSE = mSE ~ loadmat("vSE_RK_astGAS.xls");

	decl mTableAE, mTableSE;

	mTableAE = mTableSE = .NaN * ones((iModels-1)*3,(iModels-1)*3);
	

	for(decl i = 1; i<iModels;i++){
		for(decl j=0; j<i; j++){
			if(i!=j){
				//RV
				mTableAE[(i-1)*3][j*3] 		= fDieboldMariano(mAE[][i], mAE[][j]); //compute DM
				mTableSE[(i-1)*3][j*3] 		= fDieboldMariano(mSE[][i], mSE[][j]); //compute DM
				//BV
				mTableAE[(i-1)*3+1][j*3+1]	= fDieboldMariano(mAE[][i+iModels], mAE[][j+iModels]); //compute DM
				mTableSE[(i-1)*3+1][j*3+1]	= fDieboldMariano(mSE[][i+iModels], mSE[][j+iModels]); //compute DM
				//RK
				mTableAE[(i-1)*3+2][j*3+2] 	= fDieboldMariano(mAE[][i+2*iModels], mAE[][j+2*iModels]); //compute DM
				mTableSE[(i-1)*3+2][j*3+2] 	= fDieboldMariano(mSE[][i+2*iModels], mSE[][j+2*iModels]); //compute DM
			}	
		}
	
	}

/*
** ..................................................................................			
**	These two tables are the same....???
** ..................................................................................
*/	
	println("DM Absolute Error");
	print(mTableAE,"\n");
	println("DM Squared Error");
	print(mTableSE,"\n");



/*
** ..................................................................................			
**	Here I'm just fooling around to test how they can be the same ....???
** ..................................................................................
*/	
	decl v1, v2, v3, v4;

	v1 = loadmat("vAE_BV_tGAS.xls");
	v2 = loadmat("vAE_BV_astGAS.xls");
	v3 = loadmat("vSE_BV_tGAS.xls");
	v4 = loadmat("vSE_BV_astGAS.xls");

	//print((v1~v2~v3~v4)[0:9][]);
	
	println("DM-statistic AE: ", fDieboldMariano(v1,v2), "\n\n" );
	println("DM-statistic: SE", fDieboldMariano(v3,v4), "\n\n" );	

	
	
}





