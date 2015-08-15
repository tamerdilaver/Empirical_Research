/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Estimate all GARCH model parameters (gamma, omega, alpha and beta)
**		with Maximum Likelikhood for SBUX return s.t. alpha + beta <1 (Since alpha>0 and beta>0) 
**
**  Date:
**    	10/01/2015
**
**  Author:
**	  	Tamer Dilaver, Koen de Man & Sina Zolnoor
**
**	Supervisor:
**		L.F. Hoogerheide & S.J. Koopman
**
*/

#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

static decl iB;	 					//Repeats
static decl iSIZE;					//Size of time series
static decl iSTEPS;					//#Steps to divide the size
static decl iSIMS;					//# of Zt ~ N(0,1)
static decl dALPHA;
static decl dBETA;
static decl dOMEGA;
static decl dGAMMA;
static decl vSTD_NORM;				// Zt ~ N(0,1)
static decl s_vY; 					//Simulated returns
static decl s_vDate;
static decl dRATIO;
static decl dALPHA_START;
static decl dBETA_START;
static decl dOMEGA_START;
static decl dGAMMA_START;

/*
**  Function:	Transform (start)parameters	  Alpha, Beta, Omega, Gamma Startingvalues
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/

fTransform(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]);
	avThetaStar[0][1] = log(vTheta[1]);
	avThetaStar[0][2] = log(vTheta[2]);
	avThetaStar[0][3] = log(vTheta[3]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adAlpha, adBeta, aOmega, adGamma,, vTheta
**
**  Output: 	1 
*/

fGetPars(const adAlpha, const adBeta, const adOmega, const adGamma, const vTheta){

	adAlpha[0] = exp(vTheta[0]);
	adBeta[0] = exp(vTheta[1]);
	adOmega[0] = exp(vTheta[2]);
	adGamma[0] = exp(vTheta[3]);
	return 1;
}

/*
**  Function:	Calculates average value loglikelihood for GARCH given parameter values
**
**  Input: 		vTheta [parameter values], adFunc [adress function value], avScore [the score], amHessian [hessian matrix]
**
**  Output:		1
**
*/

fLogLike_Garch(const vTheta, const adFunc, const avScore, const amHessian){
	decl dAlpha, dBeta, dOmega, dGamma;
	fGetPars( &dAlpha,  &dBeta, &dOmega, &dGamma, vTheta);

	//decl dS2 = dOmega/(1-dAlpha-dBeta);		
	decl dS2 = dGamma;							//initial condition by definition
	decl vLogEta = zeros(sizerc(s_vY), 1);

	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vLogEta[i] = log(M_2PI) +log(dS2) + s_vY[i]^2 / dS2;		//Gaussian
						
			//GARCH recursion
			dS2 = dOmega + dBeta* dS2 +  dAlpha* (s_vY[i]^2 - dS2);
	}
	
	adFunc[0] = sumc(vLogEta)/(-2*sizerc(s_vY)); 									 	//Average
	return 1;
}

/*
**  Function:	Transform parameters back	Alpha, Beta, Omega, Gamma Startingvalues
**
**  Input: 		vThetaStar
**
**  Output: 	vTheta [parametervalues]
*/

fTransformBack(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0]);
	avTheta[0][1] = exp(vThetaStar[1]);
	avTheta[0][2] = exp(vThetaStar[2]);
	avTheta[0][3] = exp(vThetaStar[3]);
	return 1;
}

/*
**  Function:	calculate standard errors
**
**  Input: 		vThetaStar
**
**  Output: 	vStdErrors
*/

fSigmaStdError(const vThetaStar){

 		decl iN, mHessian, mHess, mJacobian, vStdErrors, vP;

		iN 			= sizerc(s_vY);
		Num2Derivative(fLogLike_Garch, vThetaStar, &mHessian);
		//NumJacobian(fTransformBack, vThetaStar, &mJacobian);	  	//numerical Jacobian
		//mHessian 	= mJacobian*invert(-iN*mHessian)*mJacobian';
		mHessian 	= invertgen(-iN*mHessian);
		vStdErrors 	= exp(vThetaStar).*sqrt(diagonal(mHessian)');	//analytisch

		return 	vStdErrors;
}

/*
**  Function:	calculate variance of model
**
**  Input: 		vTheta
**
**  Output: 	vH [vector with variances]
*/

fVariance(const vTheta){
	decl dAlpha, dBeta, dOmega, dGamma, vH;
	fGetPars(&dAlpha, &dBeta, &dOmega, &dGamma, vTheta);
	
	vH = zeros(sizerc(s_vY),1);
	//vH[0]= dOmega/(1-dAlpha-dBeta);
	vH[0]= dGamma;	
	
	for(decl i = 1; i < sizerc(s_vY); i++){													   //mixed 	
		vH[i] = dOmega + dAlpha*(s_vY[i-1]^2-vH[i-1]) + dBeta*vH[i-1];
	}
	
	return 	vH;
}

/*
**  Function:	Estimate Garch parameters
**
**  Input: 		vReturns, adAlpha_hat, adBeta_hat, adOmega_hat, adGamma_hat
**
**  Output: 	vTheta [estimated parametervalues]
*/

fEstimateGarch(const vReturns, const adAlpha_hat, const adBeta_hat, const adOmega_hat, const adGamma_hat, const avVariance){

	//initialise parameter values
	decl vTheta = zeros(iPARS,1);
	vTheta[0] = dALPHA_START;
	vTheta[1] = dBETA_START;
	vTheta[2] = dOMEGA_START;
	vTheta[3] = dGAMMA_START;
	
	decl vThetaStart = vTheta;

	//globalize returns and vectorize true pars
	s_vY = vReturns;

	//transform parameters
	decl vThetaStar; 
	fTransform(&vThetaStar, vTheta);

	//Maximize the LL
	decl dFunc;
	decl iA;
	iA=MaxBFGS(fLogLike_Garch, &vThetaStar, &dFunc, 0, TRUE);

	//Transform thetasStar back
  	fTransformBack(&vTheta, vThetaStar);

	//return alpha, beta, omega
	adAlpha_hat[0] = vTheta[0];
	adBeta_hat[0] = vTheta[1];
	adOmega_hat[0] = vTheta[2];
	adGamma_hat[0] = vTheta[3];

	decl vSigmaStdError = fSigmaStdError(vThetaStar);
	decl vVariance = fVariance(vThetaStar);
	avVariance[0] = vVariance;
	
	print("\n",MaxConvergenceMsg(iA));
	println("\nFunctiewaarde likelihood eindwaardes:", dFunc);
	print("\nOptimale parameters met standaarderrors \n",
          	"%r", { "A",  "B",  "omega", "gamma"},
          	"%c", {"thetaStart","theta","std.error"}, vThetaStart~vTheta~vSigmaStdError);

	return 1;
}

/*
**  Function:	Determine Forecast
**
**  Input: 		vTheta
**
**  Output: 	vH [vector of forecasts]
*/

fForecast(const vTheta){
	decl dAlpha,  dBeta, dOmega, dGamma, vH;
	fGetPars(&dAlpha,  &dBeta, &dOmega, &dGamma, vTheta);
	
	vH = zeros((sizerc(s_vY)+1),1);
	vH[0]= dGamma;
	
	for(decl i = 0; i < sizerc(s_vY); i++){			
		vH[i+1] = dOmega + dAlpha*(s_vY[i]^2 - vH[i]) + dBeta*vH[i];
	}
	
	return 	vH[sizerc(s_vY)];
}

/*
**  Function:	Compute MAE
**
**  Input: 		adMAE_OC [adress of MAE], vReturns_1 [return series], vBenchmark [Benchmark], dC [ratio]
**
**  Output: 	1
*/

fMAE(const vReturns_1, const vRV, const vBV, const vRK){

	decl iWindow = 250;
	decl iT = sizerc(vReturns_1);
	decl vH_forecast, vAbs_error_RV, vSqrd_error_RV,vAbs_error_BV, vSqrd_error_BV, vAbs_error_RK, vSqrd_error_RK; 

	vH_forecast = vAbs_error_RV = vSqrd_error_RV =vAbs_error_BV = vSqrd_error_BV = vAbs_error_RK = vSqrd_error_RK = zeros(iWindow, 1);

	dALPHA_START = 0.1;
	dBETA_START = 0.99;
	dOMEGA_START = 0.01;
	dGAMMA_START = 0.1;

	for(decl j = 0; j<iWindow; j++){
		s_vY = 	vReturns_1[j:(iT - iWindow +j-1)];

		//initialise parameter values
		decl vTheta = zeros(iPARS,1);
		vTheta[0] = dALPHA_START;
		vTheta[1] = dBETA_START;
		vTheta[2] = dOMEGA_START;
		vTheta[3] = dGAMMA_START;
	
		//transform parameters
		decl vThetaStar; 
		fTransform(&vThetaStar, vTheta);
	
		//Maximize the LL
		decl dFunc;
		decl iA;
		iA=MaxBFGS(fLogLike_Garch, &vThetaStar, &dFunc, 0, TRUE);
	
		//Transform thetasStar back
	  	fTransformBack(&vTheta, vThetaStar);

		dALPHA_START = vTheta[0];
		dBETA_START = vTheta[1];
		dOMEGA_START = vTheta[2];
		dGAMMA_START = vTheta[3];

		vH_forecast[j] = fForecast(vThetaStar);
		//print(vH_forecast[j]~vBenchmark[(iT - iWindow +j)]);
		vAbs_error_RV[j] 	= fabs(vH_forecast[j] - vRV[(iT - iWindow +j)]);
		vSqrd_error_RV[j] =  sqr(vH_forecast[j] - vRV[(iT - iWindow +j)]);

		vAbs_error_BV[j] 	= fabs(vH_forecast[j] - vBV[(iT - iWindow +j)]);
		vSqrd_error_BV[j] =  sqr(vH_forecast[j] - vBV[(iT - iWindow +j)]);

		vAbs_error_RK[j] 	= fabs(vH_forecast[j] - vRK[(iT - iWindow +j)]);
		vSqrd_error_RK[j] =  sqr(vH_forecast[j] - vRK[(iT - iWindow +j)]);

	}

	savemat("vAE_RV_GARCH.xls", vAbs_error_RV);
//	adMAE_OC[0] = meanc(vSqrd_error);
		
	savemat("vSE_RV_GARCH.xls", vSqrd_error_RV);
//	adMSE_OC[0] = meanc(vSqrd_error2);

	savemat("vAE_BV_GARCH.xls", vAbs_error_BV);
//	adMAE_OC[0] = meanc(vSqrd_error);
		
	savemat("vSE_BV_GARCH.xls", vSqrd_error_BV);
//	adMSE_OC[0] = meanc(vSqrd_error2);

	savemat("vAE_RK_GARCH.xls", vAbs_error_RK);
//	adMAE_OC[0] = meanc(vSqrd_error);
		
	savemat("vSE_RK_GARCH.xls", vSqrd_error_RK);
//	adMSE_OC[0] = meanc(vSqrd_error2);

	print("\n dMAE_CC [RV, BV, RK] = ",meanc(vAbs_error_RV)~meanc(vAbs_error_BV)~meanc(vAbs_error_RK));
	print("\n dMSE_CC [RV, BV, RK] = ",meanc(vSqrd_error_RV)~meanc(vSqrd_error_BV)~meanc(vSqrd_error_RK));
	
	return 1;

}

/*
**  Function:	Compute MSE
**
**  Input: 		adMSE_OC [adress of MAE], vReturns_1 [return series], vBenchmark [Benchmark], dC [ratio]
**
**  Output: 	1
*/

//fMSE(const adMSE_OC, const vReturns_1, const vBenchmark, const dC){
//
//	decl iWindow = 250;
//	decl iT = sizerc(vReturns_1);
//	decl vH_forecast = zeros(iWindow, 1);
//	decl vSqrd_error = zeros(iWindow, 1);
//
//	dALPHA_START = 0.1;
//	dBETA_START = 0.99;
//	dOMEGA_START = 0.01;
//	dGAMMA_START = 0.1;
//
//	for(decl j = 0; j<iWindow; j++){
//		s_vY = 	vReturns_1[j:(iT - iWindow +j)];
//
//		//initialise parametervalues
//		decl vTheta = zeros(4,1);
//		vTheta[0] = dALPHA_START;
//		vTheta[1] = dBETA_START;
//		vTheta[2] = dOMEGA_START;
//		vTheta[3] = dGAMMA_START;
//	
//		//transform parameters
//		decl vThetaStar; 
//		fTransform(&vThetaStar, vTheta);
//	
//		//Maximize the LL
//		decl dFunc;
//		decl iA;
//		iA=MaxBFGS(fLogLike_Garch, &vThetaStar, &dFunc, 0, TRUE);
//	
//		//Transform thetasStar back
//	  	fTransformBack(&vTheta, vThetaStar);
//
//		dALPHA_START = vTheta[0];
//		dBETA_START = vTheta[1];
//		dOMEGA_START = vTheta[2];
//		dGAMMA_START = vTheta[3];
//
//		vH_forecast[j] = fForecast(vThetaStar);
//		vSqrd_error[j] = (dC*vH_forecast[j] - dRATIO*vBenchmark[(iT - iWindow +j)])^2;
//
//	}
//
//	savemat("vSE_RV_GARCH.xls", vSqrd_error);
//	adMSE_OC[0] = meanc(vSqrd_error);
//
//	return 1;
//}

/*
**				MAIN PROGRAM
**
**  Purpose:	Estimate GARCH parameters alpha, beta, omega and gamma.
**
**  Input: 		dALPHA, dBETA, dOMEGA, iB, iSIZE, iSIMS, iSTEPS
**
**  Output: 	Figures
*/

main(){
	//laad SBUX returns
//	decl mData_1 = loadmat("daily_sp100.csv");
//	decl mData_2 = loadmat("ReturnsCloseToClose.csv"); 
//	decl vReturns_1 = 100*mData_1[:][1];
//	decl vReturns_2 = 100*mData_2[:][1];


//fLoadData(){
	decl mData_1;
//	mData_1 = loadmat("daily_sp100.csv");
	mData_1 = loadmat("sbux_cropped.csv");
	decl vReturns_1 = mData_1[:][6];
	vReturns_1 = reversec(vReturns_1); 
	vReturns_1 = 100*(log(vReturns_1)-log(lag0(vReturns_1,1)));
	vReturns_1 = dropr(vReturns_1,0);
//	return s_vY;
//}

//	dRATIO = (varc(vReturns_1) +varc(vReturns_2))/varc(vReturns_1);

	decl vRV = loadmat("RV.csv");
	vRV = dropr(vRV,0);
	decl vBV = loadmat("BV.csv");
	vBV = dropr(vBV,0);
	decl vRK = loadmat("RK.csv");
	vRK = dropr(vRK,0);

	//laad Dates SBUX returns
	decl vTemp_Date = mData_1[][0];
//	decl vYear 		= floor(vTemp_Date/10000);							
//	decl vMonth 	= floor((vTemp_Date-floor(vTemp_Date/10000)*10000)/100);	
//	decl vDay 		= vTemp_Date-floor(vTemp_Date/100)*100;
//	s_vDate 		= dayofcalendar(vYear, vMonth, vDay);

	dALPHA_START = 0.1;
	dBETA_START = 0.99;
	dOMEGA_START = 0.01;
	dGAMMA_START = 0.01;
	iPARS = 4;
	
	decl dAlpha_hat, dBeta_hat, dOmega_hat, dGamma_hat;
	decl vVariance_1, vVariance_2;
//	print("\nO-C");
	fEstimateGarch(vReturns_1, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dGamma_hat, &vVariance_1);
//	
//	print("\nC-C");
//	fEstimateGarch(vReturns_2, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &vVariance_2);

	//Graphs
//	SetDrawWindow("EMP_NGAS(1,1)");
//	DrawTMatrix(0, (vReturns_1~sqrt(vVariance_1))', {"S&P100"}, vTemp_Date');
////	DrawTMatrix(1, (vReturns_2~sqrt(vVariance_2))', {"Close-to-close"}, s_vDate');
//	ShowDrawWindow();

	//forecasts MAE and MSE
	//decl vBenchmark = vRV;
	//decl dMAE_OC;
	//fMAE(&dMAE_OC, vReturns_1, vBenchmark, dRATIO);
	//print("\n dMAE_OC = ",dMAE_OC);
	
//	decl dMAE_CC, dMSE_CC;
	fMAE(vReturns_1, vRV, vBV, vRK);
//	print("\n dMAE_CC = ",dMAE_CC);
//	print("\n dMSE_CC = ",dMSE_CC);
	
	//decl dMSE_OC;
	//fMSE(&dMSE_OC, vReturns_1, vBenchmark, dRATIO);
	//print("\n dMSE_OC = ",dMSE_OC);
	
//	decl dMSE_CC;
//	fMSE(&dMSE_CC, vReturns_1, vBenchmark, 1);
//	print("\n dMSE_CC = ",dMSE_CC);	
}