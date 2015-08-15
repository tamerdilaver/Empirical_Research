
/*
**	Master Thesis Econometrics
**
**  Purpose:
**  	Estimate all t-GAS model parameters 
**
**  Date:
**    	15/08/2015
**
**  Author:
**	  	Tamer Dilaver
**
**	Supervisor:
**		Fransisco Blasques
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
static decl vSTD_NORM;				// Zt ~ N(0,1)
static decl s_vY; 					//Simulated returns
static decl s_vDate;
static decl dALPHA_START;
static decl dBETA_START;
static decl dOMEGA_START;
static decl dGAMMA_START;
static decl dLAMBDA_START;
static decl iPARS;
static decl dRATIO;

/*
**  Function:	Transform (start)parameters
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
	avThetaStar[0][4] = log(vTheta[4]-2)-log(100-vTheta[4]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adAlpha, adBeta, aOmega, adGamma, adLambda, vTheta
**
**  Output: 	1 
*/

fGetPars(const adAlpha, const adBeta, const adOmega, const adGamma, const adLambda, const vTheta){

	adAlpha[0] = exp(vTheta[0]);
	adBeta[0] = exp(vTheta[1]);
	adOmega[0] = exp(vTheta[2]);
	adGamma[0] = exp(vTheta[3]);
	adLambda[0] = 2+(100-2)*exp(vTheta[4])/(1+exp(vTheta[4]));
	return 1;
}

/*
**  Function:	Calculates average value loglikelihood for t-GAS given parameter values
**
**  Input: 		vTheta [parameter values], adFunc [adress function value], avScore [the score], amHessian [hessian matrix]
**
**  Output:		1
**
*/

fLogLike_GAS(const vTheta, const adFunc, const avScore, const amHessian){
	decl dAlpha, dBeta, dOmega, dGamma, dLambda;
	fGetPars( &dAlpha,  &dBeta, &dOmega, &dGamma, &dLambda, vTheta);

	decl dS2 = dGamma;	//	decl dS2 = dOmega/(1-dBeta);					 											//initial condition by definition
	decl vLogEta = zeros(sizerc(s_vY), 1);

	//Please make these more efficient!!! delete the log()
	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vLogEta[i] = -1/2*log(M_PI)-1/2*log(dLambda-2) -1/2*log(dS2) -log(gammafact(dLambda/2))+log(gammafact((dLambda+1)/2)) -(dLambda+1)/2*log(1+ s_vY[i]^2 / ((dLambda-2)*dS2));
			//GAS recursion
			dS2 = dOmega + dAlpha*(dLambda+3)/dLambda*((dLambda+1)/(dLambda-2)*(1+sqr( s_vY[i])/((dLambda-2)*  dS2))^(-1)*sqr( s_vY[i]) -   dS2) + dBeta*dS2;
	}																	
	
	adFunc[0] = sumc(vLogEta)/sizerc(s_vY); 									 	//Average
	return 1;
}

/*
**  Function:	Transform parameters back
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
	avTheta[0][4] = 2+(100-2)*exp(vThetaStar[4])/(1+exp(vThetaStar[4]));
	//actually need to restrict dLambda_hat between (2,100)
	//otherwise there will be no convergence for small samples that occur Gaussian
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
		Num2Derivative(fLogLike_GAS, vThetaStar, &mHessian);
		NumJacobian(fTransformBack, vThetaStar, &mJacobian);	  //numerical Jacobian
		mHessian 	= mJacobian*invert(-iN*mHessian)*mJacobian';
		vStdErrors 	= sqrt(diagonal(mHessian)');

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
	decl dAlpha,  dBeta, dOmega, dGamma, dLambda, vH;

	fGetPars(&dAlpha,  &dBeta, &dOmega, &dGamma, &dLambda, vTheta);
	
	vH = zeros(sizerc(s_vY),1);
	vH[0]= dGamma;	//vH[0]= dOmega/fabs(1-dBeta);	
	
	for(decl i = 1; i < sizerc(s_vY); i++){													   //mixed 	
		vH[i] = dOmega + dAlpha*(dLambda+3)/dLambda*((dLambda+1)/(dLambda-2)*(1+sqr(s_vY[i-1])/((dLambda-2)*vH[i-1]))^(-1)*sqr(s_vY[i-1]) - vH[i-1]) + dBeta*vH[i-1];
	}	
	
	return 	vH;


}

/*
**  Function:	Estimate t-GAS parameters
**
**  Input: 		vReturns, adAlpha_hat, adBeta_hat, adOmega_hat, adGamma_hat, adLambda_hat, avVariance 
**
**  Output: 	vTheta [estimated parametervalues]
*/
fEstimate_t_GAS(const vReturns, const adAlpha_hat, const adBeta_hat, const adOmega_hat, const adGamma_hat, const adLambda_hat, const avVariance){

	//initialise parametervalues
	decl vTheta = zeros(iPARS,1);
	vTheta[0] = dALPHA_START;
	vTheta[1] = dBETA_START;
	vTheta[2] = dOMEGA_START;
	vTheta[3] = dGAMMA_START;
	vTheta[4] = dLAMBDA_START;
	decl vThetaStart = vTheta;

	//globalize returns and vectorize true pars
	s_vY = vReturns;

	//transform parameters
	decl vThetaStar; 
	fTransform(&vThetaStar, vTheta);

	//Maximize the LL
	decl dFunc;
	decl iA;
	iA=MaxBFGS(fLogLike_GAS, &vThetaStar, &dFunc, 0, TRUE);

	//Transform thetasStar back
  	fTransformBack(&vTheta, vThetaStar);

	//return alpha, beta, omega and lambda
	adAlpha_hat[0] = vTheta[0];
	adBeta_hat[0] = vTheta[1];
	adOmega_hat[0] = vTheta[2];
	adGamma_hat[0] = vTheta[3];
	adLambda_hat[0] = vTheta[4];

	decl vSigmaStdError = fSigmaStdError(vThetaStar);
	decl vVariance = fVariance(vThetaStar);
	avVariance[0] = vVariance;
	
	print("\n",MaxConvergenceMsg(iA));
	println("\nFunctiewaarde likelihood eindwaardes:", dFunc);
	print("\nOptimale parameters met standaarderrors \n",
          	"%r", { "alpha",  "beta",  "omega", "gamma", "lambda"},
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
	decl dAlpha, dBeta, dOmega, dGamma, dLambda, vH;
	fGetPars(&dAlpha, &dBeta, &dOmega, &dGamma, &dLambda, vTheta);
	
	vH = zeros((sizerc(s_vY)+1),1);
	vH[0]= dGamma;	//	vH[0]= dOmega/fabs(1-dBeta);	
	
	for(decl i = 0; i < sizerc(s_vY); i++){			
		vH[i+1] = dOmega + dAlpha*(dLambda+3)/dLambda*((dLambda+1)/(dLambda-2)*(1+sqr(s_vY[i])/((dLambda-2)*vH[i]))^(-1)*sqr(s_vY[i]) - vH[i]) + dBeta*vH[i];
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
	dLAMBDA_START = 10;

	for(decl j = 0; j<iWindow; j++){
		s_vY = 	vReturns_1[j:(iT - iWindow +j-1)];

		//initialise parameter values
		decl vTheta = zeros(iPARS,1);
		vTheta[0] = dALPHA_START;
		vTheta[1] = dBETA_START;
		vTheta[2] = dOMEGA_START;
		vTheta[3] = dGAMMA_START;
		vTheta[4] = dLAMBDA_START;
	
		//transform parameters
		decl vThetaStar; 
		fTransform(&vThetaStar, vTheta);
	
		//Maximize the LL
		decl dFunc;
		decl iA;
		iA=MaxBFGS(fLogLike_GAS, &vThetaStar, &dFunc, 0, TRUE);
	
		//Transform thetasStar back
	  	fTransformBack(&vTheta, vThetaStar);

		dALPHA_START = vTheta[0];
		dBETA_START = vTheta[1];
		dOMEGA_START = vTheta[2];
		dGAMMA_START = vTheta[3];
		dLAMBDA_START = vTheta[4];

		vH_forecast[j] = fForecast(vThetaStar);
		//print(vH_forecast[j]~vBenchmark[(iT - iWindow +j)]);
		vAbs_error_RV[j] 	= fabs(vH_forecast[j] - vRV[(iT - iWindow +j)]);
		vSqrd_error_RV[j] =  sqr(vH_forecast[j] - vRV[(iT - iWindow +j)]);

		vAbs_error_BV[j] 	= fabs(vH_forecast[j] - vBV[(iT - iWindow +j)]);
		vSqrd_error_BV[j] =  sqr(vH_forecast[j] - vBV[(iT - iWindow +j)]);

		vAbs_error_RK[j] 	= fabs(vH_forecast[j] - vRK[(iT - iWindow +j)]);
		vSqrd_error_RK[j] =  sqr(vH_forecast[j] - vRK[(iT - iWindow +j)]);

	}

	savemat("vAE_RV_tGAS.xls", vAbs_error_RV);
//	adMAE_OC[0] = meanc(vSqrd_error);
		
	savemat("vSE_RV_tGAS.xls", vSqrd_error_RV);
//	adMSE_OC[0] = meanc(vSqrd_error2);

	savemat("vAE_BV_tGAS.xls", vAbs_error_BV);
//	adMAE_OC[0] = meanc(vSqrd_error);
		
	savemat("vSE_BV_tGAS.xls", vSqrd_error_BV);
//	adMSE_OC[0] = meanc(vSqrd_error2);

	savemat("vAE_RK_tGAS.xls", vAbs_error_RK);
//	adMAE_OC[0] = meanc(vSqrd_error);
		
	savemat("vSE_RK_tGAS.xls", vSqrd_error_RK);
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
////	decl vTemp_returns = vReturns_1;
//	decl vH_forecast = zeros(iWindow, 1);
//	decl vSqrd_error = zeros(iWindow, 1);
//
//	dALPHA_START = 0.1;
//	dBETA_START = 0.85;
//	dOMEGA_START = 0.01;
//	dGAMMA_START = 10;
//	dLAMBDA_START = 10;
//	
//
//	for(decl j = 0; j<iWindow; j++){
//		s_vY = 	vReturns_1[j:(iT - iWindow +j)];
//
//		//initialise parametervalues
//		decl vTheta = zeros(iPARS,1);
//		vTheta[0] = dALPHA_START;
//		vTheta[1] = dBETA_START;
//		vTheta[2] = dOMEGA_START;
//		vTheta[3] = dGAMMA_START;
//		vTheta[4] = dLAMBDA_START;
//	
//		//transform parameters
//		decl vThetaStar; 
//		fTransform(&vThetaStar, vTheta);
//	
//		//Maximize the LL
//		decl dFunc;
//		decl iA;
//		iA=MaxBFGS(fLogLike_GAS, &vThetaStar, &dFunc, 0, TRUE);
//	
//		//Transform thetasStar back
//	  	fTransformBack(&vTheta, vThetaStar);
//
//		dALPHA_START = vTheta[0];
//		dBETA_START = vTheta[1];
//		dOMEGA_START = vTheta[2];
//		dGAMMA_START = vTheta[3];
//		dLAMBDA_START = vTheta[4];
//
//		vH_forecast[j] = fForecast(vThetaStar);
//		vSqrd_error[j] = (dC*vH_forecast[j] - dRATIO*vBenchmark[(iT - iWindow +j)])^2;
//
//	}
//
//	adMSE_OC[0] = meanc(vSqrd_error);
//
//	return 1;
//
//}

/*
**				MAIN PROGRAM
**
**				Estimate t-GAS parameters 
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
	dOMEGA_START = 0.001;
	dGAMMA_START = 0.1;
	dLAMBDA_START = 10;
	iPARS = 5;
	
	decl dAlpha_hat, dBeta_hat, dOmega_hat, dGamma_hat, dLambda_hat;
	decl vVariance_1, vVariance_2;
//	print("\nO-C");
	fEstimate_t_GAS(vReturns_1, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dGamma_hat, &dLambda_hat, &vVariance_1);

//	print("\nC-C");
//	fEstimate_t_GAS(vReturns_2, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dLambda_hat, &vVariance_2);

	//graphs
//	SetDrawWindow("EMP_t-GAS(1,1)");
//	DrawTMatrix(0, (vReturns_1~sqrt(vVariance_1))', {"S&P100"}, vTemp_Date');
////	DrawTMatrix(1, (vReturns_2~sqrt(vVariance_2))', {"Close-to-close"}, s_vDate');
//	ShowDrawWindow();

	//forecasts MAE and MSE
//	decl vBenchmark = vBV;
	//decl dMAE_OC;
	//fMAE(&dMAE_OC, vReturns_1, vBenchmark, dRATIO);
	//print("\n dMAE_OC = ",dMAE_OC);

	fMAE(vReturns_1, vRV, vBV, vRK);
//	decl dMAE_CC, dMSE_CC;
//	fMAE(&dMAE_CC, &dMSE_CC, vReturns_1, vBenchmark, 1);
//	print("\n dMAE_CC = ",dMAE_CC);
//	print("\n dMSE_CC = ",dMSE_CC);
//	
	//decl dMSE_OC;
	//fMSE(&dMSE_OC, vReturns_1, vBenchmark, dRATIO);
	//print("\n dMSE_OC = ",dMSE_OC);
	
//	decl dMSE_CC;
//	fMSE(&dMSE_CC, vReturns_1, vBenchmark, 1);
//	print("\n dMSE_CC = ",dMSE_CC);	
}
