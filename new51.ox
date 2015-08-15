#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

main()
{
	decl mTEST;
//
//	mTEST = rann(2.9,1.2);
//	print(mTEST);
//
//	print(1*mTEST[9][0]);
//	print(meanc(rann(100*1000000,1)));

//	print(ceil(0.1));
//decl vRK = loadmat("RK.csv");
//decl mData_3 = loadmat("ReturnsCloseToClose.csv");
//decl mData_2 = loadmat("ReturnsOpenToClose.csv");
//decl mData_1 = loadmat("sbux_cropped.csv");
//	decl vReturns_1 = mData_1[:][6];
//	vReturns_1 = reversec(vReturns_1); 
//	vReturns_1 = log(vReturns_1)-log(lag0(vReturns_1,1));
//	//vReturns_1 = dropr(vReturns_1,0);
//
//	savemat("test_sbux.csv", vRK/100~vReturns_1);

//	mTEST = <1;2;3;4;5;6;7;8;9;10>;
//	decl iWindow = 5;
//	decl iT = sizerc(mTEST);
//	print(mTEST[iT-1]);
	
//	for(decl j = 0; j<iWindow; j++){		
//		print(mTEST[j:(iT - iWindow +j-1)]);	
//	}

//	mTEST = <2,3;4,5>;
//	if(2!=3){
//		print(mTEST[][1]);
//	}

//	mTEST = rann(5,5);
//	while(sumc(sumr(isdotnan(mTEST)))<2){
//		mTEST[5*ranu(1,1)][5*ranu(1,1)]=.NaN;
//	}
//
//
//	///FROM HERE TO
//	print(mTEST);
//	decl vTEST;
//	vTEST = meanc(mTEST); 
//	print(vTEST);
////	print(isdotnan(meanc(mTEST)));
////	print(vecindex(isdotnan(meanc(mTEST))));
//	for(decl i=0;i<sizerc(vecindex(isdotnan(meanc(mTEST))));i++){
//		//print("test\n");
//		vTEST[vecindex(isdotnan(meanc(mTEST)))[i]] = meanc(deleter(mTEST[:][vecindex(isdotnan(meanc(mTEST)))[i]]));
//	}
//
//	//HERE WITHOUT THE PRINTS
//
//	
//		print(vTEST);
////	vTEST[vecindex(isdotnan(meanc(mTEST)))] = meanc(deleter(mTEST[vecindex(isdotnan(meanc(mTEST)))]));
////	print(quantilec(mTEST,0.5));

	mTEST = rann(5,5);

	decl vMeanTEST;

	//vMeanTEST = meanc(mTEST[][]);
	print(mTEST'[][1]);
	
}
