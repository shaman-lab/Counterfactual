/*********************************************************************
 * model.cpp
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [x*dimy+y] instead of [x][y]
 * Adapted from the code by Shawn Lankton (http://www.shawnlankton.com/2008/03/getting-started-with-mex-a-short-tutorial/)
 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm> 
using namespace std; 

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;

#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

void mexFunction(int nlmxhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
//declare variables
    mxArray *mxnl, *mxpart, *mxC, *mxCave, *mxS, *mxE, *mxIr, *mxIu, *mxpara, *mxbetamap, *mxalphamap;//input
    mxArray *mxnewS, *mxnewE, *mxnewIr, *mxnewIu, *mxdailyIr, *mxdailyIu;//output
    const mwSize *dims;
    double *nl, *part, *C, *Cave, *S, *E, *Ir, *Iu, *para, *betamap, *alphamap;//input
    double *newS, *newE, *newIr, *newIu, *dailyIr, *dailyIu;//output
    int num_mp, num_loc, num_para;

//associate inputs
    mxnl = mxDuplicateArray(prhs[0]);
    mxpart = mxDuplicateArray(prhs[1]);
    mxC = mxDuplicateArray(prhs[2]);
    mxCave = mxDuplicateArray(prhs[3]);
    mxS = mxDuplicateArray(prhs[4]);
    mxE = mxDuplicateArray(prhs[5]);
    mxIr = mxDuplicateArray(prhs[6]);
    mxIu = mxDuplicateArray(prhs[7]);
    mxpara = mxDuplicateArray(prhs[8]);
    mxbetamap = mxDuplicateArray(prhs[9]);
    mxalphamap = mxDuplicateArray(prhs[10]);
    
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);//number of subpopulation
    num_mp = (int)dims[0];
    dims = mxGetDimensions(prhs[1]);//number of locations
    num_loc = (int)dims[0]-1;
    dims = mxGetDimensions(prhs[8]);//number of parameters
    num_para = (int)dims[0];

//associate outputs
    mxnewS = plhs[0] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxnewE = plhs[1] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxnewIr = plhs[2] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxnewIu = plhs[3] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxdailyIr = plhs[4] = mxCreateDoubleMatrix(num_mp,1,mxREAL);
    mxdailyIu = plhs[5] = mxCreateDoubleMatrix(num_mp,1,mxREAL);

//associate pointers
    nl = mxGetPr(mxnl);
    part = mxGetPr(mxpart);
    C = mxGetPr(mxC);
    Cave = mxGetPr(mxCave);
    S = mxGetPr(mxS);
    E = mxGetPr(mxE);
    Ir = mxGetPr(mxIr);
    Iu = mxGetPr(mxIu);
    para = mxGetPr(mxpara);
    betamap = mxGetPr(mxbetamap);
    alphamap = mxGetPr(mxalphamap);
    
    newS = mxGetPr(mxnewS);
    newE = mxGetPr(mxnewE);
    newIr = mxGetPr(mxnewIr);
    newIu = mxGetPr(mxnewIu);
    dailyIr = mxGetPr(mxdailyIr);
    dailyIu = mxGetPr(mxdailyIu);
    ////////////////////////////////////////
    //do something
    default_random_engine generator((unsigned)time(NULL));
    //initialize auxillary variables
    int i, j;
    //para=Z,D,mu,theta,alpha1,alpha2,alpha3,beta0,beta1,...
    double Z=para[0],D=para[1],mu=para[2],theta=para[3];
    double dt1=(double)1/3, dt2=1-dt1;
    //daytime population,Ir,Iu,S,E,R
    vector<double> ND(num_loc),IrD(num_loc),IuD(num_loc),SD(num_loc),ED(num_loc),RD(num_loc);
    //daytime enter R, E, Iu
    vector<double> RentD(num_loc),EentD(num_loc),IuentD(num_loc);
    //nighttime population,Ir,Iu,S,E,R
    vector<double> NN(num_loc),IrN(num_loc),IuN(num_loc),SN(num_loc),EN(num_loc),RN(num_loc);
    //nighttime enter R, E, Iu
    vector<double> RentN(num_loc),EentN(num_loc),IuentN(num_loc);
    //total outgoing population
    vector<double> popleft(num_loc);
    //intermediate S, E, Iu, Ir, R
    vector<double> tempS(num_mp),tempE(num_mp),tempIu(num_mp),tempIr(num_mp),tempR(num_mp);
    //R and newR
    vector<double> R(num_mp), newR(num_mp);
    
    /////////////////////////////
    //change index in nl and part (0-based index)
    for (i=0; i<num_mp; i++)
        nl[i]=nl[i]-1;
    for (i=0; i<num_loc+1; i++)
        part[i]=part[i]-1;
    for (i=0; i<num_loc; i++){
        betamap[i]=betamap[i]-1;
        alphamap[i]=alphamap[i]-1;
    }
    //compute popleft
    for (i=0; i<num_loc; i++){
        for (j=part[i]+1; j<part[i+1]; j++){
            popleft[i]=popleft[i]+Cave[j];
        }
    }
    
    
    //assgn intermediate S, E, Iu, Ir, R
    for (i=0; i<num_mp; i++){
        R[i] = max(C[i]-S[i]-E[i]-Ir[i]-Iu[i],0.0);
        tempS[i] = S[i];
        tempE[i] = E[i];
        tempIr[i] = Ir[i];
        tempIu[i] = Iu[i];
        tempR[i] = R[i];
    }
    
    ///////////////////////////
    //daytime transmission
    //compute ND
    for (i=0; i<num_loc; i++){
        ND[i]=C[(int)part[i]];//i<-i
        for (j=part[i]+1; j<part[i+1]; j++){
            ND[i]=ND[i]+Ir[j];//reported infections (no mobility)
        }
    }
    for (i=0; i<num_loc; i++){
        for (j=part[i]+1; j<part[i+1]; j++){
            ND[(int)nl[j]]=ND[(int)nl[j]]+C[j]-Ir[j];//commuting with reported infections removed
        }
    }
    
    
    //comput IrD,IuD,SD,ED,RD
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            IrD[i]=IrD[i]+Ir[j];
            IuD[(int)nl[j]]=IuD[(int)nl[j]]+Iu[j];
            SD[(int)nl[j]]=SD[(int)nl[j]]+S[j];
            ED[(int)nl[j]]=ED[(int)nl[j]]+E[j];
            RD[(int)nl[j]]=RD[(int)nl[j]]+R[j];
        }
    }
    
    
    //compute RentD, EentD and IuentD
    for (i=0; i<num_loc; i++){
        for (j=part[i]+1; j<part[i+1]; j++){
            RentD[(int)nl[j]]=RentD[(int)nl[j]]+Cave[j]*RD[i]/(ND[i]-IrD[i]);
            EentD[(int)nl[j]]=EentD[(int)nl[j]]+Cave[j]*ED[i]/(ND[i]-IrD[i]);
            IuentD[(int)nl[j]]=IuentD[(int)nl[j]]+Cave[j]*IuD[i]/(ND[i]-IrD[i]);
        }
    }
    //compute for each subpopulation
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            //////////////////////////////
            double beta=para[(int)betamap[(int)nl[j]]];
            double alpha=para[(int)alphamap[i]];
            double Eexpr=beta*S[j]*IrD[(int)nl[j]]/ND[(int)nl[j]]*dt1;//new exposed due to reported cases
            double Eexpu=mu*beta*S[j]*IuD[(int)nl[j]]/ND[(int)nl[j]]*dt1;//new exposed due to unreported cases
            double Einfr=alpha*E[j]/Z*dt1;//new reported cases
            double Einfu=(1-alpha)*E[j]/Z*dt1;//new unreported cases
            double Erecr=Ir[j]/D*dt1;//new recovery of reported cases
            double Erecu=Iu[j]/D*dt1;//new recovery of unreported cases
            ///////////////////////////
            double ERenter=theta*dt1*(C[j]-Ir[j])/ND[(int)nl[j]]*RentD[(int)nl[j]];//incoming R
            double ERleft=theta*dt1*R[j]/(ND[(int)nl[j]]-IrD[(int)nl[j]])*popleft[(int)nl[j]];//outgoing R
            double EEenter=theta*dt1*(C[j]-Ir[j])/ND[(int)nl[j]]*EentD[(int)nl[j]];//incoming E
            double EEleft=theta*dt1*E[j]/(ND[(int)nl[j]]-IrD[(int)nl[j]])*popleft[(int)nl[j]];//outgoing E
            double EIuenter=theta*dt1*(C[j]-Ir[j])/ND[(int)nl[j]]*IuentD[(int)nl[j]];//incoming Iu
            double EIuleft=theta*dt1*Iu[j]/(ND[(int)nl[j]]-IrD[(int)nl[j]])*popleft[(int)nl[j]];//outgoing Iu
            
            ////////////////////
            //stochastic: poisson
            poisson_distribution<int> distribution1(Eexpr);
            int expr=min(distribution1(generator),(int)(S[j]*dt1));
            poisson_distribution<int> distribution2(Eexpu);
            int expu=min(distribution2(generator),(int)(S[j]*dt1));
            poisson_distribution<int> distribution3(Einfr);
            int infr=min(distribution3(generator),(int)(E[j]*dt1));
            poisson_distribution<int> distribution4(Einfu);
            int infu=min(distribution4(generator),(int)(E[j]*dt1));
            poisson_distribution<int> distribution5(Erecr);
            int recr=min(distribution5(generator),(int)(Ir[j]*dt1));
            poisson_distribution<int> distribution6(Erecu);
            int recu=min(distribution6(generator),(int)(Iu[j]*dt1));
            poisson_distribution<int> distribution7(ERenter);
            int Renter=distribution7(generator);
            poisson_distribution<int> distribution8(ERleft);
            int Rleft=min(distribution8(generator),(int)(R[j]*dt1));
            poisson_distribution<int> distribution9(EEenter);
            int Eenter=distribution9(generator);
            poisson_distribution<int> distribution10(EEleft);
            int Eleft=min(distribution10(generator),(int)(E[j]*dt1));
            poisson_distribution<int> distribution11(EIuenter);
            int Iuenter=distribution11(generator);
            poisson_distribution<int> distribution12(EIuleft);
            int Iuleft=min(distribution12(generator),(int)(Iu[j]*dt1));
            
            /////////////////////
            tempR[j]=max((int)(tempR[j]+recr+recu+Renter-Rleft),0);
            tempE[j]=max((int)(tempE[j]+expr+expu-infr-infu+Eenter-Eleft),0);
            tempIr[j]=max((int)(tempIr[j]+infr-recr),0);
            tempIu[j]=max((int)(tempIu[j]+infu-recu+Iuenter-Iuleft),0);
            dailyIr[j]=max((int)(dailyIr[j]+infr),0);
            dailyIu[j]=max((int)(dailyIu[j]+infu),0);
            tempS[j]=max((int)(C[j]-tempE[j]-tempIr[j]-tempIu[j]-tempR[j]),0);
        }
    }
    
    ////////////////////////////////
    //nighttime transmission
    //assgn new S, E, Iu, Ir, R
    for (i=0; i<num_mp; i++){
        newS[i] = tempS[i];
        newE[i] = tempE[i];
        newIr[i] = tempIr[i];
        newIu[i] = tempIu[i];
        newR[i] = tempR[i];
    }
    //compute NN
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            NN[i]=NN[i]+C[j];
        }
    }
    //comput IrN,IuN,SN,EN,RN
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            IrN[i]=IrN[i]+tempIr[j];
            IuN[i]=IuN[i]+tempIu[j];
            SN[i]=SN[i]+tempS[j];
            EN[i]=EN[i]+tempE[j];
            RN[i]=RN[i]+tempR[j];
        }
    }
    //compute RentN, EentN and IuentN
    for (i=0; i<num_loc; i++){
        for (j=part[i]+1; j<part[i+1]; j++){
            RentN[(int)nl[j]]=RentN[(int)nl[j]]+Cave[j]*RN[i]/(NN[i]-IrN[i]);
            EentN[(int)nl[j]]=EentN[(int)nl[j]]+Cave[j]*EN[i]/(NN[i]-IrN[i]);
            IuentN[(int)nl[j]]=IuentN[(int)nl[j]]+Cave[j]*IuN[i]/(NN[i]-IrN[i]);
        }
    }
    //compute for each subpopulation
    for (i=0; i<num_loc; i++){
        for (j=part[i]; j<part[i+1]; j++){
            //////////////////////////////
            double beta=para[(int)betamap[i]];
            double alpha=para[(int)alphamap[i]];
            double Eexpr=beta*tempS[j]*IrN[i]/NN[i]*dt2;//new exposed due to reported cases
            double Eexpu=mu*beta*tempS[j]*IuN[i]/NN[i]*dt2;//new exposed due to unreported cases
            double Einfr=alpha*tempE[j]/Z*dt2;//new reported cases
            double Einfu=(1-alpha)*tempE[j]/Z*dt2;//new unreported cases
            double Erecr=tempIr[j]/D*dt2;//new recovery of reported cases
            double Erecu=tempIu[j]/D*dt2;//new recovery of unreported cases
            ///////////////////////////
            double ERenter=theta*dt2*C[j]/NN[i]*RentN[i];//incoming R
            double ERleft=theta*dt2*tempR[j]/(NN[i]-IrN[i])*popleft[i];//outgoing R
            double EEenter=theta*dt2*C[j]/NN[i]*EentN[i];//incoming E
            double EEleft=theta*dt2*tempE[j]/(NN[i]-IrN[i])*popleft[i];//outgoing E
            double EIuenter=theta*dt2*C[j]/NN[i]*IuentN[i];//incoming Iu
            double EIuleft=theta*dt2*tempIu[j]/(NN[i]-IrN[i])*popleft[i];//outgoing Iu
            ////////////////////
            //stochastic: poisson
            poisson_distribution<int> distribution1(Eexpr);
            int expr=min(distribution1(generator),(int)(tempS[j]*dt2));
            poisson_distribution<int> distribution2(Eexpu);
            int expu=min(distribution2(generator),(int)(tempS[j]*dt2));
            poisson_distribution<int> distribution3(Einfr);
            int infr=min(distribution3(generator),(int)(tempE[j]*dt2));
            poisson_distribution<int> distribution4(Einfu);
            int infu=min(distribution4(generator),(int)(tempE[j]*dt2));
            poisson_distribution<int> distribution5(Erecr);
            int recr=min(distribution5(generator),(int)(tempIr[j]*dt2));
            poisson_distribution<int> distribution6(Erecu);
            int recu=min(distribution6(generator),(int)(tempIu[j]*dt2));
            poisson_distribution<int> distribution7(ERenter);
            int Renter=distribution7(generator);
            poisson_distribution<int> distribution8(ERleft);
            int Rleft=min(distribution8(generator),(int)(tempR[j]*dt2));
            poisson_distribution<int> distribution9(EEenter);
            int Eenter=distribution9(generator);
            poisson_distribution<int> distribution10(EEleft);
            int Eleft=min(distribution10(generator),(int)(tempE[j]*dt2));
            poisson_distribution<int> distribution11(EIuenter);
            int Iuenter=distribution11(generator);
            poisson_distribution<int> distribution12(EIuleft);
            int Iuleft=min(distribution12(generator),(int)(tempIu[j]*dt2));
            /////////////////////
            newR[j]=max((int)(newR[j]+recr+recu+Renter-Rleft),0);
            newE[j]=max((int)(newE[j]+expr+expu-infr-infu+Eenter-Eleft),0);
            newIr[j]=max((int)(newIr[j]+infr-recr),0);
            newIu[j]=max((int)(newIu[j]+infu-recu+Iuenter-Iuleft),0);
            dailyIr[j]=max((int)(dailyIr[j]+infr),0);
            dailyIu[j]=max((int)(dailyIu[j]+infu),0);
            newS[j]=max((int)(C[j]-newE[j]-newIr[j]-newIu[j]-newR[j]),0);
        }
    }
    /////////////////////
    return;
}
