function [para,paramax,paramin,betamap,alphamap]=initializepara_eakf(dailyincidence,num_ens,parafit)
%Z,D,mu,theta,alpha,beta1,...,beta400,beta401,...,beta416
Zlow=2;Zup=5;%latency period
Dlow=2;Dup=5;%infectious period
mulow=0.1;muup=1.0;%relative transmissibility
thetalow=0.01;thetaup=0.3;%movement factor
alphalow=0.03;alphaup=0.25;%reporting rate
betalow=0.01;betaup=2.5;%transmission rate

%selected areas
num_loc=size(dailyincidence,1);
ranks=zeros(num_loc,2);
ranks(:,1)=(1:num_loc)';
ranks(:,2)=sum(dailyincidence,2);
ranks=sortrows(ranks,-2);

%define alpha
selected=ranks(ranks(:,2)>=400,1);
alphamap=ones(num_loc,1)*5;%share one alpha

%define beta
betamap=zeros(num_loc,1);%beta0 for low-case counties
betamap(selected)=5+(1:length(selected))';

load popd
PD=log10(popd);
PD_median=median(PD);

factors=ones(length(selected)+16,1);
scale=0.8;
for i=1:length(selected)
    loc=selected(i);
    factors(i)=PD(loc)/PD_median*scale;
end

notselected=(1:num_loc)';
notselected(selected)=[];
total=sum(dailyincidence,2);
totalrank(:,1)=notselected;
totalrank(:,2)=total(notselected);
totalrank=sortrows(totalrank,-2);
PDrank(:,1)=notselected;
PDrank(:,2)=PD(notselected);
PDrank=sortrows(PDrank,-2);

for i=1:length(notselected)
    l=notselected(i);
    caserankid=find(totalrank(:,1)==l);
    caserankid=ceil(caserankid/(length(notselected)/4));
    PDrankid=find(PDrank(:,1)==l);
    PDrankid=ceil(PDrankid/(length(notselected)/4));
    betamap(l)=5+length(selected)+(caserankid-1)*4+PDrankid;
    
    PDave=PDrank((PDrankid-1)*floor(length(notselected)/4)+1:PDrankid*floor(length(notselected)/4),2);
    PDave=mean(PDave);
    factors(length(selected)+(caserankid-1)*4+PDrankid)=PDave/PD_median*0.8;
    
end

%Z,D,mu,theta,alpha,beta1,...,beta400,beta401,...,beta416
paramin=[Zlow;Dlow;mulow;thetalow;alphalow;ones(length(selected)+16,1)*betalow];
paramax=[Zup;Dup;muup;thetaup;alphaup;ones(length(selected)+16,1)*betaup];

para=zeros(size(paramin,1),num_ens);

%parafit:beta;mu;Z;D;alpha;theta
%Z
para(1,:)=parafit(3,ceil(rand(1,num_ens)*size(parafit,2)));
%D
para(2,:)=parafit(4,ceil(rand(1,num_ens)*size(parafit,2)));
%mu
para(3,:)=parafit(2,ceil(rand(1,num_ens)*size(parafit,2)));
%theta
para(4,:)=parafit(6,ceil(rand(1,num_ens)*size(parafit,2)));
%alpha
para(5,:)=parafit(5,ceil(rand(1,num_ens)*size(parafit,2)));

%beta
for i=6:size(paramin,1)
    para(i,:)=parafit(1,ceil(rand(1,num_ens)*size(parafit,2)))*factors(i-5);
end
