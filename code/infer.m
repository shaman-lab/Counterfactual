function infer()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for Differential Effects of Intervention Timing on COVID-19
% Spread in the United States.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The transmission model is run in cpp.
% need to mex the cpp function before use using the following command
% mex model_eakf.cpp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Td=15.4;%average reporting delay
a=2.6;%shape parameter of gamma distribution
b=Td/a;%scale parameter of gamma distribution
rnds=ceil(gamrnd(a,b,1e4,1));%pre-generate gamma random numbers
Td_death=15.4+7;%average delay of death
a=2.6;%shape parameter of gamma distribution
b=Td_death/a;%scale parameter of gamma distribution
rnds_d=ceil(gamrnd(a,b,1e4,1));%pre-generate gamma random numbers
load commutedata
load population
%%%%%%%%%%%%%%%
%Inter-county commuting is stored using neighbor list: nl, part and C. 
%nl (neighbor list) and part (partition) both describe the network structure. 
%For instance, the neighbors of location i are nl(part(i):part(i+1)-1). 
%This neighbor set include i itself. 
%C has the same size with nl. 
%The commuters from location i to location j=nl(part(i)+x) is C(part(i)+x).
%%%%%%%%%%%%%%%
num_loc=size(part,1)-1;%number of counties
num_mp=size(nl,1);%number of subpopulations
load dailyincidence%county-level incidence data unitl May 19
load dailydeaths%county-level death data unitl May 19
num_times=size(dailyincidence,2);%total length of data
obs_case=zeros(size(dailyincidence));
obs_death=zeros(size(dailydeaths));
%smooth the data: 7 day moving average
for l=1:num_loc
    for t=1:num_times
        obs_case(l,t)=mean(dailyincidence(l,max(1,t-3):min(t+3,num_times)));
        obs_death(l,t)=mean(dailydeaths(l,max(1,t-3):min(t+3,num_times)));
    end
end
load deathrate %IFR for each county, averaged based on age structure
T=73;%inference until May 3
incidence=dailyincidence(:,1:T);%incidence for seeding
%set OEV
OEV_case=zeros(size(dailyincidence));
OEV_death=zeros(size(dailydeaths));
for l=1:num_loc
    for t=1:num_times
        obs_ave=mean(obs_case(l,max(1,t-6):t));
        OEV_case(l,t)=max(25,obs_ave^2/100);
        death_ave=mean(obs_death(l,max(1,t-6):t));
        OEV_death(l,t)=max(25,death_ave^2/100);
    end
end
%adjusting inter-county movement
load MI_inter
%adjusting mobility starting from March 16, day 25
MI_inter=MI_inter(:,2:end);%remove fips code
MI_inter(:,2:end)=min(MI_inter(:,2:end),1);
MI_inter_relative=[ones(num_loc,9),MI_inter];
for t=16:size(MI_inter,2)
    MI_inter_relative(:,t+9)=MI_inter(:,t)./MI_inter(:,t-1);
    MI_inter_relative(isnan(MI_inter_relative(:,t+9)),t+9)=0;
    MI_inter_relative(isinf(MI_inter_relative(:,t+9)),t+9)=0;
    MI_inter_relative(:,t+9)=min(MI_inter_relative(:,t+9),1.01);
end
MI_inter_relative(:,1:24)=1;

C=C*ones(1,num_times);
Cave=Cave*ones(1,num_times);
%adjusting mobility starting from March 16, day 25
for t=25:num_times%day 25: March 16
    C(:,t)=C(:,t-1);
    Cave(:,t)=Cave(:,t-1);
    for l=1:num_loc
        for j=part(l)+1:part(l+1)-1
            if t<=size(MI_inter_relative,2)
                C(part(l),t)=C(part(l),t)+((1-MI_inter_relative(nl(j),t))*C(j,t));
                Cave(part(l),t)=Cave(part(l),t)+((1-MI_inter_relative(nl(j),t))*Cave(j,t));
                C(j,t)=(MI_inter_relative(nl(j),t)*C(j,t));
                Cave(j,t)=(MI_inter_relative(nl(j),t)*Cave(j,t));
            end
        end
    end
end
num_ens=100;%number of ensemble members
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[S,E,Ir,Iu,Seedc]=initialize(nl,part,C(:,1),num_ens,incidence);
%%%%%%%%%%%%%%%%%%%%
%S,E,Ir and Iu represent susceptible, exposed, reported infection,
%unreported infection in all subpopulations.
%%%%%%%%%%%%%%%%%%%%
obs_temp=zeros(num_loc,num_ens,num_times);%records of reported cases
death_temp=zeros(num_loc,num_ens,num_times);%records of death
load('parafit1')% prior parameter setting, estimated using data until March 13 in a previous work
%initialize parameters
[para,paramax,paramin,betamap,alphamaps]=initializepara_eakf(dailyincidence,num_ens,parafit);
paramax_ori=paramax;
paramin_ori=paramin;
%%%%%%%%%%%%inflate parameters
lambda=2;
para=mean(para,2)*ones(1,num_ens)+lambda*(para-mean(para,2)*ones(1,num_ens));
para=checkbound_para(para,paramax,paramin);
%fix Z, D and theta
%Z
para(1,:)=parafit(3,1:num_ens);
%D
para(2,:)=parafit(4,1:num_ens);
%mu
para(3,:)=parafit(2,1:num_ens);
%theta
para(4,:)=parafit(6,1:num_ens);
parastd=std(para,0,2);
%EAKF
lambda=1.2;%inflation to avoid ensemble collapse
num_para=size(para,1);
para_post=zeros(num_para,num_ens,T);%posterior parameters
S_post=zeros(num_loc,num_ens,T);%posterior susceptiblilty
for t=1:T%start from Feb 21
    t
    tic
    %smoothing constraint for parameters
    for i=1:num_para
        if i>4
            paramax(i)=min(mean(para(i,:))*1.3,paramax_ori(i));
            paramin(i)=max(mean(para(i,:))*0.7,paramin_ori(i));
        end
    end
    %fix Z, D and theta
    %Z
    para(1,:)=parafit(3,1:num_ens);
    %D
    para(2,:)=parafit(4,1:num_ens);
    %mu
    para(3,:)=parafit(2,1:num_ens);
    %theta
    para(4,:)=parafit(6,1:num_ens);
    %%%%%%%%%%%seeding
    if t<=size(Seedc,2)
        [S,E,Ir,Iu]=seeding(S,E,Ir,Iu,nl,part,C(:,t),Seedc,t);
    end
    %%%%%%%%%%%%%%%%%%
    %record model variables
    S_temp=S; E_temp=E; Ir_temp=Ir; Iu_temp=Iu;
    %integrate forward one step
    dailyIr_prior=zeros(num_mp,num_ens);
    dailyIu_prior=zeros(num_mp,num_ens);
    for k=1:num_ens%run for each ensemble member
        %adjust population according to change of inter-county movement
        [S(:,k),E(:,k),Ir(:,k),Iu(:,k)]=adjustmobility(S(:,k),E(:,k),Ir(:,k),Iu(:,k),nl,part,MI_inter_relative,t);
        [S(:,k),E(:,k),Ir(:,k),Iu(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t),Cave(:,t),S(:,k),E(:,k),Ir(:,k),Iu(:,k),para(:,k),betamap,alphamaps);
        dailyIr_prior(:,k)=dailyIr_temp;
        dailyIu_prior(:,k)=dailyIu_temp;
    end
    %%%%%%%%%%%%%%%%%%%%%%
    %mini-forecast based on current model state
    %integrate forward for 16 days, prepare for observation
    Tproj=min(14,num_times-t);
    obs_temp1=obs_temp;
    death_temp1=death_temp;
    for t1=t:t+Tproj
        for k=1:num_ens%run for each ensemble member
            [S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k)]=adjustmobility(S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),nl,part,MI_inter_relative,t1);
            [S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),dailyIr_temp,dailyIu_temp]=model_eakf(nl,part,C(:,t1),Cave(:,t1),S_temp(:,k),E_temp(:,k),Ir_temp(:,k),Iu_temp(:,k),para(:,k),betamap,alphamaps);
            %reporting delay
            for l=1:num_loc
                for j=part(l):part(l+1)-1
                    inci=round(dailyIr_temp(j));
                    if inci>0
                        rnd=datasample(rnds,inci);
                        for h=1:length(rnd)
                            if (t1+rnd(h)<=num_times)
                                obs_temp1(l,k,t1+rnd(h))=obs_temp1(l,k,t1+rnd(h))+1;
                            end
                        end 
                    end
                end
            end
            %death delay
            for l=1:num_loc
                for j=part(l):part(l+1)-1
                    inci=round((dailyIr_temp(j)+dailyIu_temp(j))*deathrate(l));
                    if inci>0
                        rnd=datasample(rnds_d,inci);
                        for h=1:length(rnd)
                            if (t1+rnd(h)<=num_times)
                                death_temp1(l,k,t1+rnd(h))=death_temp1(l,k,t1+rnd(h))+1;
                            end
                        end 
                    end
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update model using death 16 days ahead
    t1=min(t+Tproj,num_times);
    death_ens=death_temp1(:,:,t1);%death at t1, prior
    %loop through local observations
    for l=1:num_loc
        %%%%%%%%%%%%%%%%%%%death
        %Get the variance of the ensemble
        obs_var = OEV_death(l,t1);
        prior_var = var(death_ens(l,:));
        post_var = prior_var*obs_var/(prior_var+obs_var);
        if prior_var==0%if degenerate
            post_var=1e-3;
            prior_var=1e-3;
        end
        prior_mean = mean(death_ens(l,:));
        post_mean = post_var*(prior_mean/prior_var + obs_death(l,t1)/obs_var);
        %%%% Compute alpha and adjust distribution to conform to posterior moments
        alpha = (obs_var/(obs_var+prior_var)).^0.5;
        dy = post_mean + alpha*(death_ens(l,:)-prior_mean)-death_ens(l,:);
        %Loop over each state variable (connected to location l)
        %adjust related metapopulation
        neighbors=part(l):part(l+1)-1;%metapopulation live in l
        neighbors1=find(nl==l);%%metapopulation work in l
        neighbors=union(neighbors1,neighbors);
        for h=1:length(neighbors)
            j=neighbors(h);
            %E
            temp=E(j,:);
            A=cov(temp,death_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            E(j,:)=E(j,:)+dx;
            %Ir
            temp=Ir(j,:);
            A=cov(temp,death_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            Ir(j,:)=Ir(j,:)+dx;
            %Iu
            temp=Iu(j,:);
            A=cov(temp,death_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            Iu(j,:)=Iu(j,:)+dx;
            %dailyIr
            temp=dailyIr_prior(j,:);
            A=cov(temp,death_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            dailyIr_prior(j,:)=(max(dailyIr_prior(j,:)+dx,0));
            %dailyIu
            temp=dailyIu_prior(j,:);
            A=cov(temp,death_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            dailyIu_prior(j,:)=(max(dailyIu_prior(j,:)+dx,0));
        end
        %adjust alpha before running out of observations of case and death
        if (t+10<num_times)
            temp=para(alphamaps(l),:);
            A=cov(temp,death_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            para(alphamaps(l),:)=para(alphamaps(l),:)+dx;
            %inflation
            if std(para(alphamaps(l),:))<parastd(alphamaps(l))
                para(alphamaps(l),:)=mean(para(alphamaps(l),:),2)*ones(1,num_ens)+lambda*(para(alphamaps(l),:)-mean(para(alphamaps(l),:),2)*ones(1,num_ens));
            end
        end
        %adjust beta
        temp=para(betamap(l),:);
        A=cov(temp,death_ens(l,:));
        rr=A(2,1)/prior_var;
        dx=rr*dy;
        para(betamap(l),:)=para(betamap(l),:)+dx;
        %inflation
        if std(para(betamap(l),:))<parastd(betamap(l))
            para(betamap(l),:)=mean(para(betamap(l),:),2)*ones(1,num_ens)+lambda*(para(betamap(l),:)-mean(para(betamap(l),:),2)*ones(1,num_ens));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update model using case 9 days ahead
    t1=min(t+10,num_times);
    obs_ens=obs_temp1(:,:,t1);%observation at t1, prior
    for l=1:num_loc
        %%%%%%%%%%%%%%%%%%%case
        %Get the variance of the ensemble
        obs_var = OEV_case(l,t1);
        prior_var = var(obs_ens(l,:));
        post_var = prior_var*obs_var/(prior_var+obs_var);
        if prior_var==0%if degenerate
            post_var=1e-3;
            prior_var=1e-3;
        end
        prior_mean = mean(obs_ens(l,:));
        post_mean = post_var*(prior_mean/prior_var + obs_case(l,t1)/obs_var);
        %%%% Compute alpha and adjust distribution to conform to posterior moments
        alpha = (obs_var/(obs_var+prior_var)).^0.5;
        dy = post_mean + alpha*(obs_ens(l,:)-prior_mean)-obs_ens(l,:);
        %Loop over each state variable (connected to location l)
        %adjust related metapopulation
        neighbors=part(l):part(l+1)-1;%metapopulation live in l
        neighbors1=find(nl==l);%%metapopulation work in l
        neighbors=union(neighbors1,neighbors);
        for h=1:length(neighbors)
            j=neighbors(h);
            %E
            temp=E(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            E(j,:)=E(j,:)+dx;
            %Ir
            temp=Ir(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            Ir(j,:)=Ir(j,:)+dx;
            %Iu
            temp=Iu(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            Iu(j,:)=Iu(j,:)+dx;
            %dailyIr
            temp=dailyIr_prior(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            dailyIr_prior(j,:)=(max(dailyIr_prior(j,:)+dx,0));
            %dailyIu
            temp=dailyIu_prior(j,:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            dailyIu_prior(j,:)=(max(dailyIu_prior(j,:)+dx,0));
        end
        %adjust alpha
        if (t+10<num_times)
            temp=para(alphamaps(l),:);
            A=cov(temp,obs_ens(l,:));
            rr=A(2,1)/prior_var;
            dx=rr*dy;
            para(alphamaps(l),:)=para(alphamaps(l),:)+dx;
            %inflation
            if std(para(alphamaps(l),:))<parastd(alphamaps(l))
                para(alphamaps(l),:)=mean(para(alphamaps(l),:),2)*ones(1,num_ens)+lambda*(para(alphamaps(l),:)-mean(para(alphamaps(l),:),2)*ones(1,num_ens));
            end
        end
        %adjust beta
        temp=para(betamap(l),:);
        A=cov(temp,obs_ens(l,:));
        rr=A(2,1)/prior_var;
        dx=rr*dy;
        para(betamap(l),:)=para(betamap(l),:)+dx;
        %inflation
        if std(para(betamap(l),:))<parastd(betamap(l))
            para(betamap(l),:)=mean(para(betamap(l),:),2)*ones(1,num_ens)+lambda*(para(betamap(l),:)-mean(para(betamap(l),:),2)*ones(1,num_ens));
        end
    end
    para=checkbound_para(para,paramax,paramin);
    %update Ir and Iu posterior
    dailyIr_post=dailyIr_prior;
    dailyIu_post=dailyIu_prior;
    %fix alpha after running out of observations
    if (t+10>=num_times)
        para(5,:)=para_post(5,:,t-1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %update observations
    for k=1:num_ens
        %update obs_temp
        for l=1:num_loc
            for j=part(l):part(l+1)-1
                inci=round(dailyIr_post(j,k));
                if inci>0
                    rnd=datasample(rnds,inci);
                    for h=1:length(rnd)
                        if (t+rnd(h)<=num_times)
                            obs_temp(l,k,t+rnd(h))=obs_temp(l,k,t+rnd(h))+1;
                        end
                    end 
                end
            end
        end
        %update death_temp
        for l=1:num_loc
            for j=part(l):part(l+1)-1
                inci=round((dailyIr_post(j,k)+dailyIu_post(j,k))*deathrate(l));
                if inci>0
                    rnd=datasample(rnds_d,inci);
                    for h=1:length(rnd)
                        if (t+rnd(h)<=num_times)
                            death_temp(l,k,t+rnd(h))=death_temp(l,k,t+rnd(h))+1;
                        end
                    end 
                end
            end
        end
    end
    [S,E,Ir,Iu]=checkbound(S,E,Ir,Iu);
    para_post(:,:,t)=para;
    for l=1:num_loc
        S_post(l,:,t)=min(1,sum(S(part(l):part(l+1)-1,:))/population(l));
    end
   toc
end

%save posterior parameter (para_post).
save('PosteriorParameters.mat','para_post');


function para = checkbound_para(para,paramax,paramin)
for i=1:size(para,1)
    if (i>4)
        para(i,para(i,:)<paramin(i))=paramin(i)*(1+0.1*rand(sum(para(i,:)<paramin(i)),1));
        para(i,para(i,:)>paramax(i))=paramax(i)*(1-0.1*rand(sum(para(i,:)>paramax(i)),1));
    end
end

function [S,E,Ir,Iu]=checkbound(S,E,Ir,Iu)
for k=1:size(S,2)
    S(S(:,k)<0,k)=0; E(E(:,k)<0,k)=0; Ir(Ir(:,k)<0,k)=0; Iu(Iu(:,k)<0,k)=0;
end

function [S,E,Ir,Iu]=adjustmobility(S,E,Ir,Iu,nl,part,MI_inter_relative,t)
num_loc=size(MI_inter_relative,1);
for l=1:num_loc
    for j=part(l)+1:part(l+1)-1
        if t<=size(MI_inter_relative,2)
            S(part(l))=S(part(l))+((1-MI_inter_relative(nl(j),t))*S(j));
            S(j)=(MI_inter_relative(nl(j),t)*S(j));
            E(part(l))=E(part(l))+((1-MI_inter_relative(nl(j),t))*E(j));
            E(j)=(MI_inter_relative(nl(j),t)*E(j));
            Ir(part(l))=Ir(part(l))+((1-MI_inter_relative(nl(j),t))*Ir(j));
            Ir(j)=(MI_inter_relative(nl(j),t)*Ir(j));
            Iu(part(l))=Iu(part(l))+((1-MI_inter_relative(nl(j),t))*Iu(j));
            Iu(j)=(MI_inter_relative(nl(j),t)*Iu(j));
        end
    end
end