function [S,E,Ir,Iu,Seedc]=initialize(nl,part,C,num_ens,incidence)
num_loc=size(part,1)-1;
num_mp=size(nl,1);
S=C*ones(1,num_ens);
E=zeros(num_mp,num_ens);
Ir=zeros(num_mp,num_ens);
Iu=zeros(num_mp,num_ens);

num_times=size(incidence,2);
Seedc=zeros(num_loc,num_times);

reported=sum(incidence,2);%total incidence
reported(:,2)=(1:size(reported,1));
reported=sortrows(reported,-1);

seedlocs=reported(reported(:,1)>0,2);%cases >0

for l=1:size(seedlocs)
    seedloc=seedlocs(l);
    temp=incidence(seedloc,:);
    index=find(temp>0);
    if ~isempty(index)
        T0=index(1);%first reporting
        c=sum(temp(T0:min(T0+4,num_times)));
        Seedc(seedloc,max(1,T0-8))=c;
    end
end