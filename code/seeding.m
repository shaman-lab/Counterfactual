function [S,E,Ir,Iu]=seeding(S,E,Ir,Iu,nl,part,C,Seedc,t)
num_loc=size(Seedc,1);
num_ens=size(S,2);
for l=1:num_loc
    seedloc=l;
    seedE=(rand(1,num_ens)*20*Seedc(l,t));
    seedIu=(rand(1,num_ens)*18*Seedc(l,t));
    if Seedc(l,t)>0
        pop=sum(C(part(seedloc):part(seedloc+1)-1));
        for i=part(seedloc):part(seedloc+1)-1
            E(i,:)=(seedE*C(i)/pop);
            Iu(i,:)=(seedIu*C(i)/pop);
            S(i,:)=S(i,:)-E(i,:)-Ir(i,:)-Iu(i,:);
        end
    end
end
