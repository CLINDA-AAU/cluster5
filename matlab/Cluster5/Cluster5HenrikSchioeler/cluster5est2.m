clear

Pp=1/70;
Ns=4*350000; %number of samples 
p=12/Ns; %estimate of cluster5 probability
V=p*(1-p)/Ns; %variance of estimate of cluster5 probabilitysqrt
Ninit=ceil(p*6e6);
Nlow=ceil(Ninit-2*sqrt(V)*6e6);
Nhigh=ceil(Ninit+2*sqrt(V)*6e6);

R=1.2; %reproduction no (/5 days)
gamma=1/5;
beta=R*gamma;
Rw=R^(7/5); % repr no (/week)


M=20; %no of weeks since last observation of cluster5
Np=85;
Probs=zeros(Np,M);
Probs(Nlow:Nhigh,1)=1/(Nhigh-Nlow);
Pint1=zeros(Np,1);



% p_est=p;
% V_est=V;
clf

for i=2:M
    Pint1=zeros(Np,1);

   
    for k=1:Np
        %k
        for j=1:k
            Pint1(k)=Pint1(k)+poisspdf(k-j,7*beta*j)*Probs(j,i-1);
        end
    end
    
    for k=1:Np
        %k
        for j=k:Np
            Probs(k,i)=Probs(k,i)+poisspdf(j-k,7*gamma*j)*Pint1(j);
        end
    end
    
    for k=1:Np
        Pc=(k-1)/6e6;
        Probs(k,i)=Probs(k,i)*(1-Pc/Pp)^(1500);
    end
    
    Probs(:,i)=Probs(:,i)/sum(Probs(:,i));
    plot(0:Np-1,Probs(:,i))
    hold on
    Probs(1,i)
    pause(0.1)
end
S=sum(Probs(1:11,:));
clf
plot(S)
hold on
plot(Probs(1,:))








