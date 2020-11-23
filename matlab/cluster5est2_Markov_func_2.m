function [likelihood OT]=cluster5est2_Markov_func_2(beta,gamma,M,W,W2)

Pp=1.5e-2; %CoVid19 frequency among tests
Ns=4*350000; %approximate number of samples from mid Aug to mid Sept
p=12/Ns; %mean estimate of cluster5 probability
V=p*(1-p)/Ns; %variance of estimate of cluster5 probability
Ninit=ceil(p*6e6); %mean estimate of number of Cluster5 among population
Nlow=ceil(Ninit-2*sqrt(V)*6e6); %low limit
Nhigh=ceil(Ninit+2*sqrt(V)*6e6); %high limit
OT=M*7; %to ensure assignment of return value

R=1.2; %reproduction no (/5 days)
% gamma=1/5;
% beta=1.2*R*gamma;


%M=7; %no of weeks since last observation of cluster5
T_obs=7*W; %number of days from last obs of Cluster5 to lockdown and massive test
Conditional=1; %Flag to control whether results are conditional (1) or joint probabilities (0) - to calculate observation likelihood
Np=1000; %highest number of infected in model
Probs=zeros(Np,7*M); %matrix hold all results
Probs(Nlow:Nhigh,1)=1/(Nhigh-Nlow); %initial uniform distribution
Q=zeros(Np,Np); %rate matrix before lockdown and massive test
Q2=Q; %rate matrix after lockdown

eta=ones(7*M,1);

%constructing Q-matrix
for k=1:Np
    %k
    if(k>1)
        %recovery
        Q(k,k-1)=(k-1)*gamma; 
        Q2(k,k-1)=(k-1)*gamma*2;
    end
    if(k<Np)
        %infection
        Q(k,k+1)=(k-1)*beta;
        Q2(k,k+1)=(k-1)*beta/2;
    end
    %ensuring 0 row sum
    Q(k,k)=-sum(Q(k,:));
    Q2(k,k)=-sum(Q2(k,:));
end

%Solution of system of diff equations between 1 day sample points
EQ=expm(Q);
EQ2=expm(Q2);

%running the calculation over the entire timespan
for i=2:7*M   %no of weeks times days pr week
  
    %prediction model for 1 step (1 day)
    if(i<T_obs)
        Probs(:,i)=EQ'*Probs(:,i-1);
    else
        Probs(:,i)=EQ2'*Probs(:,i-1);
    end
    
    %observation update
    
    for k=1:Np
        Pc=(k-1)/6e6; %frequency in population given conditioned k-1 infected
        %Pc/Pp is then conditional probabilty of Cluster5 given Covid-19 positive
        %(1-Pc/Pp)^(N) is the probability of observing 0 Cluster5 in a sample of N Covid-19 positives which are sequenced
        if(i<T_obs)
            Probs(k,i)=Probs(k,i)*(1-Pc/Pp)^(1500/7); %before lockdown 1500 sequenced pr week
        else
            Probs(k,i)=Probs(k,i)*(1-Pc/Pp)^(3500/7); %after lockdown with 3500 sequenced pr week
        end
    end
    eta(i)=eta(i-1)*sum(Probs(:,i));
    if(i==W2*7)
        likelihood=eta(i);
    end
    if(Probs(1,i)>0.99)
        OT=i;
        break
    end
    
    %eta(i)=sum(Probs(:,i)); %recording normalization factors
    if(Conditional==1)
        Probs(:,i)=Probs(:,i)/sum(Probs(:,i));
    end

end




