Nit=1000; %number of iterations
R=1.2; %reproduction no (/5 days)
gamma_nom=1/5;
beta_nom=R*gamma;
beta_nom=0.2880; %nominal beta
gamma_nom=0.2; %nominal gamma
rho=0.5; %relative change in search space

Res=zeros(3,Nit);

clf
plot3(gamma_nom,beta_nom,1e-2,'*')
axis([(1-rho)*gamma_nom (1+rho)*gamma_nom (1-rho)*beta_nom (1+rho)*beta_nom 0 100])
hold on
for it=1:Nit
    beta=(rand(1)-0.5)*2*rho*beta_nom+beta_nom;
    gamma=(rand(1)-0.5)*2*rho*gamma_nom+gamma_nom;
    [likelihood OT]=cluster5est2_Markov_func_2(beta,gamma,12,4,4)
    Results(:,it)=[gamma;beta;likelihood;OT];
    if(likelihood>0.001)
        plot3(gamma,beta,OT,'.b')
    end
     if(likelihood<0.001)
       plot3(gamma,beta,OT,'.r')
    end
        
    pause(0.1)
end
    