%BirthDeadBayesFilter calculates the Bayes filter solution for estimating the discrete
% stokastic variable number of infected with a rare virus like Cluster5 Corona
%
% External input: None

% Time-stamp: <2020-11-16 07:44:48 tk>
% Version 1: 2020-11-13 13:28:25 tk Intial version based on the ideas from Henrik
%            Schiøler (HS)
% Torben Knudsen
% Aalborg University, Dept. of Electronic Systems, Section of Automation
% and Control
% E-mail: tk@es.aau.dk

%% Parameters
% Switch for using HS parameters etc to compare with he's results
HS= 1;
% TK test parameters
N= 300000; n= 10000; gamma= 1/3.4; R= 1.2;
MaxI= 100; NumDays= 4*7; IniMean= 12; UseMeas= 1; fig= 1;
TimeIntervention= inf;
n= ones(NumDays,1)*n;                   % Test size at day 1:NumDays
IniProb= 3;
% HS test parameters
if HS
  Pp=1.5e-2; %CoVid19 frequency among tests
  Ns=4*350000; %approximate number of samples from mid Aug to mid Sept
  p=12/Ns; %mean estimate of cluster5 probability
  V=p*(1-p)/Ns; %variance of estimate of cluster5 probability
  Ninit=ceil(p*6e6); %mean estimate of number of Cluster5 among population
  Nlow=ceil(Ninit-2*sqrt(V)*6e6); %low limit
  Nhigh=ceil(Ninit+2*sqrt(V)*6e6); %high limit
  N= 6e6*Pp; gamma= 1/5; R= 1.2;
  MaxI= 85; NumDays= 7*20; IniMean= 12; UseMeas= 1; fig= 1;
  MaxI= 85; NumDays= 60;   IniMean= 12; UseMeas= 1; fig= 1;
  TimeIntervention= 7*4; n1= round(1500/7); n2= round(3500/7); 
  n= ones(NumDays,1)*n1;
  n(TimeIntervention:end)= n2;
  IniProb= 2;
end

%% Definitions etc.
beta= R*gamma;
% Construct Q
Q1= diag((1:MaxI-1)*gamma,1);
Q1= Q1+diag((0:MaxI-2)*beta,-1);
Q1= Q1-diag(sum(Q1));
EQ1= expm(Q1);
% HS Q2
if HS
  Q2= diag((1:MaxI-1)*gamma*2,1);
  Q2= Q2+diag((0:MaxI-2)*beta/2,-1);
  Q2= Q2-diag(sum(Q2));
  EQ2= expm(Q2);
end;

% Measurements
Y= zeros(NumDays,1);                    % Assume no Cluster5
Y(20)= 0*1;                               % Inclue a measurements of Cluster5
if ~UseMeas;
  Y= Y*nan;
end;

% Initial probability for states
P= zeros(MaxI,1);
switch IniProb
  case 1                                % Kronecker delta
    P(IniMean)= 1;                 
  case 2                                % Uniform
    P(Nlow:Nhigh)=1/(Nhigh-Nlow+1);  
  case 3                                % Poisson
    P= poisspdf((1:MaxI)',IniMean); 
end
PP= zeros(MaxI,NumDays);

for i= 1:NumDays
  % Measurement update
  if ~isnan(Y(i))                       % If there is a measurement
    P= P.*binopdf(Y(i)*ones(MaxI,1),n(i)*ones(MaxI,1),((1:MaxI)'-1)/N);
    P= P/sum(P);
  end;
  PP(:,i)= P;                           % PP does not include the initial PDF
                                        % Time update;
  if i<TimeIntervention
    EQ= EQ1;
  else
    EQ= EQ2;
  end
  P= EQ*P;
  % Correction for the error due to truncating number of infected to MaxI
  P= P/sum(P); 
end;

figure(fig); fig= fig+1;
th= tiledlayout('flow');
nexttile
plot((1:size(PP,1))-1,PP);
xlabel('#C5');
ylabel('PDF');
grid('on');
nexttile
C= (1:NumDays)'*ones(1,MaxI);
mesh((1:size(PP,1))-1,1:size(PP,2),PP',C);
set(gca,'View',[30 30]);
xlabel('#C5');
ylabel('Time');
zlabel('PDF');
nexttile
plot(PP(1,:));
xlabel('Time');
ylabel('Prob. #C5=0');
grid('on')

if 1
  figure(fig); fig= fig+1;
  C= (1:NumDays)'*ones(1,MaxI);
  mesh((1:size(PP,1))-1,1:size(PP,2),PP',C);
  set(gca,'View',[30 30]);
  xlabel('#C5');
  ylabel('Time');
  zlabel('PDF');
end;

% Numerical results
pp= [0.9 0.95 0.99 0.999];
pp= pp(find(pp<max(PP(1,:))));          % Limit to actuap max P(C5=0)
AuxRes= [];
for i= 1:length(pp)
  d= find(PP(1,:)>pp(i));
  AuxRes= [AuxRes; pp(i) PP(1,min(d)) min(d)];
end;
ColLab= {'Min probability of zero Cluster5' 'Probability of zero Cluster5' 'Day'};
Lab= 'Probability of zero Cluster5 by day';
disp(Lab);
disp(array2table(AuxRes,'VariableNames',ColLab));

% Save plots
if 0
  print -depsc -f1 C5PlotsTK;
  print -depsc -f2 PDFC5TK;
end
if 0
  print -depsc -f1 C5PlotsTK60Days;
  print -depsc -f2 PDFC5TK60Days;
end;
