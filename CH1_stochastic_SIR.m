% Gillespie simulation for the SIR process
% reaction 1 infection S + I -> 2I at rate alpha
% reaction 2 recovery I -> R at rate beta
 
clear 

set(0,                           ...
   'defaultaxesfontsize', 20,   ...
   'defaultaxeslinewidth', 1.0, ...
   'defaultlinelinewidth', 1.2, ...
   'defaultpatchlinewidth', 0.7);

bset=[20,55];  % set of beta values to simulate

for jc=1:2
    
N = 50; %initial susceptible population size
% initial infected population = 1;
 
K = 2000;  % number of trials
alpha = 1;  % without loss of generality
beta = bset(jc);   

%set the rate constants for two reactions:
c(1) = alpha;  %S + I -> 2I  %  
c(2) = beta; %I -> R

%Specify the change matrix
%C is a three by 2 matrix; three states, two possible reactions
C = [-1,1,0;0,-1,1];  %What happens to s, i and r when a reaction occurs:
%initialize the state space

s = N* ones(K,1); % start with N susceptibles
i = ones(K,1); %start each population with one infected

  
S = s ; %This will keep track of the   trajectories
I = i ;
T = zeros(K,1);  %This will track the transition times
j = 1; % number of reaction steps
rk = [];

ndx = find(i>0);  %update only the s's that are not zero
Nt = length(ndx); %=K at first
Kt = 5000; % max number of time steps to prevent very long runs without termination
while (Nt>0&j<Kt)
    j = j+1;
    h(:,1) = c(1)*s.*i;  %reaction rate 1
    h(:,2) = c(2)*i; % reaction rate 2
     
    hc = cumsum(h')'; % the cumulative sum of h
    H = sum(h')';
   
    rn = rand(Nt,2); %find a random number for the number of not yet terminated trajectories
    
    T(:,j) = T(:,j-1);  % add the current time to T
    
    for k = 1:Nt
     % update only those for which i >0
     
    T(ndx(k),j) = - log(rn(k,1))/H(ndx(k)) +T(ndx(k),j-1); % time of next reaction
    rk  = min(find(rn(k,2) <=hc(ndx(k),:)/H(ndx(k)) )); % this determines which reaction occurs
    s(ndx(k)) = s(ndx(k)) + C(rk,1); % update s and i
    i(ndx(k)) = i(ndx(k)) + C(rk,2);
    
    end 
    % save the values of the  trajectories
    S(:,j) = s ;
    I(:,j) = i ;
    ndx = find(i>0);  %check to see which trajectories are not extinct yet
    if (isempty(ndx)==1)
        Nt = 0;
    else
        Nt = length(ndx); %Nt is the number of trajectories not yet extinct
    
    end
end
 
 kj = 10; % use the first several trials to plot   sample trajectories
  
% phase portrait of some trajectories
 figure(1+5*(jc-1))
  for k = 1:kj
 plot(S(k,:) ,I(k,:) ,'linewidth',2)
 hold on
 end
 xlabel('S', 'fontsize', 20)
 ylabel('I', 'fontsize', 20)
title('Sample Trajectories','fontsize',20)
 hold off
 
% now process the data 
% scatterplot of recovery times and number of survivors
 figure(5+5*(jc-1))
 plot(S(:,end),T(:,end),'*','linewidth',2)
 xlabel('Number of Survivors','fontsize',20)
 ylabel('Recovery Time','fontsize',20)

 % create the pdf for extinction times from the data
 
 
 [NN,TT]=hist(T(:,end),50);  % histogram of the extinction times with 50 boxes
 tt = sort(TT);
   
 dt = mean(TT(2:end)-TT(1:end-1));  %timestep increment on the histogram
 
 % NN/(K*dt) is the approximate pdf for the extinction times from data
 % it is normalized to have total integral = 1
 figure(2+5*(jc-1))
 
 NNne0=find(NN>0);
 plot(TT(NNne0), NN(NNne0)/(K*dt),'*')
 xlabel('recovery time','fontsize',20)
 ylabel('pdf of recovery times','fontsize',20)
 

 % create a histogram for number of survivors
 SNb = zeros(N+1,1);
 for j = 1:K
     jdx = S(j,end);
     SNb(jdx+1) = SNb(jdx+1)+1;
 end
 xnb = [0:N];
 figure(3+5*(jc-1))
 sne0=find(SNb>0);
 plot(xnb(sne0),SNb(sne0)/K,'*')
 xlabel('Number of survivors','fontsize',20)
 ylabel('% of trials','fontsize',20)
  
  
end
