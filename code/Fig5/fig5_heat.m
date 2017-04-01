% Simulate the phage-bacteria-immune model, calculate steady state
% bacteria and phage densities at different phage decay and adsorption rate, 
% and plot the resulting heat maps

% Figure caption: 
%Heat map showing the dependence of (a) log steady state bacteria
%density $B_S$ and (b) log steady state phage density $P_S$ on the phage decay
%rate $\omega$ and adsorption rate $\phi$. The solid line is
%$B_P=B_I^S$, the dashed line is $B_P=B_I^U$ and the dotted line is $B_P=B_I^M$. 
%The points labeled by A, B and C are sample points of the different regimes 
%with time series shown in Figs. \ref{fig6}(a), (b) and (c) respectively. 
%The initial conditions are given by $B_0=B_I^S$, $P_0=10^5$ ml$^{-1}$ and $I_0=K_I$.

close all;

%% Initialize parameters
TT=11000; tstep=0.1; % Total time of simulation and time step
ntime=int64(1+TT/tstep); % Total number of time points
mtime=int64(1+(TT-1000)/tstep); % Time point to begin sampling steady state
n_samp=ntime-mtime+1; % Number of steady state time points sampled
para=struct('r',[],'KC',[],'KD',[],'phi',[],'beta',[],'omg',[], ...
    'eps',[],'alpha',[],'KI',[],'KN',[],'thres',[]);

% Threshold below which extinction is assumed
para.thres=1;

% Bacteria parameters
para.r=1; % Growth rate of bacteria (h^-1)
para.KC=1e9; % Carrying capacity of bacteria (ml^-1)
para.KD=2.2e6; % Bacteria conc. at which immune response is half as effective (ml^-1)

% Virus parameters
para.beta=100; % Burst size of phage

% Immune response parameters
para.eps=8.2e-8; % Killing rate parameter of immune response (ml h^-1)
para.alpha=0.97; % Max growth rate of immune response (h^-1)
para.KI=2.4e7; % Max capacity of immune response (ml^-1)
para.KN=1e5; % Bacteria conc. when immune response growth rate is half its max (ml^-1)

% Calculate equilibrium densities B_I^U, B_I^S, and B_I^M as defined in Eqs.
% (4) and (10)
disc=sqrt((para.KC+para.KD)^2/4-para.KC*para.KD*para.eps*para.KI/para.r);
BIU=(para.KC-para.KD)/2-disc;
BIS=(para.KC-para.KD)/2+disc;
BIM=para.KD*(sqrt(para.KC*para.eps*para.KI/(para.KD*para.r))-1);

% Initial conditions
bpop0=BIS;
%vpop0=1e5;
vpop0=10*BIS;
ipop0=para.KI;

% Scanning parameters for phage decay rate \omega and adsorption rate \phi
omg_min=-1;omg_max=2;omg_step=0.1;
omg_range=10.^(omg_min:omg_step:omg_max);
phi_min=-12;phi_max=-8;phi_step=0.1;
phi_range=10.^(phi_min:phi_step:phi_max);
len_omg=length(omg_range); len_phi=length(phi_range);

% Arrays to store the bacteria, phage and immune steady state densities
lastbpop=zeros(len_phi,len_omg);
lastvpop=zeros(len_phi,len_omg);
lastimm=zeros(len_phi,len_omg);

tic
% Begin looping through phage parameters
for omgc=1:len_omg
omgc
para.omg=omg_range(omgc);
for iphic=1:len_phi
para.phi=phi_range(iphic);

% Initialize simulation
pop0=[bpop0 vpop0 ipop0];
    
% Solving the set of ODEs
    infect_red=@(t,y)infection_immune_bistable(t,y,para);
    rel_tol=1e-6;
    abs_tol=1e-7.*ones(3,1);
    options = odeset('RelTol',rel_tol,'AbsTol',abs_tol);
    [T,Y] = ode45(infect_red,[0:tstep:TT],transpose(pop0),options);
    
    % Output bacteria, phage and immune densities
    bpop=Y(mtime:ntime,1);
    vpop=Y(mtime:ntime,2);
    imm=Y(mtime:ntime,3);
    
    % Implement population extinction threshold
    bpop(bpop<para.thres)=0;
    vpop(vpop<para.thres)=0;
    
    % Averaging sampled values at steady state
    lastbpop(iphic,omgc)=mean(bpop);
    lastvpop(iphic,omgc)=mean(vpop);
    lastimm(iphic,omgc)=mean(imm);
end
end
toc

% Save variables
save('data_fig5')

% Plot heat maps
run('plot_heat')