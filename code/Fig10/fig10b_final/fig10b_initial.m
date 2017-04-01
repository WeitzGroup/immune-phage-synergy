% Simulate the phage-bacteria-immune model, calculate steady state
% bacteria and phage densities at different initial conditions
% and plot the resulting heat maps

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
para.omg=1;
para.phi=10^-10.5;

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
%bpop0=BIS;
%vpop0=1e5;
%vpop0=10*BIS;
ipop0=para.KI;

% Scanning parameters for phage decay rate \omega and adsorption rate \phi
B0_min=4;B0_max=12;B0_step=0.1;
B0_range=10.^(B0_min:B0_step:B0_max);
P0_min=4;P0_max=12;P0_step=0.1;
P0_range=10.^(P0_min:P0_step:P0_max);
len_B0=length(B0_range); len_P0=length(P0_range);

% Arrays to store the bacteria, phage and immune steady state densities
lastbpop=zeros(len_P0,len_B0);
lastvpop=zeros(len_P0,len_B0);
lastimm=zeros(len_P0,len_B0);

% Setting up the set of ODEs
infect_red=@(t,y)infection_immune_bistable(t,y,para);
rel_tol=1e-6;
abs_tol=1e-7.*ones(3,1);
options = odeset('RelTol',rel_tol,'AbsTol',abs_tol);

tic
% Begin looping through phage parameters
for B0c=1:len_B0
B0c
bpop0=B0_range(B0c);
for P0c=1:len_P0
vpop0=P0_range(P0c);

% Initialize simulation
pop0=[bpop0 vpop0 ipop0];
    
% Solving the set of ODEs
    [T,Y] = ode45(infect_red,[0:tstep:TT],transpose(pop0),options);
    
    % Output bacteria, phage and immune densities
    bpop=Y(mtime:ntime,1);
    vpop=Y(mtime:ntime,2);
    imm=Y(mtime:ntime,3);
    
    % Implement population extinction threshold
    bpop(bpop<para.thres)=0;
    vpop(vpop<para.thres)=0;
    
    % Averaging sampled values at steady state
    lastbpop(P0c,B0c)=mean(bpop);
    lastvpop(P0c,B0c)=mean(vpop);
    lastimm(P0c,B0c)=mean(imm);
end
end
toc

% Save variables
save('data_fig10b')

% Plot heat maps
run('plot_initial')