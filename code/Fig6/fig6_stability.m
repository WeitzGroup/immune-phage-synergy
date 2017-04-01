% Plot time series of system near feasible coexistence state (if exists) with 
%(a) phage adsorption rate \phi=10^-9.5 (b) \phi=10^-10.5 and (c) \phi=10^-11.5

close all;
clear all;

%% Initialize parameters
TT=100; % Length of time series (h)
tstep=0.1; % Time step (h)
ntime=int64(1+TT/tstep); % Total number of time points
curfold=pwd;
para=struct('r',[],'KC',[],'KD',[],'phi',[],'beta',[],'omg',[], ...
    'eps',[],'alpha',[],'KI',[],'KN',[],'thres',[]);
% Threshold below which extinction is assumed
para.thres=1; % (ml^-1)
% Bacteria parameters
para.r=1; % Growth rate of bacteria (h^-1)
para.KC=1e9; % Carrying capacity of bacteria (ml^-1)
para.KD=2.2e6; % Bacteria conc. at which immune response is half as effective (ml^-1)
%bpop0=1e5; % Initial bacteria conc. (ml^-1)

bpop=zeros(ntime,1);

% Virus parameters
omg_range=ones(1,3); % Decay rate of phage (h^-1) for the 3 subfigures
phi_range=10.^[-9.5 -10.5 -11.5]; % Adsorption rate of phage (ml h^-1) for the 3 subfigures
para.beta=100; % Burst size of phage
%vpop0=1e6; % Initial phage conc. (ml^-1) when included

vpop=zeros(ntime,1);

% Immune response parameters
para.eps=8.2e-8; % Killing rate parameter of immune response (ml h^-1)
para.alpha=0.97; % Max growth rate of immune response (h^-1)
para.KI=2.4e7; % Max capacity of immune response (ml^-1)
para.KN=1e5; % Bacteria conc. when immune response growth rate is half its max (ml^-1)
imm=zeros(ntime,1);

%imm0=2.7e6; % Initial immune response (ml^-1)
imm0=para.KI;

% Calculate equilibrium densities
BP=omg_range./(para.beta.*phi_range);
PB=(para.r./phi_range).*(1-BP./para.KC);
PBI=PB-(para.eps.*para.KI./phi_range)./(1+BP./para.KD);

% Colors of curves in time series plot
colormat=[[31 159 0];[0 0 255];[159 0 197]]/255;

% Initial conditions for the 3 subfigures
r_dev=[0.9 0.01]; % Deviation ratio of initial conditions from coexistence state (if exists)
pop0_sub=[[BP(1)*r_dev(1) PBI(1)*r_dev(1) imm0];[BP(2)*r_dev(2) PBI(2)*r_dev(2) imm0];[1e8 1e10 imm0]];

    rel_tol=1e-8;
    abs_tol=1e-10.*ones(3,1);
    options = odeset('RelTol',rel_tol,'AbsTol',abs_tol);

tic
%% Plotting each subfigure

% Setting axis for plot
xmin=0; xmax=TT;
ymin=1e4; ymax=1e12;

% Handles for curves
hand=zeros(3,4);

for subfig=1:3
% Setting \omg and \phi
para.omg=omg_range(subfig);
para.phi=phi_range(subfig);
% Initial conditions
pop0=pop0_sub(subfig,:);
% Setting up and solving the set of ODEs
    infect_red=@(t,y)infection_immune_bistable(t,y,para);
    [T,Y] = ode45(infect_red,[0:tstep:TT],transpose(pop0),options);
    bpop(:)=Y(:,1);
    vpop(:)=Y(:,2);
    imm(:)=Y(:,3);
    figure(subfig);
    run('plot_tseries')
end
toc

    % Saving the figures
    leg_size=15;
    figure(1);
    h_leg=legend([hand(1,2) hand(1,1)],{'Phage','Bacteria'},'Location','Southwest');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig6a.fig')
    print('fig6a.eps','-depsc')
    figure(2);
    h_leg=legend([hand(2,2) hand(2,1)],{'Phage','Bacteria'},'Location','Southeast');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig6b.fig')
    print('fig6b.eps','-depsc')
    figure(3);
    h_leg=legend([hand(3,1) hand(3,2)],{'Bacteria','Phage'},'Location','Southeast');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig6c.fig')
    print('fig6c.eps','-depsc')