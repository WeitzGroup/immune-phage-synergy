% Plot time series of system with (a) bacteria and phage only, (b) bacteria
% and immune response only and (c) bacteria, phage and immune response
% Figure caption: 
%Time series of bacteria density $B$ (green curve), phage density $P$
%(blue curve) and immune response $I$ (purple curve) in our model with (a) only
%bacteria and phage, (b) only bacteria and immune response, and (c) bacteria,
%phage and immune response combined.

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
bpop0=1e6; % Initial bacteria conc. (ml^-1)

bpop=zeros(ntime,1);

% Virus parameters
para.omg=1; % Decay rate of phage (h^-1)
para.beta=100; % Burst size of phage
para.phi=5e-8; % Adsorption rate of phage (ml h^-1)
vpop0=1e7; % Initial phage conc. (ml^-1) when included

vpop=zeros(ntime,1);

% Immune response parameters
para.eps=8.2e-8; % Killing rate parameter of immune response (ml h^-1)
para.alpha=0.97; % Max growth rate of immune response (h^-1)
para.KI=2.4e7; % Max capacity of immune response (ml^-1)
para.KN=1e5; % Bacteria conc. when immune response growth rate is half its max (ml^-1)
imm=zeros(ntime,1);

imm0=2.7e6; % Initial immune response (ml^-1) when included
%imm0=para.KI;

% Colors of curves in time series plot
colormat=[[31 159 0];[0 0 255];[159 0 197]]/255;

% Initial conditions for the 3 subfigures
pop0_sub=[[bpop0 vpop0 0];[bpop0 0 imm0];[bpop0 vpop0 imm0];[1e4 0 imm0]];

infect_red=@(t,y)infection_immune_bistable(t,y,para);
    rel_tol=1e-8;
    abs_tol=1e-10.*ones(3,1);
    options = odeset('RelTol',rel_tol,'AbsTol',abs_tol);

tic
%% Plotting each subfigure

% Setting axis for plot
xmin=0; xmax=TT;
ymin=1; ymax=1e10;

% Handles for curves
hand=zeros(4,3);

for subfig=1:4
% Initial conditions
pop0=pop0_sub(subfig,:);
% Solving the set of ODEs
    [T,Y] = ode45(infect_red,[0:tstep:TT],transpose(pop0),options);
    bpop(:)=Y(:,1);
    vpop(:)=Y(:,2);
    imm(:)=Y(:,3);
    run('plot_tseries')
end
toc

    % Saving the figures
    leg_size=15;
    figure(1);
    h_leg=legend([hand(1,2) hand(1,1)],{'Phage','Bacteria'},'Location','Southeast');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig3a.fig')
    print('fig3a.eps','-depsc')
    figure(2);
    h_leg=legend([hand(2,1) hand(2,3)],{'Bacteria','Immunity'},'Location','Southwest');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig3b.fig')
    print('fig3b.eps','-depsc')
    figure(3);
    h_leg=legend([hand(3,3) hand(3,2) hand(3,1)],{'Immunity','Phage','Bacteria'},'Location','Southeast');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig3c.fig')
    print('fig3c.eps','-depsc')