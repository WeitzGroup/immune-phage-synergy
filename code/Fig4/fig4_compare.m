% Simulate the phage-bacteria-immune model and plot time series and phase portraits 
% of the dynamics of the model for two different adsorption rates

% Figure caption: 
%Comparison of time series and phase portraits with and without immune
%response for two different adsorption rates. 
%(a) Time series of bacteria (green curves) and phage densities (blue curves) 
%with adsorption rate $\phi =2.0\times 10^{-10}$ ml h$^{-1}$. Solid lines are 
%for the case of $I=K_I$ while dashed lines are for the case without immune response ($I=0$). 
%(b) shows the phase portrait for the time series in (a). Black solid line 
%is the phase trajectory when $I=K_I$ while orange dashed line is the trajectory 
%when $I=0$. The initial conditions ($B_0=B_I^S$, $P_0=10^5$ ml$^{-1}$) of the 
%trajectory is marked by a filled circle and the fixed point for the case of 
%$I=0$ is denoted by a plus sign.
%(c) and (d) are the same as (a) and (b) but with an adsorption rate $\phi =3.0\times 10^{-11}$ ml h$^{-1}$. 
%The coexistence fixed point $\{B^*=B_{PI},P^*=P_{BI},I^*=K_I\}$ is marked by a cross in (d).

close all;
%% Initialize parameters
TT=500; % Length of time series (h)
tstep=0.01; % Time step (h)
ntime=int64(1+TT/tstep); % Total number of time points
para=struct('r',[],'KC',[],'KD',[],'phi',[],'beta',[],'omg',[], ...
    'eps',[],'alpha',[],'KI',[],'KN',[],'thres',[]);
para.thres=0; % Threshold below which extinction is assumed (ml^-1)
% Set to be 0 here to show the entire phase trajectory

% Bacteria parameters
para.r=1; % Growth rate of bacteria (h^-1)
para.KC=1e9; % Carrying capacity of bacteria (ml^-1)
para.KD=2.2e6; % Bacteria conc. at which immune response is half as effective (ml^-1)

bpop=zeros(ntime,1);

% Virus parameters
para.omg=1; % Decay rate of phage (h^-1)
para.beta=100; % Burst size of phage
vpop0=1e8; % Initial phage conc. (ml^-1) when included

vpop=zeros(ntime,1);

% Immune response parameters
para.eps=8.2e-8; % Killing rate parameter of immune response (ml h^-1)
para.alpha=0.97; % Max growth rate of immune response (h^-1)
para.KI=2.4e7; % Max capacity of immune response (ml^-1)
para.KN=1e5; % Bacteria conc. when immune response growth rate is half its max (ml^-1)

imm=zeros(ntime,1);

% Calculate equilibrium densities B_I^U and B_I^S as defined in Eq. (4)
disc=sqrt((para.KC+para.KD)^2/4-para.KC*para.KD*para.eps*para.KI/para.r);
BIU=(para.KC-para.KD)/2-disc;
BIS=(para.KC-para.KD)/2+disc;
%BIM=para.KD*(sqrt(para.KC*para.eps*para.KI/(para.KD*para.r))-1);

% Initial immune response

tic

imm0_values=[para.KI 0]; % Initial immune response values to test (ml^-1)
%bpop0_values=[BIS BIS]; % Initial bacteria conc. values to test (ml^-1)
bpop0_values=[1e7 1e7];
phi_values=[1e-8 5e-11]; % Adsorption rate of phage values to test (ml h^-1)

% Colors and line styles used in plotting different curves in time series
colormat=[[31 159 0];[0 0 255]]/255;
lstyle_vec={'-', '--'};
% Colors used in plotting different curves in phase portrait
phase_color_mat=[[0 0 0]; [1 165/255 0]];
% Set axis for the different subfigures
subfig_ax=[[0 100 1e-4 1e12];
    [1e-8 1e12 1e-4 1e12]; %1
    [0 100 1e2 1e12];
    [1e4 1e10 1e2 1e12]]; %1e4
log_int=[4,2];
% Duration of phase trajectory to plot for the subfigures and initial
% conditions
phase_time=[[70 70]; [20 20]];

% Set tolerance of simulation
rel_tol=1e-10;
abs_tol=1e-10.*ones(3,1);

subfig=1;

% Handles for curves
len_phi=length(phi_values);
len_imm0=length(imm0_values);
hand=zeros(len_phi,len_imm0,3);

% Loop through different phi values
for iphi=1:len_phi
para.phi=phi_values(iphi);

% Compute equilibrium levels of bacteria and phage (B_P and P_B) without
% immune response
BP=para.omg/(para.beta*para.phi);
PB=(para.r/para.phi)*(1-BP/para.KC);
PBI=PB-(para.eps*para.KI/para.phi)/(1+BP/para.KD);

% Set up system of equations
infect_red=@(t,y)infection_immune_bistable(t,y,para);
options = odeset('RelTol',rel_tol,'AbsTol',abs_tol);

% Loop through initial conditions with and without immune response
for ic=1:len_imm0;
    bpop0=bpop0_values(ic);
    imm0=imm0_values(ic);
    pop0=[bpop0 vpop0 imm0];
    
    % Solving the set of ODEs
    [T,Y] = ode45(infect_red,[0:tstep:TT],transpose(pop0),options);
    bpop(:)=Y(:,1);
    vpop(:)=Y(:,2);
    imm(:)=Y(:,3);
    pcolor=colormat(ic,:);
    lstyle=lstyle_vec{ic};
    phase_color=phase_color_mat(ic,:);
    run('plot_phase')
end
subfig=subfig+1;
end
toc

    % Saving figures
    leg_size=13;
    figure(1);
    h_leg=legend([hand(1,2,2) hand(1,2,1) hand(1,1,2) hand(1,1,1)],...
    {'Phage (I=0)','Bacteria (I=0)','Phage (I=K_I)','Bacteria (I=K_I)'},'Location','Southeast');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig4a.fig')
    print('fig4a.eps','-depsc')
    figure(2);
    h_leg=legend([hand(1,2,3) hand(1,1,3)],{'I=0','I=K_I',},'Location','Northwest');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig4b.fig')
    print('fig4b.eps','-depsc')
    figure(3);
    h_leg=legend([hand(2,2,2) hand(2,2,1) hand(2,1,2) hand(2,1,1)],...
    {'Phage (I=0)','Bacteria (I=0)','Phage (I=K_I)','Bacteria (I=K_I)'},'Location','Southeast');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig4c.fig')
    print('fig4c.eps','-depsc')
    figure(4);
    h_leg=legend([hand(2,2,3) hand(2,1,3)],{'I=0','I=K_I',},'Location','Northwest');
    set(h_leg,'box','off','FontSize',leg_size);
    saveas(gcf,'fig4d.fig')
    print('fig4d.eps','-depsc')