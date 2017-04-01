%% Plotting the functional forms of the growth and death rates of bacteria
% with max immune response and no phage
%Figure caption: 
%Per capita rates of bacteria growth $G(B)\equiv r(1-B/K_C)$ (solid line) 
%and bacteria death $D(B)\equiv\frac{\epsilon K_I}{1+B/K_D}$ (dashed line) 
%caused by the maximum immune response $K_I$ as functions of bacteria density. 
%The dotted vertical lines mark the positions of $B_I^U$ and $B_I^S$.

    close all;
    clear all;

    % Plot settings
    lwidth=3; lbsize=30; tlbsize=20;
    set(0,'DefaultLineLinewidth',lwidth)
    set(0,'DefaultAxesLinewidth',lwidth)
    set(0,'DefaultAxesFontSize',tlbsize)
    r=1; % Growth rate of bacteria (h^-1)
    KC=1e9; % Carrying capacity of bacteria (ml^-1)
    KD=2.2e6; % Bacteria conc. at which immune response is half as effective (ml^-1)
    eps=8.2e-8; % Killing rate parameter of immune response (ml h^-1)
    KI=2.4e7; % Max capacity of immune response (ml^-1)

    % Calculate equilibrium densities B_I^U and B_I^S as defined in Eq. (4)
    disc=sqrt((KC+KD)^2/4-KC*KD*eps*KI/r);
    BIU=(KC-KD)/2-disc;
    BIS=(KC-KD)/2+disc;
    
    % Rates at equilibrium densities
    rate_BIU=r*(1-BIU/KC);
    rate_BIS=r*(1-BIS/KC);
    
    % Computing data points to be plotted
    % Main figure
    axis_main=[0,1e7,0,4];
    xtick_main=axis_main(1):2e6:axis_main(2); 
    ytick_main=axis_main(3):1:axis_main(4);
    Bvec_main=axis_main(1):(axis_main(2)-axis_main(1))/100:axis_main(2);
    growth_main=r*(1-Bvec_main./KC);
    death_main=eps*KI./(1+Bvec_main./KD);
    % Inset figure
    axis_inset=[9.9e8,KC,0,0.01];
    xtick_inset=axis_inset(1):0.05e8:axis_inset(2); 
    ytick_inset=axis_inset(3):0.005:axis_inset(4);
    Bvec_inset=axis_inset(1):(axis_inset(2)-axis_inset(1))/100:axis_inset(2);
    growth_inset=r*(1-Bvec_inset./KC);
    death_inset=eps*KI./(1+Bvec_inset./KD);
    
    %Plotting the rates
    figure(1);
    plot(Bvec_main,growth_main,'-k')
    hold on
    plot(Bvec_main,death_main,'--k')
    hold on
    line([BIU BIU],[0 rate_BIU],'Color',[1 0 0],'LineStyle',':')
    axis(axis_main);
    set(gca,'Xtick',xtick_main,'Ytick',ytick_main);
    xlabel('B (ml^{-1})','FontSize',lbsize)
    ylabel('Per capita rates (h^{-1})','FontSize',lbsize)
    
    axes('Position',[.53 .58 .3 .3])
    box on
    plot(Bvec_inset,growth_inset,'-k')
    hold on
    plot(Bvec_inset,death_inset,'--k')
    hold on
    line([BIS BIS],[0 rate_BIS],'Color',[0 0 1],'LineStyle',':')
    axis(axis_inset);
    set(gca,'Xtick',xtick_inset,'Ytick',ytick_inset);
    set(gca,'FontSize',13);

    saveas(gcf,'fig2.fig')
    print('fig2.eps','-depsc')