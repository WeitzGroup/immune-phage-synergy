%% Plotting the heat map of bacteria/virus at steady state
    lwidth=3; lbsize=25; tlbsize=20;
    set(0,'DefaultLineLinewidth',lwidth)
    set(0,'DefaultAxesLinewidth',lwidth)
    set(0,'DefaultAxesFontSize',tlbsize)
    
    % Set tick labels and markers on figure
    %omg_tick=[0.5 1 1.5 2];
    omg_tick=omg_min:omg_max;
    omg_tick_lb={'10^{-1}','10^{0}','10^{1}','10^{2}'};
    %phi_tick=[0.05 0.5 1 1.5 2]*1d-10;
    phi_tick=phi_min:phi_max;
    phi_tick_lb={'10^{-12}','10^{-11}','10^{-10}','10^{-9}','10^{-8}'};
    BS_tick_lb={'10^{8}','10^{8.2}','10^{8.4}','10^{8.6}','10^{8.8}','10^9'};
    PS_tick_lb={'10^0','10^{2}','10^{4}','10^{6}','10^{8}','10^{10}','10^{12}'};
    markx=zeros(1,3); marky=[-9.5 -10.5 -11.5];
    msize=15;
    letters=strsplit(sprintf('%c\n','A':'C'));
    
% Read in data from file
%load('data_fig5')
    
% Calculate phi values to plot the different phi contour lines at threshold
phiIU_min=omg_min-log10(para.beta*BIU);
phiIS_min=omg_min-log10(para.beta*BIS);
phiIM_min=omg_min-log10(para.beta*BIM);
phiIU_max=omg_max-log10(para.beta*BIU);
phiIS_max=omg_max-log10(para.beta*BIS);
phiIM_max=omg_max-log10(para.beta*BIM);

%% Heat map of bacterial concentration
figure(1);
title('Bacteria','FontSize',lbsize);

% Take logarithm of data
loglastbpop=lastbpop; loglastbpop(loglastbpop<1)=0;
loglastbpop=log10(loglastbpop);

% Plot heat map
imagesc(log10(omg_range),log10(phi_range),loglastbpop); 
set(gca,'Ydir','normal','XTick',omg_tick,'YTick',phi_tick,...
    'XTickLabel',omg_tick_lb,'YTickLabel',phi_tick_lb);
xlabel('\omega (h^{-1})','fontsize',lbsize); ylabel('\phi (ml h^{-1})','fontsize',lbsize);
lowestValue = min(loglastbpop(loglastbpop(:)>0));
highestValue = max(loglastbpop(:));
map=parula(256);
colormap(map);
c_lim=[floor(lowestValue-2/256), ceil(highestValue)];
caxis(c_lim);
map(1,:)=[1,1,1];
colormap(map);
c=colorbar;
set(c,'Linewidth',lwidth,'Ytick',[c_lim(1):0.2:c_lim(2)],'Yticklabel',BS_tick_lb);
ylabel(c,'B_S (ml ^{-1})','fontsize',lbsize);

% Add threshold lines
xext = get(gca,'XLim');
yext = get(gca, 'YLim');
hold on
plot([xext(1), xext(2)], [yext(1)+phiIS_min-phi_min, yext(1)+phiIS_max-phi_min], 'k-');
hold on
plot([xext(1), xext(2)], [yext(1)+phiIU_min-phi_min, yext(1)+phiIU_max-phi_min], 'k--');
hold on
plot([xext(1), xext(2)], [yext(1)+phiIM_min-phi_min, yext(1)+phiIM_max-phi_min], 'k:');
text(markx, marky, letters(1:length(markx)),'Fontsize',msize,'HorizontalAlignment','center','FontWeight','bold')

% Save figure
saveas(gcf,'fig5a.fig')
print('fig5a.eps','-depsc')

%% Heat map of virus concentration

figure(2);
title('Phage','FontSize',lbsize);

% Take logarithm of data
loglastvpop=lastvpop; loglastvpop(loglastvpop<1)=0;
loglastvpop=log10(loglastvpop);

% Plot heat map
imagesc(log10(omg_range),log10(phi_range),loglastvpop); 
set(gca,'Ydir','normal','XTick',omg_tick,'YTick',phi_tick,...
    'XTickLabel',omg_tick_lb,'YTickLabel',phi_tick_lb);
xlabel('\omega (h^{-1})','fontsize',lbsize); ylabel('\phi (ml h^{-1})','fontsize',lbsize);
lowestValue = min(loglastvpop(loglastvpop(:)>0));
highestValue = max(loglastvpop(:));
map=parula(256);
colormap(map);
%c_lim=[floor(lowestValue-2/256), ceil(highestValue)];
c_lim=[0, ceil(highestValue)];
caxis(c_lim);
map(1,:)=[1,1,1];
colormap(map);
c=colorbar;
set(c,'Linewidth',lwidth,'Ytick',[c_lim(1):2:c_lim(2)],'Yticklabel',PS_tick_lb);
ylabel(c,'P_S (ml ^{-1})','fontsize',lbsize);

% Add threshold lines
xext = get(gca,'XLim');
yext = get(gca, 'YLim');
hold on
plot([xext(1), xext(2)], [yext(1)+phiIS_min-phi_min, yext(1)+phiIS_max-phi_min], 'k-');
hold on
plot([xext(1), xext(2)], [yext(1)+phiIU_min-phi_min, yext(1)+phiIU_max-phi_min], 'k--');
hold on
plot([xext(1), xext(2)], [yext(1)+phiIM_min-phi_min, yext(1)+phiIM_max-phi_min], 'k:');
text(markx, marky, letters(1:length(markx)),'Fontsize',msize,'HorizontalAlignment','center','FontWeight','bold')

% Save figure
saveas(gcf,'fig5b.fig')
print('fig5b.eps','-depsc')