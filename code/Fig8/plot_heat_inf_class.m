%% Plotting the heat map of bacteria/virus at steady state

    % Plot settings
    lwidth=3; lbsize=25; tlbsize=20;
    set(0,'DefaultLineLinewidth',lwidth)
    set(0,'DefaultAxesLinewidth',lwidth)
    set(0,'DefaultAxesFontSize',tlbsize)
    omg_tick=omg_min:omg_max;
    omg_tick_lb={'10^{-1}','10^{0}','10^{1}','10^{2}'};
    phi_tick=phi_min:phi_max;
    phi_tick_lb={'10^{-12}','10^{-11}','10^{-10}','10^{-9}','10^{-8}'};
    BS_tick_lb={'10^{8}','10^{8.2}','10^{8.4}','10^{8.6}','10^{8.8}','10^9'};
    FS_tick_lb={'10^{7}','10^{7.2}','10^{7.4}','10^{7.6}','10^{7.8}','10^8'};
    PS_tick_lb={'10^0','10^{2}','10^{4}','10^{6}','10^{8}','10^{10}','10^{12}'};

% Read in data from file
%load('data_fig8')

% Calculate phi values to plot shifted threshold lines
phiIU_min=omg_min-log10(para.beta*BIU);
phiIS_min=omg_min-log10(para.beta*BIS);
phiIM_min=omg_min-log10(para.beta*BIM);
phiIU_max=omg_max-log10(para.beta*BIU);
phiIS_max=omg_max-log10(para.beta*BIS);
phiIM_max=omg_max-log10(para.beta*BIM);

beta_shiftM=log10((para.eta+para.eps*para.KI/(1+BIM/para.KD))/para.eta);
beta_shiftS=log10((para.eta+para.eps*para.KI/(1+BIS/para.KD))/para.eta);
beta_shiftU=log10((para.eta+para.eps*para.KI/(1+BIU/para.KD))/para.eta);

% Heat map of uninfected bacteria concentration
figure(1);

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
plot([xext(1), xext(2)], [yext(1)+beta_shiftM+(phiIM_min-phi_min), yext(1)+beta_shiftM+(phiIM_max-phi_min)], 'k:');
hold on
plot([xext(1), xext(2)], [yext(1)+beta_shiftS+(phiIS_min-phi_min), yext(1)+beta_shiftS+(phiIS_max-phi_min)], 'k-');
hold on
plot([xext(1), xext(2)], [yext(1)+beta_shiftU+(phiIU_min-phi_min), yext(1)+beta_shiftU+(phiIU_max-phi_min)], 'k--');

title('Uninfected bacteria','FontSize',lbsize);

% Save figure
saveas(gcf,'fig8a.fig')
print('fig8a.eps','-depsc')

% Heat map of infected bacteria concentration
figure(2);
loglastfpop=lastfpop; loglastfpop(loglastfpop<1)=0;
loglastfpop=log10(loglastfpop);

imagesc(log10(omg_range),log10(phi_range),loglastfpop); 
set(gca,'Ydir','normal','XTick',omg_tick,'YTick',phi_tick,...
    'XTickLabel',omg_tick_lb,'YTickLabel',phi_tick_lb);
xlabel('\omega (h^{-1})','fontsize',lbsize); ylabel('\phi (ml h^{-1})','fontsize',lbsize);
lowestValue = min(loglastfpop(loglastfpop(:)>0));
highestValue = max(loglastfpop(:));
map=parula(256);
colormap(map);
c_lim=[floor(lowestValue-2/256), ceil(highestValue)];
caxis(c_lim);
map(1,:)=[1,1,1];
colormap(map);
c=colorbar;
set(c,'Linewidth',lwidth,'Ytick',[c_lim(1):0.2:c_lim(2)],'Yticklabel',FS_tick_lb);
ylabel(c,'F_S (ml ^{-1})','fontsize',lbsize);

xext = get(gca,'XLim');
yext = get(gca, 'YLim');
hold on
plot([xext(1), xext(2)], [yext(1)+beta_shiftM+(phiIM_min-phi_min), yext(1)+beta_shiftM+(phiIM_max-phi_min)], 'k:');
hold on
plot([xext(1), xext(2)], [yext(1)+beta_shiftS+(phiIS_min-phi_min), yext(1)+beta_shiftS+(phiIS_max-phi_min)], 'k-');
hold on
plot([xext(1), xext(2)], [yext(1)+beta_shiftU+(phiIU_min-phi_min), yext(1)+beta_shiftU+(phiIU_max-phi_min)], 'k--');

title('Infected bacteria','FontSize',lbsize);

saveas(gcf,'fig8b.fig')
print('fig8b.eps','-depsc')

% Heat map of virus concentration
figure(3);
loglastvpop=lastvpop; loglastvpop(loglastvpop<1)=0;
loglastvpop=log10(loglastvpop);

imagesc(log10(omg_range),log10(phi_range),loglastvpop); 
set(gca,'Ydir','normal','XTick',omg_tick,'YTick',phi_tick,...
    'XTickLabel',omg_tick_lb,'YTickLabel',phi_tick_lb);
xlabel('\omega (h^{-1})','fontsize',lbsize); ylabel('\phi (ml h^{-1})','fontsize',lbsize);
lowestValue = min(loglastvpop(loglastvpop(:)>0));
highestValue = max(loglastvpop(:));
map=parula(256);
colormap(map);
c_lim=[floor(lowestValue-2/256), ceil(highestValue)];
caxis(c_lim);
map(1,:)=[1,1,1];
colormap(map);
c=colorbar;
set(c,'Linewidth',lwidth,'Ytick',[c_lim(1):2:c_lim(2)],'Yticklabel',PS_tick_lb);
ylabel(c,'P_S (ml ^{-1})','fontsize',lbsize);

xext = get(gca,'XLim');
yext = get(gca, 'YLim');
hold on
plot([xext(1), xext(2)], [yext(1)+beta_shiftM+(phiIM_min-phi_min), yext(1)+beta_shiftM+(phiIM_max-phi_min)], 'k:');
hold on
plot([xext(1), xext(2)], [yext(1)+beta_shiftS+(phiIS_min-phi_min), yext(1)+beta_shiftS+(phiIS_max-phi_min)], 'k-');
hold on
plot([xext(1), xext(2)], [yext(1)+beta_shiftU+(phiIU_min-phi_min), yext(1)+beta_shiftU+(phiIU_max-phi_min)], 'k--');

title('Free phage','FontSize',lbsize);

saveas(gcf,'fig8c.fig')
print('fig8c.eps','-depsc')