%% Plotting the heat map of bacteria/virus at steady state
    lwidth=3; lbsize=25; tlbsize=20;
    set(0,'DefaultLineLinewidth',lwidth)
    set(0,'DefaultAxesLinewidth',lwidth)
    set(0,'DefaultAxesFontSize',tlbsize)
    
    % Set tick labels and markers on figure
    %B0_tick=[0.5 1 1.5 2];
    B0_tick=B0_min:2:B0_max;
    %P0_tick=[0.05 0.5 1 1.5 2]*1d-10;
    P0_tick=P0_min:2:P0_max;
    
% Read in data from file
%load('data_fig10c')

%% Heat map of bacterial concentration
figure(1);

% Take logarithm of data
loglastbpop=lastbpop; loglastbpop(loglastbpop<1)=0;
loglastbpop=log10(loglastbpop);

% Plot heat map
imagesc(log10(B0_range),log10(P0_range),loglastbpop); 
set(gca,'Ydir','normal','XTick',B0_tick,'YTick',P0_tick);
xlabel('log(B_0) [log(ml^{-1})]','fontsize',lbsize); ylabel('log(P_0) [log(ml^{-1})]','fontsize',lbsize);
lowestValue = min(loglastbpop(loglastbpop(:)>0));
highestValue = max(loglastbpop(:));
map=parula(256);
colormap(map);
if (highestValue>=0)
    c_lim=[floor(lowestValue-2/256), ceil(highestValue)];
else
    c_lim=[0, log10(para.KC)];
end
caxis(c_lim);
map(1,:)=[1,1,1];
colormap(map);
c=colorbar;
set(c,'Linewidth',lwidth,'Ytick',[c_lim(1):0.2:c_lim(2)]);
ylabel(c,'log(B_S) [log(ml ^{-1})]','fontsize',lbsize);

hold on
plot([log10(BIU), log10(BIU)], [P0_min,P0_max], 'k:');

% Save figure
saveas(gcf,'fig10c.fig')
print('fig10c.eps','-depsc')

%% Heat map of virus concentration

%figure(2);

% Take logarithm of data
% loglastvpop=lastvpop; loglastvpop(loglastvpop<1)=0;
% loglastvpop=log10(loglastvpop);

% Plot heat map
% imagesc(log10(B0_range),log10(P0_range),loglastvpop); 
% set(gca,'Ydir','normal','XTick',B0_tick,'YTick',P0_tick);
% xlabel('log(B_0) [log(ml^{-1})]','fontsize',lbsize); ylabel('log(P_0) [log(ml^{-1})]','fontsize',lbsize);
% lowestValue = min(loglastvpop(loglastvpop(:)>0));
% highestValue = max(loglastvpop(:));
% map=parula(256);
% colormap(map);
% if (highestValue>=0)
%     caxis([floor(lowestValue-2/256), ceil(highestValue)]);
% else
%     caxis([0, 12]);
% end
% map(1,:)=[1,1,1];
%colormap(map);
% c=colorbar;
% ylabel(c,'log(P_S) [log(ml ^{-1})]','fontsize',lbsize);

% Save figure
% saveas(gcf,'fig10b.fig')
% print('fig10b.eps','-depsc')