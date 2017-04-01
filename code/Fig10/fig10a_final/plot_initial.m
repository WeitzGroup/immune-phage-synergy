%% Plotting the heat map of bacteria/virus at steady state
    lwidth=3; lbsize=25; tlbsize=20;
    set(0,'DefaultLineLinewidth',lwidth)
    set(0,'DefaultAxesLinewidth',lwidth)
    set(0,'DefaultAxesFontSize',tlbsize)
    
    % Set tick labels and markers on figure
    %B0_tick=[0.5 1 1.5 2];
    B0_tick=B0_min:2:B0_max;
    B0_tick_lb={'10^4','10^6','10^8','10^{10}','10^{12}'};
    %P0_tick=[0.05 0.5 1 1.5 2]*1d-10;
    P0_tick=P0_min:2:P0_max;
    P0_tick_lb={'10^4','10^6','10^8','10^{10}','10^{12}'};
    BS_tick_lb={'10^0','10^2','10^4','10^{6}','10^{8}'};
    
% Read in data from file
%load('data_fig10a')

%% Heat map of bacterial concentration
figure(1);

% Take logarithm of data
loglastbpop=lastbpop; loglastbpop(loglastbpop<1)=0;
loglastbpop=log10(loglastbpop);

% Plot heat map
imagesc(log10(B0_range),log10(P0_range),loglastbpop); 
set(gca,'Ydir','normal','XTick',B0_tick,'YTick',P0_tick,...
'XTickLabel',B0_tick_lb,'YTickLabel',P0_tick_lb);
xlabel('B_0 (ml^{-1})','fontsize',lbsize); ylabel('P_0 (ml^{-1})','fontsize',lbsize);
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
set(c,'Linewidth',lwidth,'Ytick',[c_lim(1):2:c_lim(2)],'Yticklabel',BS_tick_lb);
ylabel(c,'B_S (ml ^{-1})','fontsize',lbsize);

hold on
plot([log10(BIU), log10(BIU)], [P0_min,P0_max], 'k:');
%text(log10(BIU),P0_max+0.5,'B_I^U','Fontsize',lbsize,'HorizontalAlignment','center');

% Save figure
saveas(gcf,'fig10a.fig')
print('fig10a.eps','-depsc')