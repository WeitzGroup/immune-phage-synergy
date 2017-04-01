%% Plotting the time series of bacteria and virus populations, and immune response
    lwidth=3; lbsize=25; tlbsize=20;
    set(0,'DefaultLineLinewidth',lwidth)
    set(0,'DefaultAxesLinewidth',lwidth)
    set(0,'DefaultAxesFontSize',tlbsize)
    markx=5; marky=3e11;
    msize=20;
    letters=strsplit(sprintf('%c\n','A':'C'));
    bpop(bpop<para.thres)=0;
    vpop(vpop<para.thres)=0;
    % Time series
    hand(subfig,1)=semilogy(T,bpop,'color',colormat(1,:));
    hold on
    hand(subfig,2)=semilogy(T,vpop,'color',colormat(2,:));
    %hold on
    %semilogy(T,imm,'color',colormat(3,:));
    if (PBI(subfig)>0)
    hold on
    hand(subfig,3)=plot([0 TT],[BP(subfig) BP(subfig)],':','color',colormat(1,:));
    hold on
    hand(subfig,4)=plot([0 TT],[PBI(subfig) PBI(subfig)],':','color',colormat(2,:));
    text(TT+1,BP(subfig),'B_P','Fontsize',lbsize,'HorizontalAlignment','left');
    text(TT+1,PBI(subfig),'P_{BI}','Fontsize',lbsize,'HorizontalAlignment','left');
    end
    text(markx, marky, strcat({'Case '},letters(subfig)),...
        'Fontsize',msize,'HorizontalAlignment','left')
    xlabel('Time (h)','FontSize',lbsize)
    ylabel('Density (ml ^{-1})','FontSize',lbsize)
    axis([xmin xmax ymin ymax]);
    set(gca,'Xtick',[xmin:20:xmax]); 
    logymin=log10(ymin); logymax=log10(ymax);
    set(gca,'Ytick',[10.^(logymin:2:logymax)]);