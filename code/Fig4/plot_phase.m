%% Plotting the time series and phase portraits of bacteria and virus populations
    lwidth=3; lbsize=25; tlbsize=20;
    set(0,'DefaultLineLinewidth',lwidth)
    set(0,'DefaultAxesLinewidth',lwidth)
    set(0,'DefaultAxesFontSize',tlbsize)
    bpop(bpop<para.thres)=0;
    vpop(vpop<para.thres)=0;
    
    % Time series
    figure(2*subfig-1);
    hand(iphi,ic,1)=semilogy(T,bpop,lstyle,'color',colormat(1,:));
    hold on
    hand(iphi,ic,2)=semilogy(T,vpop,lstyle,'color',colormat(2,:));
    %hold on
    %semilogy(T,imm,lstyle,'color',colormat(3,:));
    hold on
    plot([0 TT],[BIU BIU],':','color',[0 0 0])
    xlabel('Time (h)','FontSize',lbsize)
    ylabel('Density (ml^{-1})','FontSize',lbsize)
    axis(subfig_ax(2*subfig-1,:));
    xl=xlim; 
    set(gca,'Xtick',[xl(1):20:xl(2)]); 
    yl=ylim;
    logyl=log10(yl);
    set(gca,'Ytick',[10.^(logyl(1):log_int(subfig):logyl(2))]);
    
    % Phase portrait
    figure(2*subfig);
    t_len=phase_time(iphi,ic);
    hand(iphi,ic,3)=loglog(bpop(T<t_len),vpop(T<t_len),lstyle,'color',phase_color);
    hold on
    plot([BIU BIU],[1e-5 1e12],':','color',[0 0 0])
    hold on
    scatter(BP,PB,200,[1 165/255 0],'+','linewidth',lwidth)
    if (PBI>0)
    hold on
    scatter(BP,PBI,200,[0 0 0],'x','linewidth',lwidth)
    end
    hold on
    scatter(bpop0,vpop0,100,[0 0 0],'filled')
    xlabel('B (ml^{-1})','FontSize',lbsize)
    ylabel('P (ml^{-1})','FontSize',lbsize)
    axis(subfig_ax(2*subfig,:));
    xl=xlim;
    logxl=log10(xl);
    set(gca,'Xtick',[10.^(logxl(1):log_int(subfig):logxl(2))]);
    yl=ylim;
    logyl=log10(yl);
    set(gca,'Ytick',[10.^(logyl(1):log_int(subfig):logyl(2))]);