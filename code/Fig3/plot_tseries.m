%% Plotting the time series of bacteria and virus populations, and immune response
    lwidth=3; lbsize=25; tlbsize=20;
    set(0,'DefaultLineLinewidth',lwidth)
    set(0,'DefaultAxesLinewidth',lwidth)
    set(0,'DefaultAxesFontSize',tlbsize)
    bpop(bpop<para.thres)=0;
    vpop(vpop<para.thres)=0;
    % Time series
    if (subfig~=4)
    figure(subfig);
    hand(subfig,1)=semilogy(T,bpop,'color',colormat(1,:));
    hold on
    hand(subfig,2)=semilogy(T,vpop,'color',colormat(2,:));
    hold on
    hand(subfig,3)=semilogy(T,imm,'color',colormat(3,:));
    xlabel('Time (h)','FontSize',lbsize)
    ylabel('Density (ml ^{-1})','FontSize',lbsize)
    else
    figure(2);
    axes('Position',[.55 .3 .3 .3])
    box on
    hand(subfig,1)=semilogy(T,bpop,'color',colormat(1,:));
    hold on
    hand(subfig,2)=semilogy(T,vpop,'color',colormat(2,:));
    hold on
    hand(subfig,3)=semilogy(T,imm,'color',colormat(3,:));
    set(gca,'FontSize',15)
    end
    axis([xmin xmax ymin ymax]);
    set(gca,'Xtick',[xmin:20:xmax]); 
    logymin=log10(ymin); logymax=log10(ymax);
    set(gca,'Ytick',[10.^(logymin:2:logymax)]);