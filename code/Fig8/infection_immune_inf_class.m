% Model of bacteria-phage-immune system with explicit infected class to be integrated
    function dy=infection(t,y,para)
    dy=zeros(4,1);
    bpop=y(1);
    fpop=y(2);
    vpop=y(3);
    imm=y(4);
    
    % Uninfected bacteria population change
    bpopt=bpop+fpop;
    dy(1)=para.r*bpop*(1-bpopt/para.KC)- ...
    bpop*para.phi*vpop- ...
    para.eps*imm*bpop/(1+bpopt/para.KD); 

    % Infected bacteria population change
    dy(2)=para.phi*bpop*vpop-para.eta*fpop-para.eps*imm*fpop/(1+bpopt/para.KD);

    % Viruses population change
    dy(3)=para.beta*para.eta*fpop-para.phi*bpop*vpop-para.omg*vpop;

    % Immune response change
    dy(4)=para.alpha*imm*(1-imm/para.KI)* ...
    bpopt/(bpopt+para.KN);
