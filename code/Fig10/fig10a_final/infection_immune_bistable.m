
% Model of bacteria-phage-immune system to be integrated
function dy=infection(t,y,para)
    dy=zeros(3,1);
    bpop=y(1);
    bpop(bpop<para.thres)=0;
    vpop=y(2);
    vpop(vpop<para.thres)=0;
    imm=y(3);
    
    % Bacteria population change
    dy(1)=para.r*bpop*(1-bpop/para.KC)- ...
    bpop*para.phi*vpop- ...
    para.eps*imm*bpop/(1+bpop/para.KD); 
    
    % Viruses population change
    dy(2)=para.beta*para.phi*vpop*bpop ...
    -para.omg*vpop;
    
    % Immune response change
    dy(3)=para.alpha*imm*(1-imm/para.KI)* ...
    bpop/(bpop+para.KN);
