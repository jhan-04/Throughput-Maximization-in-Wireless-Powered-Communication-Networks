function [tou0 ,tou1 ,tou21, tou22, t21, t22 ,lam]=find_t_lam(rho,w)


syms tou0 tou1 tou21 tou22 t21 t22 lam1 lam2 lam3 lam4
    var=[tou0 ,tou1 ,tou21, tou22, t21, t22, lam1 ,lam2, lam3, lam4];
    eqn0=w(1)-lam3-lam4==0;
    eqn1=lam3*rho(1)/(1+rho(1)*tou0/tou1)+lam4*rho(2)/(1+rho(2)*tou0/tou1)==log(2)*(lam1-lam2);
    eqn2=lam3*(log(1+rho(1)*tou0/tou1+rho(1)*tou0/tou1/(1+rho(1)*tou0/tou1)))+lam4*(log(1+rho(2)*tou0/tou1)-(rho(2)*tou0/tou1)/(1+rho(2)*tou0/tou1))==log(2)*lam1;
   
    eqn3=lam3*(log(1+rho(2)*t21/tou21)-(rho(2)*t21/tou21)/(1+rho(2)*t21/tou21))==lam1*log(2);
    eqn4=w(2)*(log(1+rho(2)*t22/tou22)-(rho(2)*t22/tou22)/(1+rho(2)*t22/tou22))==lam1*log(2);
    eqn5=lam3* rho(2)/(1+rho(2)*t21/tou21)==lam2*log(2);
    eqn6=w(2)* rho(2)/(1+rho(2)*t22/tou22)==lam2*log(2);
    eqn7=tou0+tou1+tou21+tou22==1;
    eqn8=t21+t22==tou0;
    eqn9=tou1*log2(1+rho(1)*tou0/tou1)+tou21*log2(1+rho(3)*t21/tou21)==tou1*log2(1+rho(2)*tou0/tou1);
    eqns=[eqn0 eqn1 eqn2  eqn3 eqn4 eqn5 eqn6 eqn7 eqn8 eqn9];
    S=vpasolve(eqns,var);
%     tou0=real(double(S.tou0));
%     tou1=real(double(S.tou1));
%     tou21 =real(double(S.tou21));
%     tou22=real(double(S.tou22));
%     t21=real(double(S.t21));
%     t22=real(double(S.t22));

    tou0=(double(S.tou0));
    tou1=(double(S.tou1));
    tou21 =(double(S.tou21));
    tou22=(double(S.tou22));
    t21=(double(S.t21));
    t22=(double(S.t22));

    lam1=double(S.lam1);
    lam2=double(S.lam2);
    lam3=double(S.lam3);
    lam4=double(S.lam4);
    lam=[ lam1 ,lam2, lam3, lam4];
 
    
    tou_new=[tou0 tou1 tou21 tou22];
    t_new=[tou0 tou1 tou21 tou22 t21 t22];



end