%function [R,t_new]=t_co_opt(w,rho)





g_min=1000;
k=0;
tou0=0.25;
tou1=0.25;
tou21=0.25;
tou22=0.25;
t21=0.125;
t22=0.125;
t_new=[tou0,tou1, tou21, tou22, t21,t22];

R1_10=tou1*log2(1+rho(1)*tou0/tou1);%%%%2)4.
R1_12=tou1*log2(1+rho(2)*tou0/tou1);
R1_20=tou21*log2(1+rho(3)*t21/tou21);
R1_co=min(R1_10+R1_20,R1_12);
R2_co=tou22*log2(1+rho(3)*t22/tou22);
R=[R1_co, R2_co];

R_min=1000;
t_min=t_new;



v=[1-(tou0+tou1+tou21+tou22), tou0-(t21+t22),R1_10+R1_20-R1_co,R1_12-R1_co];

n=0;
lam=[w(1)*0.6  w(1)*0.4 w(1)/2 w(1)/2];
g_lam_new=sum(w.*R)+sum(lam.*v);
g_lam=100;

while abs(g_lam_new - g_lam) > abs(0.01*g_lam)
    
    a=log(2)*(lam(1)-lam(2))*rho(1)*rho(2);
    b=log(2)*(lam(1)-lam(2))*(rho(1)+rho(2))-w(1)*rho(1)*rho(2);
    c=log(2)*(lam(1)-lam(2))-lam(3)*rho(1)-lam(4)*rho(2);

    syms z1 z21 z22
    
    eq1=lam(3)*(log(1+rho(1)*z1)-rho(1)*z1/(1+rho(1)*z1))+lam(4)*(log(1+rho(2)*z1)-rho(2)*z1/(1+rho(2)*z1))==lam(1)*log(2);
    eq2=log(1+z21)-z21/(1+z21)==lam(1)*log(2)/lam(3);
    eq3=log(1+z22)-z22/(1+z22)==lam(1)*log(2)/w(2);
    
    vars = [z1 z21 z22 ];
    eqns = [eq1 eq2 eq3];
    S = vpasolve(eqns,vars);
    z1=real(double(S.z1));
    z21=real(double(S.z21));
    z22=real(double(S.z22));
    
 
    t=100;
        
        while  sum(abs(t-t_new))>=0.001

            tou1=tou0/z1;
            tou21=rho(3)*t21/z21;
            tou22= rho(3)*t22/z22;
            tou0=(sqrt(b^2-4*a*c)-b)*tou1/(2*a) ;
            t21=max(0, lam(3)*tou21/(lam(2)*log(2))-tou21/rho(3));
            t22=max(0, w(2)*tou22/(lam(2)*log(2))-tou22/rho(3));
            

            tou_sum=(tou0+tou1+tou21+ tou22);
            tou0=tou0/tou_sum;
            tou1=tou1/tou_sum;
            tou21= tou21/tou_sum;
            tou22=tou22/tou_sum;
            t_sum=t21+t22;
            
            t21=t21*tou0/t_sum;
            t22=t22*tou0/t_sum;
            
            
            t=t_new;
            t_new=[tou0,tou1, tou21, tou22, t21,t22];
            
        end
t_new
  
    a =(sum(t_new<0)==0)&& (sum(lam<-0.0000001)==0)&&min((real(t_new)==t_new));%&&min((real(lam)==lam));
    if a
        R1_10=tou1*log2(1+rho(1)*tou0/tou1);%%%%2)4.
        R1_12=tou1*log2(1+rho(2)*tou0/tou1);
        R1_20=tou21*log2(1+rho(3)*t21/tou21);
        R1_co=min(R1_10+R1_20,R1_12);
        R2_co=tou22*log2(1+rho(3)*t22/tou22);
        R=[R1_co, R2_co];
        
        v=[1-(tou0+tou1+tou21+tou22), tou0-(t21+t22),R1_10+R1_20-R1_co,R1_12-R1_co];
        
%         
%         s=lam(3)+lam(4);
%          lam(3)=lam(3)*w(1)/s;
%         lam(4)=lam(4)*w(1)/s;
%         
        k=k+1;
        step=1/sqrt(k+1);
        lam=lam-step.*v;
        
        g_lam=g_lam_new;
        g_lam_new=sum(w.*R)+sum(lam.*v)
       
    else
        lam=[w(2)*rand()  w(2)*rand() w(2)*rand() w(2)*rand()]
    end
    
    
    
    
end
end
