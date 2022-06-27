function [R_t, t,lam] = t_opt(lam_new, gam)

% 
% lam=10;
% t_new=10;
% a=0;
n=0;
lam=lam_new;
% while 1




while 1
%%%%%(15)
syms a b m
vars = [a b m];
eqns = [(log(1+a)-(a/(1+a)))== (m*log(2))/lam(1), log(1+b)-b/(1+b) == m*log(2)/lam(2),((lam(1)*gam(1))/(1+a))+(lam(2)*gam(2))/(1+b)==m*log(2)];
S = vpasolve(eqns,vars);
z(1)=eval(S.a);%%%eval==double
z(2)=double(S.b);
mu=double(S.m);

t0=1./(1+sum(gam./z));
t12=(gam./z)./(1+sum(gam./z));
t=[t0 t12];
R_t(1)=(t12(1).*log2(1+gam(1)*(t0./t12(1))));
R_t(2)=(t12(2).*log2(1+gam(2)*(t0./t12(2))));

% 
    a =(sum(t<0)==0)&& (sum(lam<-0.00000001)==0)&&min((real(t)==t));
    if a
        t=[t0 t12];
        break;
        
    else
        
       lam=abs(lam_new*rand()*1.2)
    end
end

end