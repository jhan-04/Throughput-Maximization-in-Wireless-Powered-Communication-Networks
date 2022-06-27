function [R_t, t] = t_opt_K(lam, gam,k)

%%%%%(15)
% syms a b m
% vars = [a b m];

syms z [1 k] 
syms m
vars = [z m];

eqns = [(log(1+z)-(z./(1+z)))== (m*log(2))./lam,sum((lam.*gam)./(1+z))==m*log(2)];
S=vpasolve(eqns,vars);
a=struct2cell(S);%%struct-> cell

for i=1:k

    zi(i)=a{i};
    
    
end
mu=a{k+1};
    
    
    
t0=1./(1+sum(gam./zi));
t=(gam./zi)./(1+sum(gam./zi));
R_t=(t.*log2(1+gam.*(t0./t)));


end