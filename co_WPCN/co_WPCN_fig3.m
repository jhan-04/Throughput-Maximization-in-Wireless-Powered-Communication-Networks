clear all
clc
D10=10;
D20=5;
D12=5;
P0=10^(20/10);
eff=0.67;%energy conversion efficiency ay receiver of Ui
band=10^6;%
p=0.75;%a fixed portion of harvested energy

sigma = -160 + 10*log10(10^6); % noise power,

for a=2:3% channel pathloss
    
    
    h10=10^(-3)*D10^(-a);
    h20=10^(-3)*D20^(-a);
    h12=10^(-3)*D12^(-a);
    %t=[tou0 tou1 tou21 tou22 t21 t22 ];

    
    rho1_10=h10^2*eff*p*P0/10^((sigma)/10);
    rho1_12=h10*h12*eff*p*P0/10^((sigma)/10);
    rho2=h20^2*eff*p*P0/10^((sigma)/10);
    
    
    n=200;
    m=1;

    for i=1:n
        tou1=i/n;
        for j=1:n-i
            tou21=j/n;
            for k=1:n-i-j
                tou22=k/n;
                
                
                tou0=1-tou1 -tou21- tou22;
                for t21=0:1/n:tou0
                    t22=tou0-t21;
                    
                    
                    
                    R1_10=tou1*log2(1+rho1_10*tou0/tou1)*(tou1 -tou21- tou22<1)*(t21+t22<=tou0);
                    R1_12=tou1*log2(1+rho1_12*tou0/tou1)*(tou1 -tou21- tou22<1)*(t21+t22<=tou0);
                    R1_20=tou21*log2(1+rho2*t21/tou21)*(tou1 -tou21- tou22<1)*(t21+t22<=tou0);
                    R1_co(m)=min(R1_10+R1_20,R1_12);
                    R2_co(m)=tou22*log2(1+rho2*t22/tou22)*(tou1 -tou21- tou22<1)*(t21+t22<=tou0);
                    m=m+1;
                end
                
            end
            
        end
    end
    max(R1_co)
    max(R2_co)
    r=1;
    r2_co=0;
    r1_co=0;
    R2_co=round(R2_co,2);
    R1_co=round(R1_co,2);
    for i=0:1/10^2:max(R2_co)
        if isempty(max(R1_co(find(R2_co==i))))
        else
            r1_co(r)=max(R1_co(find(R2_co==i)));
            r2_co(r)=i;
            r=r+1;
        end
    end
    hold on
    plot(r2_co,r1_co)
     
    
    %%%%%%%%%%%%%%5
    
    
    n=1000;
    gam1 = eff*p*h10^2*P0/(10^((sigma)/10));%%==rho
    gam2 = eff*p*h20^2*P0/(10^((sigma)/10));
    T1=1/n:1/n:1;
    T2=T1';
    T0=(1-T1-T2);
    R1=(T1.*log2(1+gam1*(T0./T1))).*(T1+T2<1)+0;
    R2=(T2.*log2(1+gam2*(T0./T2))).*(T1+T2<1)+0;
    
    r=1;
   r1=0;
   r2=0;
    R2=round(R2,2);
    R1=round(R1,2);
    for i=0:1/10^2:max(R2)
        if isempty(max(R1(find(R2==i))))
        else
            r1(r)=max(R1(find(R2==i)));
            r2(r)=i;
            r=r+1;
        end
    end
    
    plot(r2,r1,':')
    
    
    
end

hold off

