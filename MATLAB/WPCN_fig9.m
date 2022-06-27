clear all
clc
m=0;
mt1=0;
mt2=0;

gam1=10^(2.2);
gam2=10;
n=1000;

T1=1/n:1/n:1;
T2=T1';
T0=(1-T1-T2);
R1=(T1.*log2(1+gam1*(T0./T1))).*(T1+T2<1)+0;
R2=(T2.*log2(1+gam2*(T0./T2))).*(T1+T2<1)+0;
R=R1+R2;

for i=1:n
%     m=max(R(:,i));
%     row=find(R(:,i)==m);
% if m==0
%     row=1;
% end
%     U1(i)=R1(row,i);
% U2(i)=R2(row,i);
   m=max(max(R(:,1:i)));

[row col]=find(R==m);
U1(i)=R1(row,col);
U2(i)=R2(row,col);

 

    
end

hold on
plot(U1,U2)


for i=1:n

  m=max(max(R(1:i,:)));

[row col]=find(R==m);
U1(i)=R1(row,col);
U2(i)=R2(row,col);
    
end
axis([0,4.5,0,4.5])
grid on
hold on
plot(U1,U2)



hold off
