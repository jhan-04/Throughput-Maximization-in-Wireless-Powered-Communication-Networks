%%Fig4. f(z)=z.*log(z)-z+1
z=[0:0.01:5];
y=z.*log(z)-z+1;
plot(z,y)
hold on
line(xlim,[1 1])
xlabel('z')
ylabel('f(z)')