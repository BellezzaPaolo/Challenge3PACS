close all
clear
clc

Data=load("plotData.m");

n_col=4;
dati=zeros(length(Data)/n_col,n_col);

pointer=1;
for ii=1:n_col:length(Data)
    dati(pointer,:)=Data(ii:ii+n_col-1);
    pointer=pointer+1;
end

h=1./(dati(1:n_col,2)-1);
figure(1)
loglog(h,h,h,h.^2)
hold on
loglog(h,dati(1:n_col,4),'LineWidth',5)
loglog(h,dati(end-n_col+1:end,4),'LineWidth',2)
legend('h','h^2','err sequential','err parallel')
xlabel('h')
ylabel('error in L2 norm')

figure(2)
plot(h,dati(1:n_col,3),'LineWidth',2)
hold on
plot(h,dati(end-n_col+1:end,3),'LineWidth',2)
legend('time seq','time par')
xlabel('h')
ylabel('time')