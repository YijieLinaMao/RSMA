
% This is a code package related to the following conference paper:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y. Mao, B. Clerckx and V. O. K. Li, "Rate-splitting multiple access for 
% downlink communication systems: bridging, generalizing, and outperforming 
% SDMA and NOMA." EURASIP Journal on Wireless Communications and Networking 
% 2018.1 (2018): 133.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% The code is written by Yijie (Lina) Mao
%
% The code is implemented in Matlab environment with CVX toolbox 
% assisted. 
%
% Fig. 6(b) of the above paper will be reproduced by running the Matlab 
% script 'main.m'. By changing the variables 'bias' (channel gain 
% difference between the users), 'NT'(number of transmit antenna), you 
% can reproduce Fig. 5--Fig. 6.





clc;clear all; clf

%number of transmit antenna
NT=4; 

%channel gain difference 
bias=0.3;


%SNR in dB
SNRdB=20; 

%accuracy of convergence 
tolerance = 10^-6;


%number of random channel
N_channel=100;

%user weights
weight=[-3 -1:0.05:1 3];
u2=10.^weight;
u1=ones(1,length(u2));


capacity_DPC_UE1_average =zeros(length(u1),1);
capacity_DPC_UE2_average =zeros(length(u1),1);
capacity_MULP_UE1_average=zeros(length(u1),1);
capacity_MULP_UE2_average=zeros(length(u1),1);
capacity_NOMA_order1_UE1_average=zeros(length(u1),1);
capacity_NOMA_order1_UE2_average=zeros(length(u1),1);
capacity_NOMA_order2_UE1_average=zeros(length(u1),1);
capacity_NOMA_order2_UE2_average=zeros(length(u1),1);
capacity_RS_order1_UE1_average=zeros(length(u1),1);
capacity_RS_order1_UE2_average=zeros(length(u1),1);
capacity_RS_order2_UE1_average=zeros(length(u1),1);
capacity_RS_order2_UE2_average=zeros(length(u1),1);
   
for i_weight=1:length(u1)
    
    for i_channel=1:N_channel
        weights=[u1(i_weight),u2(i_weight)]; 
        [i_weight, i_channel]
        randn('seed',i_channel*4)
        H_BC(:,:,1)=1/sqrt(2)*(randn(1,NT)+1i*randn(1,NT));
        H_BC(:,:,2)=1/sqrt(2)*(randn(1,NT)+1i*randn(1,NT))*bias;

        H_MAC(:,:,1)=H_BC(:,:,1)';
        H_MAC(:,:,2)=H_BC(:,:,2)';
        
        Capacity_DPC = DPC_rateRegion(weights,H_MAC,SNRdB,tolerance);
        
        [Capacity_RS_order1, Capacity_RS_order2, P_common1, P_common2] = RS_rateRegion(weights,H_BC,SNRdB,tolerance);  
        
        [Capacity_NOMA_order1, Capacity_NOMA_order2 ]= NOMA_rateRegion(weights,H_BC,SNRdB,tolerance); 
        
        Capacity_MULP = MULP_rateRegion(weights,H_BC,SNRdB,tolerance);
        

          capacity_DPC_UE1(i_channel)=Capacity_DPC(1);
          capacity_DPC_UE2(i_channel)=Capacity_DPC(2);           
          
          capacity_MULP_UE1(i_channel)=Capacity_MULP(1);
          capacity_MULP_UE2(i_channel)=Capacity_MULP(2);

          capacity_NOMA_order1_UE1(i_channel)=Capacity_NOMA_order1(1);
          capacity_NOMA_order1_UE2(i_channel)=Capacity_NOMA_order1(2);
          
          capacity_NOMA_order2_UE1(i_channel)=Capacity_NOMA_order2(1);
          capacity_NOMA_order2_UE2(i_channel)=Capacity_NOMA_order2(2);

          capacity_RS_order1_UE1(i_channel) = Capacity_RS_order1(1);
          capacity_RS_order1_UE2(i_channel) = Capacity_RS_order1(2);
        
          capacity_RS_order2_UE1(i_channel) = Capacity_RS_order2(1);
          capacity_RS_order2_UE2(i_channel) = Capacity_RS_order2(2);
        
    end %end user weights
    
    capacity_DPC_UE1_average(i_weight)=mean(capacity_DPC_UE1);
    capacity_DPC_UE2_average(i_weight)=mean(capacity_DPC_UE2);
    
    capacity_MULP_UE1_average(i_weight)=mean(capacity_MULP_UE1);
    capacity_MULP_UE2_average(i_weight)=mean(capacity_MULP_UE2);
    
    capacity_NOMA_order1_UE1_average(i_weight)=mean(capacity_NOMA_order1_UE1);
    capacity_NOMA_order1_UE2_average(i_weight)=mean(capacity_NOMA_order1_UE2);
    
    capacity_NOMA_order2_UE1_average(i_weight)=mean(capacity_NOMA_order2_UE1);
    capacity_NOMA_order2_UE2_average(i_weight)=mean(capacity_NOMA_order2_UE2);
    
    capacity_RS_order1_UE1_average(i_weight)=mean(capacity_RS_order1_UE1);
    capacity_RS_order1_UE2_average(i_weight)=mean(capacity_RS_order1_UE2);
    
    capacity_RS_order2_UE1_average(i_weight)=mean(capacity_RS_order2_UE1);
    capacity_RS_order2_UE2_average(i_weight)=mean(capacity_RS_order2_UE2);
    
    clear capacity_DPC_UE1
    clear capacity_DPC_UE2
    clear capacity_MULP_UE1
    clear capacity_MULP_UE2
    clear capacity_NOMA_order1_UE1
    clear capacity_NOMA_order1_UE2
    clear capacity_NOMA_order2_UE1
    clear capacity_NOMA_order2_UE2
    clear capacity_RS_order1_UE1
    clear capacity_RS_order1_UE2
    clear capacity_RS_order2_UE1
    clear capacity_RS_order2_UE2
end %end theta
 

figure(1)
plot(reshape(capacity_DPC_UE1_average,[length(u1),1]),reshape(capacity_DPC_UE2_average,[length(u2),1]),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on,grid on

x=[capacity_RS_order1_UE1_average capacity_RS_order2_UE1_average];
y=[capacity_RS_order1_UE2_average capacity_RS_order2_UE2_average];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind),y1(ind_ini(1):ind),'-+','LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on,grid on

x=[capacity_NOMA_order1_UE1_average capacity_NOMA_order2_UE1_average];
y=[capacity_NOMA_order1_UE2_average capacity_NOMA_order2_UE2_average];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
ind=ind(find(ind>ind_ini));
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'--','LineWidth',2,'Color',[0.9290    0.6940    0.1250 ])
hold on,grid on

x=capacity_MULP_UE1_average;
y=capacity_MULP_UE2_average;
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])
legend('DPC','RS','SC-SIC','MU-LP')
xlabel('R_1 (bit/s/Hz)')
ylabel('R_2 (bit/s/Hz)')

 
