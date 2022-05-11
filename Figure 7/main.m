 
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
% Fig. 7 of the above paper will be reproduced by running the Matlab 
% script 'main.m'. By changing the variable 'bias' and the channel 
% realizations, you can reproduce Fig. 7--Fig. 10.


clc;clear all; clf
%accuracy of convergence 
tolerance = 10^-6;

%channel angle between user1 and user2
theta1=0; 
theta2=[1,2,3,4]*pi/9;

%SNR in dB
SNRdB=20; 

%channel bias
bias=1;

%user weights
weight=[-3 -1:0.05:1 3];
u2=10.^weight;
u1=ones(1,length(u2));




%channel realization (NT=4)
H_BC(:,:,1)=[1 exp(1i*theta1) exp(1i*2*theta1) exp(1i*3*theta1)];

%channel realization (NT=2)
%H_BC(:,:,1)=[1 exp(1i*theta1)]; 

 
 
for i_theta2=1:length(theta2)

   %channel realization (NT=4)
   H_BC(:,:,2)=[1 exp(1i*theta2(i_theta2)) exp(1i*2*theta2(i_theta2)) exp(1i*3*theta2(i_theta2))]*bias;
%   %channel realization (NT=2)
%   H_BC(:,:,2)=[1 exp(1i*theta2(i_theta2)) ]*bias;
   
   H_MAC(:,:,1)=H_BC(:,:,1)';
   H_MAC(:,:,2)=H_BC(:,:,2)';
   
   for i_weight=1:length(u1)

        weights=[u1(i_weight),u2(i_weight)]
        Capacity_DPC = DPC_rateRegion(weights,H_MAC,SNRdB,tolerance);
        
        [Capacity_RS_order1,Capacity_RS_order2, P_common1, P_common2] = RS_rateRegion(weights,H_BC,SNRdB,tolerance);
        
        [Capacity_NOMA_order1, Capacity_NOMA_order2 ]= NOMA_rateRegion(weights,H_BC,SNRdB,tolerance);
        
        Capacity_MULP = MULP_rateRegion(weights,H_BC,SNRdB,tolerance);
        

          capacity_DPC_UE1(i_theta2,i_weight)=Capacity_DPC(1);
          capacity_DPC_UE2(i_theta2,i_weight)=Capacity_DPC(2);           
          
          capacity_MULP_UE1(i_theta2,i_weight)=Capacity_MULP(1);
          capacity_MULP_UE2(i_theta2,i_weight)=Capacity_MULP(2);

          capacity_NOMA_order1_UE1(i_theta2,i_weight)=Capacity_NOMA_order1(1);
          capacity_NOMA_order1_UE2(i_theta2,i_weight)=Capacity_NOMA_order1(2);
          
          capacity_NOMA_order2_UE1(i_theta2,i_weight)=Capacity_NOMA_order2(1);
          capacity_NOMA_order2_UE2(i_theta2,i_weight)=Capacity_NOMA_order2(2);

          capacity_RS_order1_UE1(i_theta2,i_weight) = Capacity_RS_order1(1);
          capacity_RS_order1_UE2(i_theta2,i_weight) = Capacity_RS_order1(2);
        
          capacity_RS_order2_UE1(i_theta2,i_weight) = Capacity_RS_order2(1);
          capacity_RS_order2_UE2(i_theta2,i_weight) = Capacity_RS_order2(2);


    end %end user weights
    
end %end theta
 
figure(1)

subplot(2,2,1)
plot(reshape(capacity_DPC_UE1(1,:),[length(u1),1]),reshape(capacity_DPC_UE2(1,:),[length(u2),1]),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on,grid on


x=[capacity_RS_order1_UE1(1,:) capacity_RS_order2_UE1(1,:)];
y=[capacity_RS_order1_UE2(1,:) capacity_RS_order2_UE2(1,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind),y1(ind_ini(1):ind),'-+','LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on,grid on

x=[capacity_NOMA_order1_UE1(1,:) capacity_NOMA_order2_UE1(1,:)];
y=[capacity_NOMA_order1_UE2(1,:) capacity_NOMA_order2_UE2(1,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'--','LineWidth',2,'Color',[0.9290    0.6940    0.1250 ])
hold on,grid on

x=capacity_MULP_UE1(1,:);
y=capacity_MULP_UE2(1,:);
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])

title('(a) \theta=\pi/9')
xlabel('R_{1,tot} (bit/s/Hz)')
ylabel('R_{2,tot} (bit/s/Hz)')

 
subplot(2,2,2)
plot(reshape(capacity_DPC_UE1(2,:),[length(u1),1]),reshape(capacity_DPC_UE2(2,:),[length(u2),1]),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on,grid on

x=[capacity_RS_order1_UE1(2,:) capacity_RS_order2_UE1(2,:)];
y=[capacity_RS_order1_UE2(2,:) capacity_RS_order2_UE2(2,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-+','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on,grid on

x=[capacity_NOMA_order1_UE1(2,:) capacity_NOMA_order2_UE1(2,:)];
y=[capacity_NOMA_order1_UE2(2,:) capacity_NOMA_order2_UE2(2,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini:ind(1)),y1(ind_ini:ind(1)),'--','LineWidth',2.5,'Color',  [0.9290    0.6940    0.1250])
hold on,grid on

x=capacity_MULP_UE1(2,:);
y=capacity_MULP_UE2(2,:);
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini:ind(1)),y1(ind_ini:ind(1)),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])
title('(b) \theta=2\pi/9')
xlabel('R_{1,tot} (bit/s/Hz)')
ylabel('R_{2,tot} (bit/s/Hz)')
 
subplot(2,2,3)
plot(reshape(capacity_DPC_UE1(3,:),[length(u1),1]),reshape(capacity_DPC_UE2(3,:),[length(u2),1]),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on,grid on


x=[capacity_RS_order1_UE1(3,:) capacity_RS_order2_UE1(3,:)];
y=[capacity_RS_order1_UE2(3,:) capacity_RS_order2_UE2(3,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-+','LineWidth',2,'Color',[0.8500    0.3250    0.0980] )
hold on,grid on

x=[capacity_NOMA_order1_UE1(3,:) capacity_NOMA_order2_UE1(3,:)];
y=[capacity_NOMA_order1_UE2(3,:) capacity_NOMA_order2_UE2(3,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
ind=ind(find(ind>ind_ini));
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'--','LineWidth',2,'Color', [0.9290    0.6940    0.1250])
hold on,grid on

x=capacity_MULP_UE1(3,:);
y=capacity_MULP_UE2(3,:);
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])

legend('DPC','RS','SC-SIC','MU-LP')
title('(c) \theta=\pi/3')
xlabel('R_{1,tot} (bit/s/Hz)')
ylabel('R_{2,tot} (bit/s/Hz)')
 
subplot(2,2,4)
plot(reshape(capacity_DPC_UE1(4,:),[length(u1),1]),reshape(capacity_DPC_UE2(4,:),[length(u2),1]),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on,grid on

x=[capacity_RS_order1_UE1(4,:) capacity_RS_order2_UE1(4,:)];
y=[capacity_RS_order1_UE2(4,:) capacity_RS_order2_UE2(4,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-+','LineWidth',2,'Color',[0.8500    0.3250    0.0980] )
hold on,grid on

x=[capacity_NOMA_order1_UE1(4,:) capacity_NOMA_order2_UE1(4,:)];
y=[capacity_NOMA_order1_UE2(4,:) capacity_NOMA_order2_UE2(4,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'--','LineWidth',2,'Color', [0.9290    0.6940    0.1250])
hold on,grid on

x=capacity_MULP_UE1(4,:);
y=capacity_MULP_UE2(4,:);
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])

title('(d) \theta=4\pi/9')
xlabel('R_{1,tot} (bit/s/Hz)')
ylabel('R_{2,tot} (bit/s/Hz)')



