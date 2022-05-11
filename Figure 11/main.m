 
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
% Fig. 11 of the above paper will be reproduced by running the Matlab 
% script 'main.m'. By changing the variable 'bias' and the channel 
% realizations, you can reproduce Fig. 11--Fig. 12.


clc;clear all; clf
%accuracy of convergence 
tolerance = 10^-6;

%number of transmit antenna
NT=4; 

%number of receive antenna at each user
A_U=1; 

%SNR in dB
SNRdB=20; 

%channel angle between user1 and user2
theta1=0; 
theta2=[1,2,3,4]*pi/9;

%Number of channel samples in imperfect CSIT
M=1000;

%channel bias
bias=1;

%user weights
weight=[-3 -1:0.05:1 3];
u2=10.^weight;
u1=ones(1,length(u2));




%Estimated channel at BS (NT=4)
H_BC_estimate(:,:,1)=[1 exp(1i*theta1) exp(1i*2*theta1) exp(1i*3*theta1)];

%Estimated channel at BS (NT=2)
%H_BC_estimate(:,:,1)=[1 exp(1i*theta1)]; 


%Channel error
SNR = 10^(SNRdB/10);
P_t=SNR; %total transmission power
P_e=P_t^(-0.6); %error variance

for i=1:M
    randn('seed',28+2*i);
    %channel error (NT=4)
    H_BC_error_1(:,:,i)=((randn(A_U,NT)+j*randn(A_U,NT))/sqrt(2))*sqrt(P_e);
    H_BC_error_2(:,:,i)=((randn(A_U,NT)+j*randn(A_U,NT))/sqrt(2))*sqrt(P_e)*bias;
    
end


 
 
for i_theta2=1:length(theta2)

   %Estimated channel at BS (NT=4)
   H_BC_estimate(:,:,2)=[1 exp(1i*theta2(i_theta2)) exp(1i*2*theta2(i_theta2)) exp(1i*3*theta2(i_theta2))]*bias;
   
   %Estimated channel at BS (NT=2)
%    H_BC_estimate(:,:,2)=[1 exp(1i*theta2(i_theta2))]*bias;
   
   

   
   for i_weight=1:length(u1)
        i_theta2
        weights=[u1(i_weight),u2(i_weight)]

        [Capacity_RS_order1, Capacity_RS_order2, P_common1, P_common2] = RS_rateRegion(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB,tolerance)  
        
        [Capacity_NOMA_order1, Capacity_NOMA_order2 ]= NOMA_rateRegion(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB,tolerance)       

        Capacity_MULP = MULP_rateRegion(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB,tolerance)
            
          
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
% subplot(2,2,1)

subplot(2,2,1)


x=[capacity_RS_order1_UE1(1,:) capacity_RS_order2_UE1(1,:)];
y=[capacity_RS_order1_UE2(1,:) capacity_RS_order2_UE2(1,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-+','LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on,grid on

x=[capacity_NOMA_order1_UE1(1,:) capacity_NOMA_order2_UE1(1,:)];
y=[capacity_NOMA_order1_UE2(1,:) capacity_NOMA_order2_UE2(1,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
ind=find(x1==0);
[~,ind_ini]=max(x1);
ind=ind(find(ind>ind_ini));
plot(x1(ind_ini(1):ind(1)),y1(ind_ini:ind(1)),'--','LineWidth',2.5,'Color',[0.9290    0.6940    0.1250 ])
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


title('(a) \theta=\pi/9, \gamma=0.3')
xlabel('R_1 (bit/s/Hz)')
ylabel('R_2 (bit/s/Hz)')

subplot(2,2,2)
x=[capacity_RS_order1_UE1(2,:) capacity_RS_order2_UE1(2,:)];
y=[capacity_RS_order1_UE2(2,:) capacity_RS_order2_UE2(2,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-+','LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on,grid on

x=[capacity_NOMA_order1_UE1(2,:) capacity_NOMA_order2_UE1(1,:)];
y=[capacity_NOMA_order1_UE2(2,:) capacity_NOMA_order2_UE2(1,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
ind=find(x1==0);
[~,ind_ini]=max(x1);
ind=ind(find(ind>ind_ini));
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'--','LineWidth',2,'Color',[0.9290    0.6940    0.1250 ])
hold on,grid on

x=capacity_MULP_UE1(2,:);
y=capacity_MULP_UE2(2,:);
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(x1);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])


title('(b) \theta=2\pi/9, \gamma=0.3')
xlabel('R_1 (bit/s/Hz)')
ylabel('R_2 (bit/s/Hz)')

subplot(2,2,3)
x=[capacity_RS_order1_UE1(3,:) capacity_RS_order2_UE1(3,:)];
y=[capacity_RS_order1_UE2(3,:) capacity_RS_order2_UE2(3,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-+','LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on,grid on

x=[capacity_NOMA_order1_UE1(3,:) capacity_NOMA_order2_UE1(3,:)];
y=[capacity_NOMA_order1_UE2(3,:) capacity_NOMA_order2_UE2(3,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
ind=find(x1==0);
[~,ind_ini]=max(x1);
ind=ind(find(ind>ind_ini));
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'--','LineWidth',2,'Color',[0.9290    0.6940    0.1250 ])
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

 legend('RS','SC--SIC','MU--LP')
title('(c) \theta=\pi/3, \gamma=0.3')
xlabel('R_1 (bit/s/Hz)')
ylabel('R_2 (bit/s/Hz)')
subplot(2,2,4)
x=[capacity_RS_order1_UE1(4,:) capacity_RS_order2_UE1(4,:)];
y=[capacity_RS_order1_UE2(4,:) capacity_RS_order2_UE2(4,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
xx=floor(x1.*10^(5))./(10^(5));
ind=find(xx==0);
[~,ind_ini]=max(xx);
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'-+','LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on,grid on

x=[capacity_NOMA_order1_UE1(4,:) capacity_NOMA_order2_UE1(4,:)];
y=[capacity_NOMA_order1_UE2(4,:) capacity_NOMA_order2_UE2(4,:)];
k=convhull(x,y);
x1 = x(k);
y1 = y(k);
ind=find(x1==0);
[~,ind_ini]=max(x1);
ind=ind(find(ind>ind_ini));
plot(x1(ind_ini(1):ind(1)),y1(ind_ini(1):ind(1)),'--','LineWidth',2,'Color',[0.9290    0.6940    0.1250 ])
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

 
title('(d) \theta=4\pi/9, \gamma=0.3')
xlabel('R_1 (bit/s/Hz)')
ylabel('R_2 (bit/s/Hz)')

