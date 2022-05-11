

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
% Fig. 13 of the above paper will be reproduced by running the Matlab 
% script 'main.m'. By changing the variable 'weight' and 'NT', you can 
% reproduce Fig. 13, Fig. 14, Fig. 16.




clc;clear all; clf
SNRdBs=[0:5:30 ]; 

%accuracy of convergence 
tolerance = 10^-6;


%channel strength difference among users
bias1=1;    %equal channel strength between user1 and user2
bias2=0.3;  %10 dB channel gain difference between user1 and user3

%channel angle between user1 and user2
theta1=0; 
theta2=[1,2,3,4]*pi/9;


%weight allocated to each user
weights=[0.2,0.3,0.5];

%channel realization
H_BC(:,:,1)=[1 exp(1i*theta1) exp(1i*2*theta1) exp(1i*3*theta1)];

for i_SNR=1:length(SNRdBs)
    SNRdB=SNRdBs(i_SNR);
    
    for i_theta=1:length(theta2)

        [i_SNR i_theta]

        H_BC(:,:,2)=[1 exp(1i*theta2(i_theta)) exp(1i*2*theta2(i_theta)) exp(1i*3*theta2(i_theta))]*bias1;
        H_BC(:,:,3)=[1 exp(1i*2*theta2(i_theta)) exp(1i*2*2*theta2(i_theta)) exp(1i*2*3*theta2(i_theta))]*bias2;

        H_MAC(:,:,1)=H_BC(:,:,1)';
        H_MAC(:,:,2)=H_BC(:,:,2)';
        H_MAC(:,:,3)=H_BC(:,:,3)';

        WSR_DPC_eachChannel(i_SNR,i_theta) = DPC_rateRegion(weights,H_MAC,SNRdB,tolerance);
        WSR_RS_eachChannel(i_SNR,i_theta)   = RS_rateRegion(weights,H_BC,SNRdB,tolerance);  
        WSR_NOMA_eachChannel(i_SNR,i_theta) = NOMA_rateRegion(weights,H_BC,SNRdB,tolerance);
        WSR_MULP_eachChannel(i_SNR,i_theta) = MULP_rateRegion(weights,H_BC,SNRdB,tolerance); 
        WSR_RS_oneLayer_eachChannel(i_SNR,i_theta)   = RS_oneLayer_rateRegion1(weights,H_BC,SNRdB,tolerance);  
    end %end i_theta

end %end SNR
 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,2,1)
plot(SNRdBs,WSR_DPC_eachChannel(:,1),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on;
grid on;
plot(SNRdBs,WSR_RS_eachChannel(:,1),'-+','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on;
plot(SNRdBs,WSR_NOMA_eachChannel(:,1),'--','LineWidth',2,'Color',  [0.9290    0.6940    0.1250])
hold on;
plot(SNRdBs,WSR_MULP_eachChannel(:,1),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])
hold on;
plot(SNRdBs,WSR_RS_oneLayer_eachChannel(:,1),'diamond','LineWidth',2,'Color', [1 0 1])
hold on;
xlabel('SNR (dB)')
ylabel('Weighted Sum Rate (bit/s/Hz)')
title('\theta_1=\pi/9, \theta_2=2\pi/9')

subplot(2,2,2)
plot(SNRdBs,WSR_DPC_eachChannel(:,2),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on;
grid on;
plot(SNRdBs,WSR_RS_eachChannel(:,2),'-+','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on;
plot(SNRdBs,WSR_NOMA_eachChannel(:,2),'--','LineWidth',2,'Color',  [0.9290    0.6940    0.1250])
hold on;
plot(SNRdBs,WSR_MULP_eachChannel(:,2),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])
hold on;
plot(SNRdBs,WSR_RS_oneLayer_eachChannel(:,2),'diamond','LineWidth',2,'Color', [1 0 1])
hold on;
xlabel('SNR (dB)')
ylabel('Weighted Sum Rate (bit/s/Hz)')
title('\theta_1=2\pi/9, \theta_2=4\pi/9')

subplot(2,2,3)
plot(SNRdBs,WSR_DPC_eachChannel(:,3),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on;
grid on;
plot(SNRdBs,WSR_RS_eachChannel(:,3),'-+','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on;
plot(SNRdBs,WSR_NOMA_eachChannel(:,3),'--','LineWidth',2,'Color',  [0.9290    0.6940    0.1250])
hold on;
plot(SNRdBs,WSR_MULP_eachChannel(:,3),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])
hold on
plot(SNRdBs,WSR_RS_oneLayer_eachChannel(:,3),'diamond','LineWidth',2,'Color', [1 0 1])
hold on;;
legend('DPC','RS','SC-SIC','MU-LP','1-layer RS')
xlabel('SNR (dB)')
ylabel('Weighted Sum Rate (bit/s/Hz)')
title('\theta_1=\pi/3, \theta_2=2\pi/3')

subplot(2,2,4)
plot(SNRdBs,WSR_DPC_eachChannel(:,4),'-o','LineWidth',2,'Color', [  0    0.4470    0.7410])
hold on;
grid on;
plot(SNRdBs,WSR_RS_eachChannel(:,4),'-+','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on;
plot(SNRdBs,WSR_NOMA_eachChannel(:,4),'--','LineWidth',2,'Color',  [0.9290    0.6940    0.1250])
hold on;
plot(SNRdBs,WSR_MULP_eachChannel(:,4),'-.','LineWidth',2,'Color', [0.4660    0.6740    0.1880])
hold on;
plot(SNRdBs,WSR_RS_oneLayer_eachChannel(:,4),'diamond','LineWidth',2,'Color', [1 0 1])
hold on;
xlabel('SNR (dB)')
ylabel('Weighted Sum Rate (bit/s/Hz)')
title('\theta_1=8\pi/9, \theta_2=16\pi/9')