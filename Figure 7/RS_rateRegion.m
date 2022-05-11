function [ Capacity1,Capacity2,P_common1,P_common2] = RS_rateRegion(weights,H,SNRdB,tolerance)

% Order 1: The common rate is allocated to user-1 only
[Capacity1,P_common1 ] = RS_rateRegion_order1(weights,H,SNRdB,tolerance);



% Order 2: The common rate is allocated to user-2 only
H1(:,:,1)=H(:,:,2);
H1(:,:,2)=H(:,:,1);
weights1(1)=weights(2);
weights1(2)=weights(1);
[Capacity22,P_common2] = RS_rateRegion_order1(weights1,H1,SNRdB,tolerance);
Capacity2(2)=Capacity22(1);
Capacity2(1)=Capacity22(2);
