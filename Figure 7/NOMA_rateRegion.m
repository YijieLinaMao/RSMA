function [Capacity1 Capacity2 ] = NOMA_rateRegion(weights,H,SNRdB,tolerance)

%decoding order 1: message of UE1-->message of UE2
Capacity1 = NOMA_rateRegion_order1(weights,H,SNRdB,tolerance);


%switch the order of UE1 and UE2
H1(:,:,1)=H(:,:,2);
H1(:,:,2)=H(:,:,1);

weights1(1)=weights(2);
weights1(2)=weights(1);


%decoding order 2: message of UE2-->message of UE1
Capacity22 = NOMA_rateRegion_order1(weights1,H1,SNRdB,tolerance);

Capacity2(1)=Capacity22(2);
Capacity2(2)=Capacity22(1);
