function [Capacity1 Capacity2 ] = NOMA_rateRegion(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB,tolerance) 

Capacity1 = NOMA_rateRegion_order1(M,weights,H_BC_estimate,H_BC_error_1,H_BC_error_2,SNRdB,tolerance) ;

%switch the order of UE1 and UE2
H_BC_estimate1(:,:,1)=H_BC_estimate(:,:,2);
H_BC_estimate1(:,:,2)=H_BC_estimate(:,:,1);


weights1(1)=weights(2);
weights1(2)=weights(1);

Capacity22 = NOMA_rateRegion_order1(M,weights1,H_BC_estimate1,H_BC_error_2,H_BC_error_1,SNRdB,tolerance) ;

Capacity2(1)=Capacity22(2);
Capacity2(2)=Capacity22(1);
