function WSR = RS_rateRegion(weights,H,SNRdB,tolerance)

u1=weights(1);
u2=weights(2);
u3=weights(3);

h1=H(:,:,1);
h2=H(:,:,2);
h3=H(:,:,3);
    
% There are three messages to be decoded at each user 
% Decoding order of messages at layer 2: 12-->13-->23
WSR_order(1) = RS_rateRegion_order1(weights,H,SNRdB,tolerance);

% Decoding order of messages at layer 2: 12-->23-->13 
weights_order2(1)=u2;
weights_order2(2)=u1;
weights_order2(3)=u3;
H_order2(:,:,1)=h2;
H_order2(:,:,2)=h1;
H_order2(:,:,3)=h3;
WSR_order(2) = RS_rateRegion_order1(weights_order2,H_order2,SNRdB,tolerance);

% Decoding order of messages at layer 2: 13-->12-->23 
weights_order3(1)=u1;
weights_order3(2)=u3;
weights_order3(3)=u2;
H_order3(:,:,1)=h1;
H_order3(:,:,2)=h3;
H_order3(:,:,3)=h2;
WSR_order(3) = RS_rateRegion_order1(weights_order3,H_order3,SNRdB,tolerance);

% Decoding order of messages at layer 2: 13-->23-->12
weights_order4(1)=u3;
weights_order4(2)=u1;
weights_order4(3)=u2;
H_order4(:,:,1)=h3;
H_order4(:,:,2)=h1;
H_order4(:,:,3)=h2;
WSR_order(4) = RS_rateRegion_order1(weights_order4,H_order4,SNRdB,tolerance);

% Decoding order of messages at layer 2: 23-->12-->13 
weights_order5(1)=u2;
weights_order5(2)=u3;
weights_order5(3)=u1;
H_order5(:,:,1)=h2;
H_order5(:,:,2)=h3;
H_order5(:,:,3)=h1;
WSR_order(5) = RS_rateRegion_order1(weights_order5,H_order5,SNRdB,tolerance);

% Decoding order of messages at layer 2: 23-->13-->12
weights_order6(1)=u3;
weights_order6(2)=u2;
weights_order6(3)=u1;
H_order6(:,:,1)=h3;
H_order6(:,:,2)=h2;
H_order6(:,:,3)=h1;
WSR_order(6) = RS_rateRegion_order1(weights_order6,H_order6,SNRdB,tolerance);

WSR=max(WSR_order);
% WSR_order




