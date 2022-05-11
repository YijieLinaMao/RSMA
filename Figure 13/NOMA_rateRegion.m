function WSR = NOMA_rateRegion(weights,H,SNRdB,tolerance)
u1=weights(1);
u2=weights(2);
u3=weights(3);

h1=H(:,:,1);
h2=H(:,:,2);
h3=H(:,:,3);

%SC-SIC  6 different decoding order
% decoding order 1: 1-->2-->3
WSR_order(1) = NOMA_rateRegion_order1(weights,H,SNRdB,tolerance);
% [cap_tot2 Capacity2   p_12 p_22 ] = Baseline_rateRegion_order2(weights,H,SNRdB,tolerance);

% decoding order 2: 2-->1-->3
weights_order2(1)=u2;
weights_order2(2)=u1;
weights_order2(3)=u3;
H_order2(:,:,1)=h2;
H_order2(:,:,2)=h1;
H_order2(:,:,3)=h3;
WSR_order(2) = NOMA_rateRegion_order1(weights_order2,H_order2,SNRdB,tolerance);

% decoding order 3: 1-->3-->2
weights_order3(1)=u1;
weights_order3(2)=u3;
weights_order3(3)=u2;
H_order3(:,:,1)=h1;
H_order3(:,:,2)=h3;
H_order3(:,:,3)=h2;
WSR_order(3) = NOMA_rateRegion_order1(weights_order3,H_order3,SNRdB,tolerance);

% decoding order 4: 3-->1-->2 
weights_order4(1)=u3;
weights_order4(2)=u1;
weights_order4(3)=u2;
H_order4(:,:,1)=h3;
H_order4(:,:,2)=h1;
H_order4(:,:,3)=h2;
WSR_order(4) = NOMA_rateRegion_order1(weights_order4,H_order4,SNRdB,tolerance);

% decoding order 5: 2-->3-->1 
weights_order5(1)=u2;
weights_order5(2)=u3;
weights_order5(3)=u1;
H_order5(:,:,1)=h2;
H_order5(:,:,2)=h3;
H_order5(:,:,3)=h1;
WSR_order(5) = NOMA_rateRegion_order1(weights_order5,H_order5,SNRdB,tolerance);

% decoding order 6: 3-->2-->1 (user location: 3, 2, 1)
weights_order6(1)=u3;
weights_order6(2)=u2;
weights_order6(3)=u1;
H_order6(:,:,1)=h3;
H_order6(:,:,2)=h2;
H_order6(:,:,3)=h1;
WSR_order(6) = NOMA_rateRegion_order1(weights_order6,H_order6,SNRdB,tolerance);

WSR=max(WSR_order);
% WSR_order