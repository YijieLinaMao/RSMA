function [cap_tot, p_1, p_2, p_3, p_12, p_13, p_23,  p_123, ro_1_123,ro_2_123,ro_3_123,ro_1_12,ro_2_12,ro_1_13,ro_3_13,ro_2_23,ro_3_23]=RS_update_P(weights,H,SNR,g,W)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Then we try to use CVX to solve the problem        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u1=weights(1);
    u2=weights(2);
    u3=weights(3);
    
    [A_U,NT,N_user] = size(H); 
    
    h1=H(:,:,1);
    h2=H(:,:,2);
    h3=H(:,:,3);
    
    g_1_123=g{1};
    g_1_12=g{2};
    g_1_13=g{3};
    g_1=g{4};
    
    g_2_123=g{5};
    g_2_12=g{6};
    g_2_23=g{7};
    g_2=g{8};
    
    g_3_123=g{9};
    g_3_13=g{10};
    g_3_23=g{11};
    g_3=g{12}; 
    
    W_1_123=W{1};
    W_1_12=W{2};
    W_1_13=W{3};
    W_1=W{4};
    
    W_2_123=W{5};
    W_2_12=W{6};
    W_2_23=W{7};
    W_2=W{8};
    
    W_3_123=W{9};
    W_3_13=W{10};
    W_3_23=W{11};
    W_3=W{12};
    
    
    
    cvx_begin quiet
    
    variable p_1(NT,A_U) complex
    variable p_2(NT,A_U) complex
    variable p_3(NT,A_U) complex
    variable p_12(NT,A_U) complex
    variable p_13(NT,A_U) complex
    variable p_23(NT,A_U) complex
    variable p_123(NT,A_U) complex
    variable ro_1_123 
    variable ro_2_123 
    variable ro_3_123 
    variable ro_1_12
    variable ro_2_12 
    variable ro_1_13 
    variable ro_3_13 
    variable ro_2_23
    variable ro_3_23
    expression constraints(1,19);
    
    
    
    %T_k calculation
    T_1=square_abs(h1*p_23)+square_abs(h1*p_1)+square_abs(h1*p_2)+square_abs(h1*p_3)+1;
    T_2=square_abs(h2*p_13)+square_abs(h2*p_1)+square_abs(h2*p_2)+square_abs(h2*p_3)+1;
    T_3=square_abs(h3*p_12)+square_abs(h3*p_1)+square_abs(h3*p_2)+square_abs(h3*p_3)+1;
    
    T_1_13=square_abs(h1*p_13)+T_1;
    T_1_12=square_abs(h1*p_12)+T_1_13;  %s_12 is decoded before s_13
    T_1_123=square_abs(h1*p_123)+T_1_12;
    
    T_2_23=square_abs(h2*p_23)+T_2;
    T_2_12=square_abs(h2*p_12)+T_2_23;
    T_2_123=square_abs(h2*p_123)+T_2_12;
    
    T_3_23=square_abs(h3*p_23)+T_3;
    T_3_13=square_abs(h3*p_13)+T_3_23;
    T_3_123=square_abs(h3*p_123)+T_3_13;
    
    %MSE calculation
    E_1_123=abs(g_1_123)^2*T_1_123-2*real(g_1_123*h1*p_123)+1;
    E_1_12=abs(g_1_12)^2*T_1_12-2*real(g_1_12*h1*p_12)+1;
    E_1_13=abs(g_1_13)^2*T_1_13-2*real(g_1_13*h1*p_13)+1;
    E_1=abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;
    
    E_2_123=abs(g_2_123)^2*T_2_123-2*real(g_2_123*h2*p_123)+1;
    E_2_12=abs(g_2_12)^2*T_2_12-2*real(g_2_12*h2*p_12)+1;
    E_2_23=abs(g_2_23)^2*T_2_23-2*real(g_2_23*h2*p_23)+1;
    E_2=abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;
    
    E_3_123=abs(g_3_123)^2*T_3_123-2*real(g_3_123*h3*p_123)+1;
    E_3_13=abs(g_3_13)^2*T_3_13-2*real(g_3_13*h3*p_13)+1;
    E_3_23=abs(g_3_23)^2*T_3_23-2*real(g_3_23*h3*p_23)+1;
    E_3=abs(g_3)^2*T_3-2*real(g_3*h3*p_3)+1;
    
    %Rate-WMMSE relationship
    Xi_1_123=W_1_123*E_1_123-log2(W_1_123);
    Xi_1_12=W_1_12*E_1_12-log2(W_1_12);
    Xi_1_13=W_1_13*E_1_13-log2(W_1_13);
    Xi_1=W_1*E_1-log2(W_1);
    
    Xi_2_123=W_2_123*E_2_123-log2(W_2_123);
    Xi_2_12=W_2_12*E_2_12-log2(W_2_12);
    Xi_2_23=W_2_23*E_2_23-log2(W_2_23);
    Xi_2=W_2*E_2-log2(W_2);
    
    Xi_3_123=W_3_123*E_3_123-log2(W_3_123);
    Xi_3_13=W_3_13*E_3_13-log2(W_3_13);
    Xi_3_23=W_3_23*E_3_23-log2(W_3_23);
    Xi_3=W_3*E_3-log2(W_3);
    
    %common rate part:
    WMMSE_1=ro_1_123+ro_1_12+ro_1_13+Xi_1;
    WMMSE_2=ro_2_123+ro_2_12+ro_2_23+Xi_2;
    WMMSE_3=ro_3_123+ro_3_13+ro_3_23+Xi_3;
    
    object_func=u1*(WMMSE_1)+u2*(WMMSE_2)+u3*(WMMSE_3);

    
    
    %objective function
    minimize(object_func)
    
    
    
    %constraints    
    constraints(1)=Xi_1_123-1-ro_1_123-ro_2_123-ro_3_123;
    constraints(2)=Xi_2_123-1-ro_1_123-ro_2_123-ro_3_123;
    constraints(3)=Xi_3_123-1-ro_1_123-ro_2_123-ro_3_123;
    constraints(4)=Xi_1_12-1-ro_1_12-ro_2_12;
    constraints(5)=Xi_2_12-1-ro_1_12-ro_2_12;
    constraints(6)=Xi_1_13-1-ro_1_13-ro_3_13;
    constraints(7)=Xi_3_13-1-ro_1_13-ro_3_13;
    constraints(8)=Xi_2_23-1-ro_2_23-ro_3_23;
    constraints(9)=Xi_3_23-1-ro_2_23-ro_3_23;
    constraints(10)=ro_1_123;
    constraints(11)=ro_2_123 ;
    constraints(12)=ro_3_123; 
    constraints(13)=ro_1_12;
    constraints(14)=ro_2_12 ;
    constraints(15)=ro_1_13 ;
    constraints(16)=ro_3_13 ;
    constraints(17)=ro_2_23;
    constraints(18)=ro_3_23;   
    constraints(19)=trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_3'*p_3)+trace(p_12'*p_12)+trace(p_13'*p_13)+trace(p_23'*p_23)+trace(p_123'*p_123)-SNR;

    
    subject to
        constraints<=zeros(size(constraints))
        
        
    cvx_end
    cap_tot=object_func;
    
    
end