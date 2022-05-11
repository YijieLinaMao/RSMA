
    function [cap_tot, p_1, p_2, p_3,  p_123, ro_1_123,ro_2_123,ro_3_123]=RS_oneLayer_update_P1(weights,H,SNR,g,W)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Then we try to use CVX to solve the problem         % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u1=weights(1);
    u2=weights(2);
    u3=weights(3);
    
    [A_U,NT,N_user] = size(H); 
    
    h1=H(:,:,1);
    h2=H(:,:,2);
    h3=H(:,:,3);
    
    g_1_123=g{1};
    g_1=g{2};
    
    g_2_123=g{3};
    g_2=g{4};
    
    g_3_123=g{5};
    g_3=g{6}; 
    
    W_1_123=W{1};
    W_1=W{2};
    
    W_2_123=W{3};
    W_2=W{4};
    
    W_3_123=W{5};
    W_3=W{6};
    
    
    
    cvx_begin quiet
    p_12=zeros(NT,A_U);
    p_13=zeros(NT,A_U);
    p_23=zeros(NT,A_U);
    ro_1_12=0;
    ro_2_12=0; 
    ro_1_13=0; 
    ro_3_13=0; 
    ro_2_23=0;
    ro_3_23=0;
    
    variable p_1(NT,A_U) complex
    variable p_2(NT,A_U) complex
    variable p_3(NT,A_U) complex
    variable p_123(NT,A_U) complex
    variable ro_1_123 
    variable ro_2_123 
    variable ro_3_123 

    expression constraints(1,7);
    
    
    
    %T_k calculation
    T_1=square_abs(h1*p_1)+square_abs(h1*p_2)+square_abs(h1*p_3)+1;
    T_2=square_abs(h2*p_1)+square_abs(h2*p_2)+square_abs(h2*p_3)+1;
    T_3=square_abs(h3*p_1)+square_abs(h3*p_2)+square_abs(h3*p_3)+1;
    

    T_1_123=square_abs(h1*p_123)+T_1;
    T_2_123=square_abs(h2*p_123)+T_2;
    T_3_123=square_abs(h3*p_123)+T_3;
    
    %MSE calculation
    E_1_123=abs(g_1_123)^2*T_1_123-2*real(g_1_123*h1*p_123)+1;
    E_2_123=abs(g_2_123)^2*T_2_123-2*real(g_2_123*h2*p_123)+1;
    E_3_123=abs(g_3_123)^2*T_3_123-2*real(g_3_123*h3*p_123)+1;
    
    E_1=abs(g_1)^2*T_1-2*real(g_1*h1*p_1)+1;
    E_2=abs(g_2)^2*T_2-2*real(g_2*h2*p_2)+1;
    E_3=abs(g_3)^2*T_3-2*real(g_3*h3*p_3)+1;
    
    %Rate-WMMSE relationship
    Xi_1_123=W_1_123*E_1_123-log2(W_1_123);

    Xi_1=W_1*E_1-log2(W_1);
    
    Xi_2_123=W_2_123*E_2_123-log2(W_2_123);

    Xi_2=W_2*E_2-log2(W_2);
    
    Xi_3_123=W_3_123*E_3_123-log2(W_3_123);

    Xi_3=W_3*E_3-log2(W_3);
    
    %common rate part:
    WMMSE_1=ro_1_123+Xi_1;
    WMMSE_2=ro_2_123+Xi_2;
    WMMSE_3=ro_3_123+Xi_3;
    
    object_func=u1*(WMMSE_1)+u2*(WMMSE_2)+u3*(WMMSE_3);

    
    
    %objective function
    minimize(object_func)
    
    
    
    %constraints    
    constraints(1)=Xi_1_123-1-ro_1_123-ro_2_123-ro_3_123;
    constraints(2)=Xi_2_123-1-ro_1_123-ro_2_123-ro_3_123;
    constraints(3)=Xi_3_123-1-ro_1_123-ro_2_123-ro_3_123;
    constraints(4)=ro_1_123;
    constraints(5)=ro_2_123 ;
    constraints(6)=ro_3_123;  
    constraints(7)=trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_3'*p_3)+trace(p_123'*p_123)-SNR;

    
    subject to
        constraints<=zeros(size(constraints))
        
        
    cvx_end

    cap_tot=object_func;
    
    
end