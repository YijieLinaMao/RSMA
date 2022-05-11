function [cap_tot  p_1 p_2 ]=NOMA_update_P_order1(H_BC_estimate,M,weights,SNR,W_1_1,W_2_1,W_2_2,t_1_1,t_2_1,t_2_2,...
                                     phi_1_1,phi_2_1,phi_2_2,f_1_1,f_2_1,f_2_2,v_1_1,v_2_1,v_2_2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %       Then we try to use CVX to solve the problem        % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u1=weights(1);
    u2=weights(2);
   
    [A_U,NT,N_user] = size(H_BC_estimate) ;
    
    cvx_begin quiet
    
   
    variable p_2(NT,A_U) complex
    variable p_1(NT,A_U) complex

    expression constraints;

    r_1_1=quad_form(p_1,phi_1_1)+quad_form(p_2,phi_1_1)+t_1_1-2*real(f_1_1'*p_1)+W_1_1-v_1_1;
    r_2_1=quad_form(p_1,phi_2_1)+quad_form(p_2,phi_2_1)+t_2_1-2*real(f_2_1'*p_1)+W_2_1-v_2_1;
    r_2_2=quad_form(p_2,phi_2_2)+t_2_2-2*real(f_2_2'*p_2)+W_2_2-v_2_2;

    object_1=u1*max(r_1_1,r_2_1);
    object_2=u2*r_2_2;
        
   %objective function
    object_func=object_1+object_2;
    
   
    minimize(object_func)
    
    
    
    %constraints
    constraints=trace(p_2'*p_2)+trace(p_1'*p_1)-SNR;

    
    subject to
        constraints<=zeros(size(constraints))
        
        
    cvx_end

    cap_tot=object_func;
    
    
end