function [cap_tot,p_1,p_2,p_c ]=RS_update_P_order1(H_BC_estimate,M,weights,SNR,W_c_1,W_c_2,W_1,W_2,t_c_1,t_c_2,t_1,t_2,...
                         phi_c_1,phi_c_2,phi_1,phi_2,f_c_1,f_c_2,f_1,f_2,v_c_1,v_c_2,v_1,v_2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Then we try to use CVX to solve the problem         % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u1=weights(1);
    u2=weights(2);
    
    [A_U,NT,N_user] = size(H_BC_estimate) ;
    
    cvx_begin quiet
    
    variable p_1(NT,A_U) complex
    variable p_2(NT,A_U) complex
    variable p_c(NT,A_U) complex
 
    expression constraints;


    %objective function
    
        common_1=quad_form(p_c,phi_c_1)+quad_form(p_1,phi_c_1)+quad_form(p_2,phi_c_1)+t_c_1-2*real(f_c_1'*p_c)+W_c_1-v_c_1;
        common_2=quad_form(p_c,phi_c_2)+quad_form(p_1,phi_c_2)+quad_form(p_2,phi_c_2)+t_c_2-2*real(f_c_2'*p_c)+W_c_2-v_c_2;
       
        c=max(common_1,common_2);
  
        %Private rate part
        private_1=quad_form(p_1,phi_1)+quad_form(p_2,phi_1)+t_1-2*real(f_1'*p_1)+W_1-v_1;
        private_2=quad_form(p_1,phi_2)+quad_form(p_2,phi_2)+t_2-2*real(f_2'*p_2)+W_2-v_2;
        
        object_1=u1*c+u1*private_1;
        object_2=u2*private_2;
        
   
    object_func=object_1+object_2;
    minimize(object_func)
    
    
    
    %constraints 
    constraints=trace(p_1'*p_1)+trace(p_2'*p_2)+trace(p_c'*p_c)-SNR;

    subject to
        constraints<=zeros(size(constraints))
        
        
    cvx_end

    cap_tot=object_func;
    
    
end