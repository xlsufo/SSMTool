function [fnl,dfnl] = get_fint(K,model,Null, Nullf, ud, u0, disp_idx)

mphmesh(model)

fnl = @(input) -internalforce(model,input,Null, Nullf, ud, u0, disp_idx) - K*input; % 

dfnl = @(input) internalforceJacobian(model,input,Null,Nullf,ud, u0, disp_idx) - K ;

end

