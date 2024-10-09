function [fnl,dfnl] = get_fint(K,C,model,Null, Nullf, ud)
figure
mphmesh(model,'mesh1')

fnl = @(input) internalforce(model,input,Null, Nullf, ud, K, C); % 
dfnl = @(input) internalforceJacobian(model,input,Null,Nullf,ud) - K -C;
end

