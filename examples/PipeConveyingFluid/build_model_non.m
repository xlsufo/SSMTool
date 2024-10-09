function [mass,damp,stiff,fnl,fext] = build_model_non(n,u,beta,Gamma,alpha,BC,varargin)

l = 1;
disp('Getting linear and nonlinearity coefficients')
fileName = [BC,'_tensors_',num2str(n),'_',num2str(1),'.mat'];
try 
    load(fileName,'a_ijkl','b_ijkl','c_ijkl','g_ijkl','a','b','delta','delta1','intphi');   
    disp('Loaded coefficients from storage');
catch
    disp('Calculating coefficients');
    [a,b,delta,delta1,a_ijkl,b_ijkl,c_ijkl,g_ijkl,intphi] = cal_parameters(n,l,BC);
    disp('Saving coefficients');
    save(fileName,'a_ijkl','b_ijkl','c_ijkl','g_ijkl','a','b','delta','delta1','intphi','-v7.3')
end

mass  = delta;
stiff = delta1+(u^2-Gamma)*b;
alpha_ijkl = u^2*c_ijkl+a_ijkl;
beta_ijkl  = 2*u*sqrt(beta)*b_ijkl;
gamma_ijkl = g_ijkl;

% visoelastic damping
damp = alpha*delta1+2*u*sqrt(beta)*a;
f3 = @(input) cubic(input,n,alpha_ijkl,beta_ijkl,gamma_ijkl);
        
fnl = f3;
fext = intphi;
end

function y = cubic(input,n,alpha_ijkl,beta_ijkl,gamma_ijkl)

x = input(1:n);
v = input(n+1:end);
y = zeros(n,1);
for j=1:n
    for k=1:n
        for l=1:n
            y = y+alpha_ijkl(:,j,k,l)*x(j)*x(k)*x(l)+...
                beta_ijkl(:,j,k,l)*x(j)*x(k)*v(l)+...
                gamma_ijkl(:,j,k,l)*x(j)*v(k)*v(l);
        end
    end
end

end


