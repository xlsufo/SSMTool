function [B,A,Fnl,Fext] = build_model_none(om1,om2,om3,zeta1,zeta2,zeta3,f1,type)

n   = 3;
mass = eye(n,n);
damp = [2*zeta1*om1, 0, 0;
    0, 2*zeta2*om2, 0;
    0, 0, 2*zeta3*om3];
stiff = [om1^2 0 0;
    0 om2^2 0;
    0 0 om3^2];

switch type
    case 'cubic'
        B = [damp mass;mass zeros(n)]; A = [-stiff zeros(n); zeros(n) mass];
        B = [B zeros(2*n,1); zeros(1,2*n+1)];
        A = [A [0;0;-1;0;0;0]; [0 0 1 0 0 0 0]];
        Fnl  = @(input) fnon(input,om1,om2,om3,'cubic');
        Fext = [f1;0;0;0;0;0;0];
        
    case 'sphere'
        B = [damp mass;mass zeros(n)]; A = [-stiff zeros(n); zeros(n) mass];
        B = [B zeros(2*n,1); zeros(1,2*n+1)];
        A = [A [0;0;2;0;0;0]; [0 0 -2 0 0 0 0]];
        Fnl  = @(input) fnon(input,om1,om2,om3,'sphere');
        Fext = [f1;0;0;0;0;0;0];
        
    case 'none'
        B = [damp mass;mass zeros(n)]; A = [-stiff zeros(n); zeros(n) mass];
        Fnl  = @(input) fnon(input,om1,om2,om3,'none');
        Fext = [f1;0;0;0;0;0];      
        
    otherwise
        error('type should be cubic/circle/none');
end

end


function y = fnon(z,om1,om2,om3,type)

switch type
    case 'cubic'
        x1 = z(1);
        x2 = z(2);
        x3 = z(3);
        mu = z(7);
        dgnl = [-3*x1^2; -3*x2^2; 0];
        y = [-fnl_spring(x1,x2,x3,om1,om2,om3)-dgnl*mu
            zeros(3,1)
            -x1^3-x2^3];

    case 'sphere'
        x1 = z(1);
        x2 = z(2);
        x3 = z(3);
        mu = z(7);
        dgnl = 2*[x1; x2; x3];
        y = [-fnl_spring(x1,x2,x3,om1,om2,om3)-dgnl*mu
            zeros(3,1)
            x1^2+x2^2+x3^2];

    case 'none'
        x1 = z(1);
        x2 = z(2);
        x3 = z(3);
        y = [-fnl_spring(x1,x2,x3,om1,om2,om3)
            zeros(3,1)];

end
end


function y = fnl_spring(x1,x2,x3,om1,om2,om3)
% nonlinearity due to spring

om1_2 = om1^2;
om2_2 = om2^2;
om3_2 = om3^2;
tmp1  = (om1_2+om2_2+om3_2)/2;
tmp2  = x1^2+x2^2+x3^2;
y     = zeros(3,1);
y(1)  = 0.5*om1_2*(3*x1^2+x2^2+x3^2)+om2_2*x1*x2+om3_2*x1*x3+tmp1*x1*tmp2;
y(2)  = 0.5*om2_2*(3*x2^2+x1^2+x3^2)+om1_2*x1*x2+om3_2*x2*x3+tmp1*x2*tmp2;
y(3)  = 0.5*om3_2*(3*x3^2+x1^2+x2^2)+om1_2*x1*x3+om2_2*x2*x3+tmp1*x3*tmp2;

end