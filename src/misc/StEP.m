function [v_out] = StEP(F,vectors, m, nl_order)

[multis,multi_idcs,multi_multiplicities] = unique(m.','rows');
num_vecs = size(multis,1);

% order of nonlinear function
switch nl_order
    
    case 1
        
        v1 = vectors{1};
        v_out = 1/2 * ( F(v1) - F(-v1));
        
    case 2
        
        v1 = vectors{1};
        v2 = vectors{2};
        % Process for getting quadratic coefficients from nonlinearity
        
        % Evaluate complex input via real vector inputs to nonlinearity
        F2 = @(v) (1i*F(real(v)+imag(v)) +(1-1i)*F(real(v)) - (1+1i)*F(imag(v)));
        
        if num_vecs == 1            
            % Same vector is input twice
            v_out = 1/2 * ( F2(v1) + F2(-v1));
            
        else
            % Two different vectors are input            
            v_out = 1/4 * ( F2(v1+v2) + F2(-v1-v2) - F2(v1-v2) - F2(-v1+v2) );
            
        end
    case 3
        
        % Evaluate complex input via real vector inputs to nonlinearity
        F3 = @(v)  2*F(real(v)) + (-1+1i)/2*F(real(v)+imag(v)) - (1+1i)/2*F(real(v)-imag(v)) -2i*F(imag(v));
        
        H = @(v)  1/2 * (F3(v)-F3(-v)); % Cubic function
                             
        if num_vecs == 1
            v1 = vectors{1};

            % All three input vectors are the same
            v_out = H(v1);
            
        elseif num_vecs == 2
            % Two of the input vectors are the same
            
            % check which vector is input twice
            if sum(1==multi_multiplicities ) == 2
                v1 = vectors{multi_idcs(1)};
                v2 = vectors{multi_idcs(2)};
            else
                v1 = vectors{multi_idcs(2)};
                v2 = vectors{multi_idcs(1)};                
            end
            
            v_out = 1/2 * (H(v1+v2) - H(v1-v2)) - H(v2);
            
            
            
        else
            v1 = vectors{1};
            v2 = vectors{2};
            v3 = vectors{3};
            % All input vectors are the same            
            v_out = H(v1+v2+v3) - H(v1+v2) - H(v1+v3) - H(v2+v3) + H(v1) + H(v2) + H(v3);
            
        end

end