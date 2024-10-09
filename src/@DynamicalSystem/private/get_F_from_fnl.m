function F = get_F_from_fnl(obj)
% GET_F_FROM_FNL
% get first order force from second order force
F = {};
F_flag = false;
if ~isempty(obj.fnl)
    for j = 1:numel(obj.fnl)
        if size(obj.fnl{j},1) == obj.N
            F_flag = true; % F was set under the "fnl" property
            continue
        end
    end
end

if F_flag
    % F is set in fnl
    for j = 1:numel(obj.fnl)
        if isempty(obj.fnl{j})
            sizej = obj.N*ones(1,j+2);
            F{j+1} = sptensor(sizej);
        else
            F{j+1} = obj.fnl{j};
        end
    end
else
    % F needs to be set from fnl
    for j = 1:numel(obj.fnl)
        sizej = obj.N*ones(1,j+2);
        if isempty(obj.fnl{j})
            F{j+1} = sptensor(sizej);
        else
            F{j+1} = sptensor(obj.fnl{j}.subs,-obj.fnl{j}.vals,sizej);
        end
    end
end

end