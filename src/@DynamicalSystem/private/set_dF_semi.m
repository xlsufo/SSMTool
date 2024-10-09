function dF_semi = set_dF_semi(obj,dF_semi)
dF_semi = cell(numel(dF_semi),1);
for j = 2:numel(dF_semi)
    % ensure double output
    if ~isempty(dF_semi{j})
        dF_semi{j} = @(invecs) double(dF_semi{j}(invecs));
    end
end
end