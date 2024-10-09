function F_semi = set_F_semi(obj,F_semi)
F_semi = cell(numel(F_semi),1);
for j = 2:numel(F_semi)
    % ensure double output
    if ~isempty(F_semi{j})
        F_semi{j} = @(input) double(F_semi{j}(input));
    end
end
end