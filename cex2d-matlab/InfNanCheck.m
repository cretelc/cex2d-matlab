function InfNan = InfNanCheck(line)
    for j = 1:length(line)
        if ismember('Infinity', line(j))
            InfNan = Inf;
        elseif ismember('NaN', line(j))
            InfNan = NaN;
        end
    end
end
