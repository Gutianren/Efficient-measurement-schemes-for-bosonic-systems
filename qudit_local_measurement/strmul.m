function [coes, strs] = strmul(coe, str, coe1, str1)
    % Multiply two coe-str table
    % Supposing the two parts are unique in themselves
    cnt = 0;
    for i = 1:length(coe)
        for j = 1:length(coe1)
            cnt = cnt + 1;
            coes(cnt) = coe(i) * coe1(j);
            if isempty(str)
                strs(cnt,:) = str1(j,:);
            else
                strs(cnt,:) = [str(i,:), str1(j,:)];
            end
        end
    end
    coes = coes';
end