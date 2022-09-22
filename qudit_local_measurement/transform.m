function [coes, strs] = transform(d,mat)
    % transform a d*d matrix into len=d Pauli strings
    coes = [];
    strs = [];
    len = 0;
    for i = 1:d
        for j = 1:d
            if mat(i,j) ~= 0
                % transform string into
                % 0: I+Z 1:10=X-iY 2:01=X+iY 3: 11= I-Z
                coe = [mat(i,j)/2^d];
                str = [-1];
                for k = 1:d
                    if k ~= i && k ~= j
                        coev = [1;1];   % Corresponding coefficient
                        vec = [0;3];    % Corresponding string
                    end
                    if k == i && k == j
                        coev = [1;-1];
                        vec = [0;3];
                    end
                    if k == i && k ~= j
                        coev = [1;-1i];
                        vec = [1;2];
                    end
                    if k ~= i && k == j
                        coev = [1;1i];
                        vec = [1;2];
                    end
                    [coe,str] = strmul(coe,str,coev,vec);
                end
                [coes, strs] = mergestr(coes,strs,coe,str);
            end
        end
    end
    strs = strs(:,2:end);
end
