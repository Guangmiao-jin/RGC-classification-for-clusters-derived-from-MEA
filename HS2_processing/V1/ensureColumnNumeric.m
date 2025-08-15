 function x = ensureColumnNumeric(x)
        if isempty(x)
            x = [];
        elseif ~isnumeric(x)
            x = [];
        else
            x = x(:);  % force column vector
        end
 end