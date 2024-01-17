function [err] = err_on_ff(w,Vec)
err = sum(sum( (w.^2).*sum(abs(Vec).^2,1)))^(1/2);
end

