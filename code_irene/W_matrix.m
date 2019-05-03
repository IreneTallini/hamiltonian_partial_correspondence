function W = W_matrix(par, r)
n = [1 r/par.k]'; %[2,1]
p = [1 1]'; % [2,1]
nn = n ./ norm(n);
W = zeros(par.k, par.k); %[k,k]
for i = 1:par.k
    for j = 1:par.k
        a = [i j]' - p;
        W(i,j) = exp(- par.sig * norm([i j])) * ...
            norm(det([nn a]));
    end
end
end
