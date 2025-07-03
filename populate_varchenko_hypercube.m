function Vq = populate_varchenko_hypercube(d) %Creates q-Varchenko matrix for the hyperplane arrangement forming the d-cube
    syms q
    assume(q, 'integer')
    Vq = sym(zeros(3^d,3^d));
    for i = 1:3^d
        for j = i:3^d
            Vq(i,j) = q^count_sep(i,j,d);
            Vq(j,i) = Vq(i,j);
        end
    end
end

function sep = count_sep(n,m,d) %Counts the degree of separation between regions n and m on the hyerplane arrangement corresponding to the d-cube
    sep = 0;
    for i = 1:d
        a = floor(mod(n-1,3^i)/3^(i-1));
        b = floor(mod(m-1,3^i)/3^(i-1));
        sep = sep + abs(a-b);
    end 
end