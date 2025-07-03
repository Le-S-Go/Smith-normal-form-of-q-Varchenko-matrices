%Creates a q-Varchenko matrix for the hyperplane arrangement corresponding
%to the d-dimensional simplex
function simplex_varchenko = populate_varchenko_simplex(d)
    syms q
    assume(q,'integer')
    simplex_varchenko = sym(zeros(2^(d+1)-1)); 
    for i = 1:2^(d+1)-1
        for j = i:2^(d+1)-1
            simplex_varchenko(i,j) = q^count_sep_simp(i,j,d);
            simplex_varchenko(j,i) = simplex_varchenko(i,j); 
        end
    end
end

function region_label = label_region(n)
    base_cases = [0,2,3,4,6,7,5]; 
    if n<=7 %base cases for 1 or 2 dimensions
        region_label = base_cases(n);
        return; 
    end
    dimension = floor(log2(n)); %determines the dimension
    if n == (2^(dimension+1)) - 1 %determines if n is the extra region that gets created each time we move into a new dimension
        region_label = 2^dimension + 1; %corresponds with 1 a bunch of 0s and then another 1
        return;
    end
    region_label = label_region(n+1-(2^dimension)) + 2^dimension; %recursively defines the region label by adding 1 to the front of the binary index from one lower dimension
end

function sep = count_sep_simp(n,m,dim) %calculates the bitwise difference between two binary strings
        sep = 0;
        for i = 1:(dim+1)
            a = floor(mod(label_region(n),2^i)/2^(i-1)); %calculates 
            b = floor(mod(label_region(m),2^i)/2^(i-1));
            sep = sep + abs(a-b);
        end 
end

