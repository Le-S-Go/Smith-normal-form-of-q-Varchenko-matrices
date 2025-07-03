function SNF = smith_normalize_hypercube(d,varargin)
    syms q
    assume(q, 'integer')
    transition = generate_transition_matrix(d);
    permutation = generate_final_permutation(d); 
    %permutation = eye(3^d); %For getting non-permuted matrix
    varchenko = populate_varchenko_hypercube(d);
    SNF = simplify(permutation * transition * varchenko * transition' * permutation');
    if nargin >1
        disp("The q-Varchenko corresponding to the hyperplane arrangement forming the " + d + "-dimensional hypercube is:")
        disp(varchenko)
        disp("The transition matrix needed to diagonalize this matrix is:")
        disp(transition)
        disp("The permutation matrix needed to put this matrix into Smith normal form is:")
        disp(permutation)
        disp("And the Smith normal form is:")
    end
end

%These transition matrices will put the q-Varchenko matrix into a diagonal
%form
function transition_matrix = generate_transition_matrix(d)
    syms q
    assume(q,'integer')
    if d == 0 %Base case
        transition_matrix = [1];
        return;
    end
    highest_level = [eye(3^(d-1)) -q*eye(3^(d-1)) zeros(3^(d-1),3^(d-1)) ; zeros(3^(d-1),3^(d-1)) eye(3^(d-1)) -q*eye(3^(d-1)) ; zeros(3^(d-1),3^(d-1)) zeros(3^(d-1),3^(d-1)) eye(3^(d-1))];
    transition_matrix = simplify([generate_transition_matrix(d-1) zeros(3^(d-1),3^(d-1)) zeros(3^(d-1),3^(d-1)) ; zeros(3^(d-1),3^(d-1)) generate_transition_matrix(d-1) zeros(3^(d-1),3^(d-1)) ; zeros(3^(d-1),3^(d-1)) zeros(3^(d-1),3^(d-1)) generate_transition_matrix(d-1)] * highest_level); %Recursively define the transition matrix  NOT TIME EFFICIENT
end

%This permutation matrix will put the diagonal form into SNF
function permutation_matrix = generate_final_permutation(n) 
    index_list = generate_indices(n);
    id = eye(3^n);
    permutation_matrix = zeros(3^n,3^n);
    for i = 1:3^n
        permutation_matrix(1:3^n,i) = id(1:3^n,index_list(i));
    end
end

%Support function for generating the final permutation matrix
%Counts the number of m-cubes present in an n-cube
function count = subspace_count(n,m) 
    count = (2^(n-m))*nchoosek(n,m);
end 

%Support function for generating the final permutation matrix
%Converts a base10 number into base3
function base3 = convert_to_base3(input)
    if input == 0
        base3 = 0;
        return;
    end
    base3='';
    while input > 0
        r = string(mod(input,3));
        base3 = r + string(base3);
        input = floor(input / 3);
    end
end

%Support function for generating the final permutation matrix
%Calls convert_to_base3 and then counts the number of twos
function count = count2s(input)
    input = char(convert_to_base3(input));
    count = sum(input == '2');
end 

%Support function for generating the final permutation matrix
function index_list = generate_indices(n)
    dict_counter = 0;
    d = dictionary; %This dictionary stores the number of twos found in a ternary representation as a key and the indices corresponding to those number of twos as the value
    i = n;
    while i >= 0 %Populate dictionary
        d(i) = {[dict_counter+1:dict_counter+subspace_count(n,i)]};
        dict_counter = dict_counter + subspace_count(n,i) ;
        i = i - 1;
    end
    index_list = []; %The indices of this list should be ordered where their values are
    index_counter = 1; 
    for i = 1:3^n
        cell = d(count2s(i-1)); %Cell is now a cell containing a list of indices that i will need to go into
        index_list(i) = cell{1,1}(find(cell{1,1} ~= 0, 1)); %Assigns the i-th element of the list to be the first nonzero element in the cell list
        cell{1,1}(find(cell{1,1} ~= 0, 1)) = 0; %Sets the first nonzero element in the cell list to 0 so that the index will not be reassigned
        d(count2s(i-1)) = cell; %Updates the dictionary with the updated cell
    end
end



