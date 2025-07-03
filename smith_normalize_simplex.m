function SNF = smith_normalize_simplex(d,varargin)
    syms q
    assume(q,'integer')
    Vq = populate_varchenko_simplex(d);
    tm1 = first_transition_matrix(d);
    R = reindex(d);
    %R = eye(2^(d+1)-1); %test
    tm2 = second_transition_matrix(d);
    %tm2 = eye(2^(d+1)-1); %test
    P = generate_final_simplex_permutation(d);
    SNF = simplify(P*tm2*R*tm1 * Vq * tm1' * R' * tm2' * P'); 
    if nargin >1
        disp("The q-Varchenko corresponding to the hyperplane arrangement forming the " + d + "-dimensional simplex is:")
        disp(Vq)
        disp("The transition matrix needed to block diagonalize this matrix is:")
        disp(tm1)
        disp("The permutation matrix needed to reindex this block diagonalized matrix is:")
        disp(R)
        disp("The transition matrix needed to diagonalize the reindexed matrix is:")
        disp(tm2)
        disp("The permutation matrix needed to put this matrix into Smith normal form is:")
        disp(P)
        disp("And the Smith normal form is:")
    end
end

%These are the transition matrices used before the reindexing step
%After these matrices, the q-Varchenko matrix for the d simplex will be in
%block diagonal form
function tm1 = first_transition_matrix(d) %tm is called p
    syms q
    assume(q,'integer')
    tm1 = eye(2^(d+1)-1);
    for k = 1:d-1
        p = [eye(2^(d-k+1)-1) -q*eye(2^(d-k+1)-1) zeros(2^(d-k+1)-1,1+4*(2^(d-1)-2^(d-k))) ; zeros(2^(d-k+1)-1,2^(d-k+1)-1) eye(2^(d-k+1)-1) zeros(2^(d-k+1)-1, 1+4*(2^(d-1)-2^(d-k))) ; zeros(1,2*(2^(d-k+1)-1)) 1 zeros(1,4*(2^(d-1)-2^(d-k))) ; zeros(4*(2^(d-1)-2^(d-k)),2*(2^(d-k+1)-1)+1) eye(4*(2^(d-1)-2^(d-k)))];
        tm1 = p * tm1;
    end
end

%Reindexing step
function R = reindex(d)
    R = zeros(2^(d+1)-1,2^(d+1)-1);
    R(1:3,1:3) = eye(3);
    counter = 3;
    for k = 2:d
        R(counter+1:counter+2^k,counter+1:counter+2^k) = one_layer_reindexing(k);
        counter = counter + 2^k;
    end
end

%Recursive support function for reindexing step 
function R_k = one_layer_reindexing(n)
    if n == 2
        R_k = [1 0 0 0 ; 0 0 0 1 ; 0 1 0 0 ; 0 0 1 0 ];
        return;
    end
    R_k = zeros(2^n,2^n);
    R_k(1,1) = 1;
    R_k(2,2^n) = 1;
    R_k(3,2) = 1;
    R_k(4,3) = 1;
    %R_k(n+2:2^n,n+1:2^n-1) = one_layer_reindexing(n-1);
    counter = 3;
    for i=2:n-1
        R_k(2+counter:1+counter+2^i,1+counter:counter+2^i) = one_layer_reindexing(i);
        counter = counter + 2^i;
    end
end

%These are the transition matrices used after the reindexing step
%They will put the block diagonal matrix into a fully diagonal form
function tm2 = second_transition_matrix(d)
    syms q
    assume(q,'integer')
    tm2 = sym(eye(2^(d+1)-1));
    for k=1:d-1
        S = sym(zeros(2^(d+1)-1,2^(d+1)-1));
        S(1:2^(k+1)-1,1:2^(k+1)-1) = eye(2^(k+1)-1);
        counter = 2^(k+1)-1;
      %  for i = 1:d-k
          %  S(1+counter:counter+2^(i+1),1+counter:counter+2^(i+1)) = [eye(2^i) -q*eye(2^i) ; zeros(2^i,2^i) eye(2^i)];
          %  counter = counter + 2^(i+1);
    %    end
        for i = 1:d-k %Generate S_1 through S_(d-1)
            for j = 1:2^(k-1)
                S(1+counter:counter+2^(i+1),1+counter:counter+2^(i+1)) = [eye(2^i) -q*eye(2^i) ; zeros(2^i,2^i) eye(2^i)];
                counter = counter + 2^(i+1);
            end
        end
        tm2 = S * tm2;
    end
    S_d = sym(zeros(2^(d+1)-1,2^(d+1)-1)); %Generate the final S_d
    S_d(1:3,1:3) = [ 1 -q 0 ; 0 1 -q ; 0 0 1];
    counter2 = 3;
    for i = 1:2^d-2
        S_d(1+counter2:2+counter2,1+counter2:2+counter2) = [1 -q ; 0 1];
        counter2 = counter2 + 2;
    end
    tm2 = S_d * tm2;
end

%This permutation matrix puts the diagonal form into SNF
function final_simplex_permutation = generate_final_simplex_permutation(d)
    last_index = 2^(d+1)-1;
    %index_list = generate_simplex_indices_v2(d); %Use this line to change
    %the algorithm for generating the final permutation matrix
    index_list = generate_simplex_indices(d);
    id = eye(last_index);
    final_simplex_permutation = zeros(last_index,last_index);
    for i = 1:last_index
        final_simplex_permutation(1:last_index,i) = id(1:last_index,index_list(i));
    end
end

%Support function for generating the final permutation matrix
%Counts the number of (1-q^2)^d terms present in the SNF for a simplex
function count = count_of_degrees(n,d) 
    count = nchoosek(d+1,n);
end 

%Support function for generating the final permutation matrix
%Converts a base10 number into base2
function base2 = convert_to_base2(input)
    if input == 0
        base2 = 0;
        return;
    end
    base2='';
    while input > 0
        r = string(mod(input,2));
        base2 = r + string(base2);
        input = floor(input / 2);
    end
end

%Support function for generating the final permutation matrix
%Calls convert_to_base2 and then counts the number of 1's
function count = count1s(input)
    input = char(convert_to_base2(input));
    count = sum(input == '1');
end 

% Support function for generating the final permutation matrix
function index_list = generate_simplex_indices(d)
    dict_counter = 0;
    dict = dictionary; %This dictionary stores the number of 1s found in a binary representation as a key and the indices corresponding to those number of twos as the value
    i = d+1;
    while i >= 1 %Populate dictionary
        dict(i) = {[dict_counter+1:dict_counter+nchoosek(d+1,i)]};
        dict_counter = dict_counter + nchoosek(d+1,i) ;
        i = i - 1;
    end
    index_list = []; %The indices of this list should be ordered where their values are
    for i = 1:2^(d+1)-1
        cell = dict(count1s(i)); %Cell is now a cell containing a list of indices that i will need to go into
        index_list(i) = cell{1,1}(find(cell{1,1} ~= 0, 1)); %Assigns the i-th element of the list to be the first nonzero element in the cell list
        cell{1,1}(find(cell{1,1} ~= 0, 1)) = 0; %Sets the first nonzero element in the cell list to 0 so that the index will not be reassigned
        dict(count1s(i)) = cell; %Updates the dictionary with the updated cell
    end
end

%Alternative support function
function index_list = generate_simplex_indices_v2(d) % Basically the same algorithm but does it backward
    dict_counter = 0;
    dict = dictionary; %This dictionary stores the number of 1s found in a binary representation as a key and the indices corresponding to those number of twos as the value
    i = d+1;
    while i >= 1 %Populate dictionary
        dict(i) = {[dict_counter+1:dict_counter+nchoosek(d+1,i)]};
        dict_counter = dict_counter + nchoosek(d+1,i) ;
        i = i - 1;
    end
    index_list = []; %The indices of this list should be ordered where their values are
    for i = 1:2^(d+1)-1
        cell = dict(count1s(i)); %Cell is now a cell containing a list of indices that i will need to go into
        index_list(i) = cell{1,1}(end); %Assigns the i-th element of the list to be the last element in the cell list
        cell{1,1}(end) = []; %Removes the last element from the cell so that the index will not be reassigned
        dict(count1s(i)) = cell; %Updates the dictionary with the updated cell
    end
end
