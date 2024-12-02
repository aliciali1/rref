
function [M, piv] = reduce(M, verbose)
% REDUCE Given a matrix, perform row echelon reduction to reduced echelon
% form
%
% function [M, piv] = reduce(M, verbose)
%
% Inputs: 
%   M: matrix to perform row echelon reduction on
%   verbose: prints actions to the console when true
% Outputs: 
%   M: matrix in reduced echelon form
%   piv: stores the indices of the pivot columns in the matrix

arguments 
    M 
    verbose = true
end
piv = [];
c = 1;
r = 1;
submatrix = M;
while size(submatrix,1) > 1
    [col_maxes,max_indices] = max(submatrix,[],1);
    col = find(col_maxes,1,'first');
    to_swap = max_indices(col);
    if 1 ~= to_swap
        M = exchange(M, r, r+to_swap-1, verbose);
    end
    submatrix = M(r:end,c:end);
    piv = [piv, c+col-1];
    for row = 2:size(submatrix,1)
        scale = (submatrix(row,col)/submatrix(1,col))*-1;
        if scale ~= 0
            M = add(M, scale,r, r+row-1, verbose);
            submatrix = M(r:end,c:end);
        end
    end
    M(r:end,c:end) = submatrix;
    r = r + 1;
    c = col + 1;
    submatrix = M(r:end, c:end);
end
last_piv = find(M(end,:),1,'first');
piv = [piv,last_piv];
% now, all the pivots are found and the matrix is in echelon form
% use a loop to work backwards with all the rows
for rr = size(M,1):-1:1
    pivot = find(M(rr,:),1,'first');
    if ~isempty(pivot)
        scale = 1/M(rr,pivot);
        M = mult(M,scale,rr,verbose);
        for row_to_reduce = 1:rr-1
            if M(row_to_reduce, pivot) ~= 0
                scale2 = -1*(M(row_to_reduce,pivot)/M(rr,pivot));
                if scale2 ~= 0
                    M = add(M, scale2, rr, row_to_reduce,verbose);
                end
            end
        end
    end
end
end %End of function reduce




function M = exchange(M, row1, row2, verbose)
% EXCHANGE Swap two rows in a matrix
%
% function M = exchange(M, row1, row2, verbose)
%
% Inputs: 
%   M: matrix to perform exchange on
%   row1: one row to be swapped
%   row2: the other row to be swapped
%   verbose: prints actions to the console when true
% Outputs: 
%   M: matrix after performing exchange
arguments 
    M
    row1 {mustBeInteger, mustBePositive}
    row2 {mustBeInteger, mustBePositive}
    verbose = true
end

if row1 > size(M,1)
    error('Row %d not valid in this matrix\n', row1)
elseif row2 > size(M,1)
    error('Row %d not valid in this matrix\n', row2)
end

M([row1 row2], :) = M([row2 row1], :);

if verbose
    fprintf('Exchanging rows %d and %d\n', row1, row2)
end
end %End of function exchange




function M = mult(M, d, row, verbose)
% MULT multiply every element in a row by a value
%
% function M = mult(M, d, row, verbose)
%
% Inputs: 
%   M: matrix to perform mult on
%   d: factor to multiply elements by
%   row: row to be multiplied
%   verbose: prints actions to the console when true
% Outputs: 
%   M: matrix after performing mult
arguments 
    M
    d {mustBeNonzero}
    row {mustBeInteger, mustBePositive}
    verbose = true
end

if row > size(M,1)
    error('Row %d not valid in this matrix\n', row)
end

M(row,:) = d*M(row,:);

if verbose
    fprintf('Multiplying row %d by %f\n', row, d)
end
end %End of function mult



function M = add(M, r, row1, row2, verbose)
% ADD set row2 to be the sum of r times row1 and row2
%
% function M = add(M, r, row1, row2, verbose)
%
% Inputs: 
%   M: matrix to perform add on
%   r: factor to multiply row1 by
%   row1: the row that gets multiplied
%   row2: the row that gets replaced with the sum of itself and r times
%   row1
%   verbose: prints actions to the console when true
% Outputs: 
%   M: matrix after performing add
arguments 
    M
    r
    row1 {mustBeInteger, mustBePositive}
    row2 {mustBeInteger, mustBePositive}
    verbose = true
end

if row1 > size(M,1)
    error('Row %d not valid in this matrix', row1)
elseif row2 > size(M,1)
    error('Row %d not valid in this matrix', row2)
end

M(row2,:) = M(row2,:)+r*M(row1,:);

if verbose
    fprintf('Adding %f times row %d to row %d \n', r, row1, row2)
end
end %End of function add


% Test results

% Test 1
% A = 3     1     5     5    -1     2     2     2    -2     2
%     4    -4     5     0     5    -5     3    -4    -5    -2
%    -4    -2     -4    3     3     4     3     2    -4     5
%     5     1      5   -4     5     5    -1     -5    4    -5

% [R, piv] = reduce(A)
% Exchanging rows 1 and 4
% Adding -0.800000 times row 1 to row 2 
% Adding 0.800000 times row 1 to row 3 
% Adding -0.600000 times row 1 to row 4 
% Exchanging rows 2 and 4
% Adding 3.000000 times row 2 to row 3 
% Adding 12.000000 times row 2 to row 4 
% Exchanging rows 3 and 4
% Adding -0.240000 times row 3 to row 4 
% Multiplying row 4 by -12.500000
% Adding 4.000000 times row 4 to row 1 
% Adding -7.400000 times row 4 to row 2 
% Adding -92.000000 times row 4 to row 3 
% Multiplying row 3 by 0.040000
% Adding -5.000000 times row 3 to row 1 
% Adding -2.000000 times row 3 to row 2 
% Multiplying row 2 by 2.500000
% Adding -1.000000 times row 2 to row 1 
% Multiplying row 1 by 0.200000
% 
% R =
% 
%   Columns 1 through 7
% 
%     1.0000         0         0         0 -350.2500 -563.2500  -91.5000
%          0    1.0000         0         0    7.2500   14.2500    1.5000
%          0         0    1.0000         0  287.0000  461.0000   75.0000
%          0         0         0    1.0000  -78.5000 -125.5000  -20.0000
% 
%   Columns 8 through 10
% 
%    75.2500  -33.0000  -66.5000
%    -1.2500    2.0000    1.5000
%   -62.0000   27.0000   54.0000
%    17.5000   -8.0000  -14.0000
% 
% 
% piv =
% 
%      1     2     3     4




% Test 2
% A =
    % -4    -1     4    -1
    %  5     5     1    -5
    % -5    -3    -2     4
    %  3    -3     0     5
    %  3    -4    -1     0
    %  4    -4    -5     0
    % -5     4    -3    -2
    % -1     1    -4     4
    % -3     1    -3    -1
    %  3    -4    -3    -4
%
% [R, piv] = reduce(A)
% Exchanging rows 1 and 2
% Adding 0.800000 times row 1 to row 2 
% Adding 1.000000 times row 1 to row 3 
% Adding -0.600000 times row 1 to row 4 
% Adding -0.600000 times row 1 to row 5 
% Adding -0.800000 times row 1 to row 6 
% Adding 1.000000 times row 1 to row 7 
% Adding 0.200000 times row 1 to row 8 
% Adding 0.600000 times row 1 to row 9 
% Adding -0.600000 times row 1 to row 10 
% Exchanging rows 2 and 7
% Adding -0.222222 times row 2 to row 3 
% Adding 0.666667 times row 2 to row 4 
% Adding 0.777778 times row 2 to row 5 
% Adding 0.888889 times row 2 to row 6 
% Adding -0.333333 times row 2 to row 7 
% Adding -0.222222 times row 2 to row 8 
% Adding -0.444444 times row 2 to row 9 
% Adding 0.777778 times row 2 to row 10 
% Exchanging rows 3 and 7
% Adding 0.353659 times row 3 to row 4 
% Adding 0.577236 times row 3 to row 5 
% Adding 1.386179 times row 3 to row 6 
% Adding 0.101626 times row 3 to row 7 
% Adding 0.613821 times row 3 to row 8 
% Adding 0.276423 times row 3 to row 9 
% Adding 0.943089 times row 3 to row 10 
% Adding 4.000000 times row 4 to row 6 
% Exchanging rows 5 and 6
% Adding 1.093750 times row 5 to row 6 
% Adding -0.078125 times row 5 to row 7 
% Adding -0.801339 times row 5 to row 8 
% Adding 0.446429 times row 5 to row 9 
% Adding 2.459821 times row 5 to row 10 
% Multiplying row 5 by 0.274554
% Adding 5.000000 times row 5 to row 1 
% Adding 7.000000 times row 5 to row 2 
% Adding 2.666667 times row 5 to row 3 
% Adding -2.390244 times row 5 to row 4 
% Multiplying row 4 by 4503599627370496.000000
% Adding -1.000000 times row 4 to row 1 
% Adding 2.000000 times row 4 to row 2 
% Adding -5.466667 times row 4 to row 3 
% Multiplying row 2 by 0.111111
% Adding -5.000000 times row 2 to row 1 
% Multiplying row 1 by 0.200000
% 
% R =
% 
%      1     0     0     0
%      0     1     0     0
%      0     0     0     0
%      0     0     1     0
%      0     0     0     1
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0
%      0     0     0     0


% Test 3
% Exchanging rows 2 and 3
% Multiplying row 3 by 0.500000
% Adding -1.000000 times row 3 to row 2 
% Multiplying row 2 by 0.500000
% Multiplying row 1 by 1.000000
% 
% R =
% 
%     1.0000    2.0000         0         0         0
%          0         0    1.0000         0    0.7500
%          0         0         0    1.0000    1.5000
% 
% 
% piv =
% 
%      1     3     4


% Test 4
% Adding -1.000000 times row 1 to row 2 
% Adding -1.000000 times row 1 to row 3 
% Adding -1.000000 times row 1 to row 4 
% Multiplying row 1 by 1.000000
% 
% R =
% 
%      1     2     4     3     5
%      0     0     0     0     0
%      0     0     0     0     0
%      0     0     0     0     0
% 
% 
% piv =
% 
%      1

