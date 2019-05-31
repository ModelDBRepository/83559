function listing = list_of_pairs(v1,v2,varargin)

% LIST_OF_PAIRS generate all pairs of numbers
%   LIST_OF_PAIRS(V1,V2,FLAG) generates a list of all pairs of the numbers in
%   vectors V1 and V2. Optional FLAG is a sequence of letters setting options:
%       's' to remove same-number pairs from the returned list
%       'u' to ensure all pairs are unique (e.g. so that only one of [1 2]
%       & [2 1] are retained)
%
%
%   NOT FINISHED: nchoosek(V,N) does similar thing with nchoosek(    
%
%   Mark Humphries 15/12/2005

[r1 c1] = size(v1);
n1 = max(r1,c1);

[r2 c2] = size(v2);
n2 = max(r2,c2);

if r1 == 1
    first_col = repmat(v1,n2,1);
else
    first_col = repmat(v1',n2,1);
end

if r2 == 1
    second_col = repmat(v2',n1,1);
else
    second_col = repmat(v2,n1,1);
end

listing = [first_col(:) second_col];

if nargin >= 2
    if findstr(varargin{1},'u')
        isect = intersect(v1,v2);   % the numbers that appear in both vectors
        
        if length(isect) > 1
            del_list = combs(isect,2);      % returns all unique combinations of the numbers in isect
            [r c] = size(del_list);
            % delete all occurences of those pairs in the list
            for loop = 1:r
                ix = find(listing(:,1) == del_list(loop,1) & listing(:,2) == del_list(loop,2));
                listing(ix,:) = [];
            end
        end
    end
    if findstr(varargin{1},'s')
        ix = find(listing(:,1) == listing(:,2));
        listing(ix,:) = [];
    end
end

%----------------------------------------
function P = combs(v,m)
%COMBS  All possible combinations.
%   COMBS(1:N,M) or COMBS(V,M) where V is a row vector of length N,
%   creates a matrix with N!/((N-M)! M!) rows and M columns containing
%   all possible combinations of N elements taken M at a time.
%
%   This function is only practical for situations where M is less
%   than about 15.

if nargin~=2, error('MATLAB:nchoosek:WrongInputNum', 'Requires 2 input arguments.'); end

v = v(:).'; % Make sure v is a row vector.
n = length(v);
if n == m
   P = v;
elseif m == 1
   P = v.';
else
   P = [];
   if m < n && m > 1
      for k = 1:n-m+1
         Q = combs(v(k+1:n),m-1);
         P = [P; [v(ones(size(Q,1),1),k) Q]];
      end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

