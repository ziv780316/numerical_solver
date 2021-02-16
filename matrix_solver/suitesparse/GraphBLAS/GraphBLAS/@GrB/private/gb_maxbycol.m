function C = gb_maxbycol (op, A)
%GB_MAXBYCOL max, by column
% Implements C = max (A, [ ], 1)

% C = max (A, [ ], 1) reduces each col to a scalar; C is 1-by-n
desc.in0 = 'transpose' ;
C = gbvreduce (op, A, desc) ;

% if C(j) < 0, but if A(:,j) is sparse, then assign C(j) = 0.
ctype = gbtype (C) ;

if (gb_issigned (ctype))
    % d (j) = number of entries in A(:,j); d (j) not present if A(:,j) empty
    [m, n] = gbsize (A) ;
    d = gbdegree (A, 'col') ;
    % d (j) is an explicit zero if A(:,j) has 1 to m-1 entries
    d = gbselect (d, '<', m) ;
    zero = gbnew (0, ctype) ;
    if (gbnvals (d) == n)
        % all columns A(:,j) have between 1 and m-1 entries
        C = gbapply2 (op, C, zero) ;
    else
        d = gbapply2 (['1st.' ctype], zero, d) ;
        % if d (j) is between 1 and m-1 and C (j) < 0 then C (j) = 0
        C = gbeadd (op, C, d) ;
    end
end

C = gbtrans (C) ;

