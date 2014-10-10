%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% coefs = fit_closed_b_spline(s,x,L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function fits a closed B-spline with L knots to the input
% data (s,x), and returns the control points.  The B-spline chosen
% is optimal in the least-squares sense.  The control points coefs
% chosen minimize \sum_{i=1}^N trace((x_i - f(s_i ; coefs)') * (x_i
% - f(s_i ; coefs))).  
%
% Input:
% s: An array of size N x 1, indicating where on the B-spline the
% data points x lie.
% x: An N x ndims matrix of data points to be fit by the B-spline,
% that is x_i \approx \sum_{i=0}^{L-1} B_i(s_i) * coefs(i+1,:).
% L: Number of curve pieces in the B-spline.  
%
% Output:
% coefs: An (L+3) x ndims matrix containing the control points of the
% curve.  coefs(i,:) is the i-1 th control point.  Note that the
% coefs(L+1,:) = coefs(1,:), coefs(L+2,:) = coefs(2,:), and
% coefs(L+3,:) = coefs(3,:).
%
% This function first constructs a N x L matrix A s.t. 
% A(i,j) = B_{j-1}(s(i)), so the value of the B-spline with control
% points p is f(s(i) ; p) = A(i,:)*p, and f(s ; p) = A * p.  To
% find the set of control points that minimizes \sum_{i=1}^N
% trace((x_i - f(s_i ; coefs)') * (x_i - f(s_i ; coefs)))
% = trace((x - A * coefs)'*(x-A*coefs)), we can take the derivative
% w.r.t. coefs and set it to 0:
% zeros(L,ndims) = - 2 * A' * (x - A * coefs)
% A' * x = A' * A * coefs
% inv(A'*A) * A' * x = coefs
%
% The matrix coefs is currently only L x ndims, and contains only
% control points 0 through L - 1.  We augment it with control
% points L, L+1, and L+2, which, because the spline is cyclic, are
% equivalent to control points 0, 1, and 2, respectively.
%
function [coefs,A] = fit_closed_b_spline(s,x,L)

% Hard coded to be cubic B-spline
k = 4;

% x is N x ndims
[N,ndims] = size(x);

s = mod(s,L);

%%% Create a matrix A s.t. f(s) = A * coefs.
% A is N x L, and entry A(i,j) = B_{j-1}(s(i)), f(s(i)) = sum(A(i,:)*coefs).
% coefs is L x ndims, coefs(i) is the (i-1)th control point.  Note
% that the Lth control point is coefs(1), the L+1th control point
% is coefs(2) and the L+2th control point is coefs(3).  

% oneinds(i)-1 <= s(i) < oneinds(i)-1+1
oneinds = floor(s)+1;
% twoinds(i)-1+1 <= s(i) < twoinds(i)-1+2
twoinds = mod(oneinds-2,L)+1;
% threeinds(i)-1+2 <= s(i) < threeinds(i)-1+3
threeinds = mod(oneinds-3,L)+1;
% fourinds(i)-1+3 <= s(i) < fourinds(i)-1+4
fourinds = mod(oneinds-4,L)+1;

% The input to the basis function
d = mod(s,1);

% Five different functions, depending on where s lies in relation
% to sigma

% if j <= s(i) < j+1
A1 = d.^3/6;
% if j+1 <= s(i) < j+2
A2 = (-3*d.^3 + 3*d.^2 + 3*d + 1)/6;
% if j+2 <= s(i) < j+3
A3 = (3*d.^3 - 6*d.^2 + 4)/6;
% if j+2 <= s(i) < j+3
A4 = (1-d).^3/6;
% otherwise
%A0 = zeros(N,L);

A = zeros(N,L);
A(sub2ind([N,L],(1:N)',oneinds)) = A1;
A(sub2ind([N,L],(1:N)',twoinds)) = A2;
A(sub2ind([N,L],(1:N)',threeinds)) = A3;
A(sub2ind([N,L],(1:N)',fourinds)) = A4;

% find coefs that minimizes trace((x - f(s))'*(x - f(s)))
% = trace((x - A * coefs)'*(x - A * coefs)).

% take the derivative w.r.t. coefs, and set to 0:
% zeros(L,ndims) = - 2 * A' * (x - A * coefs)
% A' * x = A' * A * coefs
% inv(A'*A) * A' * x = coefs
coefs = pinv(A) * x;

% Augment with control points L, L+1, and L+2
coefs = [coefs;coefs(1:3,:)];
