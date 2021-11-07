function C = EP2C(q)
% EP2C	
%
%	C = EP2C(Q) returns the direction cosine 
%	matrix in terms of the 4x1 Euler parameter vector
%	Q.  The first element is the non-dimensional Euler
%	parameter, while the remain three elements form 
%	the Eulerparameter vector.
%
qs = q(1); qv = q(2:4);
C = (qs^2 - qv'*qv)*eye(3) + 2*qv*(qv') - 2*qs*crossProductMatrix(qv);
