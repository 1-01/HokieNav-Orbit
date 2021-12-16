function [q, R] = davenportq(vA, vB, w)
% 
% Matt Werner (m.werner@vt.edu) - Nov 3, 2021
% 
% Given a set of measurements in frame B and corresponding equivalent
% representations of the measurements (i.e. from a model) in frame A, find
% the optimal rotation matrix that provides a coordinate rotation from
% frame B to frame A such that vA = R*vB for each vector contained in vB.
% 
% That is, minimize the cost function
%                             __
%                             \                    2
%                        J =  /_  w  ||vA - R*vB ||
%                              k   k     k      k
% 
% given a set of k measurements in frame B and equivalent representations
% in frame A for the rotation matrix R. Here, each vector vAk, vBk must be
% unit normalized so that vAk'*vAk = vBk'*vBk = 1.
% 
% This procedure is done by Davenport's q-method which parameterizes R
% using a quaternion. The optimal rotation matrix R corresponds to the 
% eigenvector paired with the most positive eigenvalue of the 4-by-4
% Davenport matrix K defined
%                              _               _
%                             |   s       z'    |
%                        K =  |                 |
%                             |_  z    S - s*I _|
% with
% 
%       s = trace B 
% 
%       S = B + B'
%            _           _
%       z = |  B   - B    |
%           |   23    32  |
%           |             |
%           |  B   - B    |
%           |   31    13  |
%           |             |
%           |  B   - B    |
%           |_  12    21 _|
%            __
%            \
%       B =  /_  wk*vAk*vBk'
%             k
% 
%      wk = chosen weights for each measurement/model pair.
%
% Such an eigenvector is the quaternion that parameterizes the optimal
% rotation matrix R. Note that the quaternion is expressed
%                                      _    _
%                                     |  qr  |
%       q = qr + qi*i + qj*j + qk*k = |      |,
%                                     |  qi  |
%                                     |      |
%                                     |  qj  |
%                                     |      |
%                                     |_ qk _|
% 
% where i*i = j*j = k*k = i*j*k = -1.
% 
% As such, the scalar part of the quaternion occupies the first element
% while the vector part of the quaternion occupies the last 3 elements.
% 
%    Inputs:
% 
%                vA - A set of 3-by-1 vectors expressed in the fixed frame
%                     (frame A). Each vector (k of them) must be positioned 
%                     adjacent to each other so that the resulting input is 
%                     3-by-k with k > 1. Each vector should be unit
%                     normalized.
%                     Size: 3-by-k
%                     Units: SI
% 
%                vB - A set of 3-by-1 measurement vectors taken in the body
%                     frame (frame B). Each measurement vector (k of them)
%                     must be positioned adjacent to each other so that the
%                     resulting input is 3-by-k with k > 1. Each vector
%                     should be unit normalized.
%                     Size: 3-by-k
%                     Units: SI
% 
%                 w - (Optional) Weights to be associated with each
%                     measurement. The weights should be given such that
%                     the sum of weights is unity; if not true, then the
%                     weights are divided by their sum. DEFAULT behavior is
%                     that each measurement receives equal weighting
%                     corresponding to
%                                   w = [1, 1, ..., 1]/N
%                     where N > 1 is the amount of vectors given in each
%                     set vA and vB.
% 
%    Outputs:
% 
%                 q - The optimal quaternion parameterizing the rotation
%                     matrix R such that the cost function J is minimized
%                     (see introduction description for more information).
%                     Size: 4-by-1 (vector)
%                     Units: - (N/A)
% 
%                 R - The rotation matrix resulting from expressing the
%                     action of quaternions in the form of a linear
%                     transformation that transforms from the body frame
%                     (frame B) to the fixed frame (frame A).
%                     Size: 3-by-3 (matrix)
%                     Units: - (N/A)
% 
% Sources: 
%   - https://ahrs.readthedocs.io/en/latest/filters/davenport.html
%   - https://math.stackexchange.com/questions/1634113/davenports-q-method-finding-an-orientation-matching-a-set-of-point-samples
%   - S.D. Ross
% 

%% Checks
% Check that the number of inputs is at least 2
narginchk(2, 3)
% Ensure that vA and vB are the same size
assert(all(size(vA) == size(vB)), "Input must be same size.")
% Check that vA contains 3 rows
assert(size(vA, 1) == 3, "Row dimension should be 3 but is %1.0f.", size(vA, 1))
% Check that N > 1
N = size(vA, 2);
assert(N > 1, "Provided sample size needs to be greater than 1.")
% Check that each vector in the sets are unit normalized. If not, then
% normalize them. (Only alter the inputs if necessary.)
for k = 1:N
    % Check each vector in frame A
    if (vA(:,k)'*vA(:,k) ~= 1)
        vA(:,k) = vA(:,k)/norm(vA(:,k));
    end
    % Check each vector in frame B
    if (vB(:,k)'*vB(:,k) ~= 1)
        vB(:,k) = vB(:,k)/norm(vB(:,k));
    end
end
% Check if the weights are given and are normalized. If not given, then
% provide the default setting that each measurement receives equal weight.
if (nargin == 2)
    % Create the default weights wk
    w = ones(1, N)/N;
else
    % Ensure that the provided weights are normalized
    if (sum(w) ~= 1)
        w = w/sum(w);
    end
end

%% Calculations
% Form the attitude profile matrix (B). Note that vA and vB are expressed
% in different frames; we treat them only as a tuple of numbers here
% independent of frames so that it makes sense to combine them in
% calculations directly.
B = zeros(3,3);
for k = 1:N
    B = B + w(k)*vA(:,k)*vB(:,k)';
end

% Form components of the Davenport matrix (K)
s = trace(B);
S = B + B';
z = [B(2,3) - B(3,2); B(3,1) - B(1,3); B(1,2) - B(2,1)];

% Calculate K and its eigenvalues/eigenvectors. Note that Matlab
% automatically normalizes its eigenvectors which satisfies the quaternion
% magnitude condition.
K = [s, z'; z, S-s*eye(3)];
[eigenvecs, eigenvals] = eig(K);

% Pick the eigenvector corresponding to the most positive eigenvalue. Since
% the Davenport matrix is symmetric, all the eigenvalues are real-valued,
% so simply take the maximum eigenvalue.
eigenvals = [eigenvals(1,1), eigenvals(2,2), eigenvals(3,3), eigenvals(4,4)];
[~, qcol] = max(eigenvals);
q = eigenvecs(:, qcol);
% Calculate the rotation matrix in terms of this quaternion
R = EP2C(q);