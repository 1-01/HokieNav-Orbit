% Example 3.16, Davenport's q-method example

yaw_pitch_roll = [30,20,-10]*pi/180;

BN = Euler3212C(yaw_pitch_roll) ;

s_N = [1 0 0]';
m_N = [0 0 1]';

s_B_true = BN*s_N ;
m_B_true = BN*m_N ;

s_B = [ 0.8190 -0.5282 0.2242]';
s_B = s_B/norm(s_B);

m_B = [-0.3138 -0.1584 0.9362]';
m_B = m_B/norm(m_B);

v1_N=s_N;
v2_N=m_N;

v1_B=s_B;
v2_B=m_B;

w1=1;
w2=1;

B = w1*v1_B*v1_N' + w2*v2_B*v2_N' ;

S = B + B' ;

sigma = trace(B) ;

Z = [ B(2,3) - B(3,2) ; ...
      B(3,1) - B(1,3) ; ...
      B(1,2) - B(2,1) ] ;
  
K = [ sigma        Z'       ; ...
        Z    S-sigma*eye(3) ];

[eigvec,eigval]=eigs(K) ;

beta = eigvec(:,2);

barBN=EP2C(beta) ;

barBB = barBN*BN';

ErrorPhi = acos(0.5*(trace(barBB)-1))*180/pi ;