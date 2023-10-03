function out = Mh(h,D)
% Setup finite element matrix Mh
% h are element lengths
% D is true for dirichlet at the right
% M is true for use f mass lumping
    k=length(h);
    %out= spdiags([[h/6;0],[h/3;0]+[0;h/3], [0;h/6]], [-1,0,1], k+1,k+1);
    out= spdiags([[h/2;0]+[0;h/2]], [0], k+1,k+1);
end