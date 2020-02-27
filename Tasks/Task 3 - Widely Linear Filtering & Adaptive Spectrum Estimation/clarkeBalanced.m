function Transform = clarkeBalanced(V,N)
%Project three-phase system onto two orthogonal axis alpha and beta
ClarkeMatrix = sqrt(2/3) * [sqrt(2)/2 sqrt(2)/2 sqrt(2)/2; 1 -1/2 -1/2; 0 sqrt(3)/2 -sqrt(3)/2];
Transform = zeros(3,N);
for i =1:N
    Transform(:,i) = ClarkeMatrix * V(:,i);
end

