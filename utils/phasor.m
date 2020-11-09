function R = phasor(omega, D)

if nargin < 2
    D = 2;
end

if D < 2
    disp('Dimension of phasor must be at least 2.');
    R = cos(omega);
    return;
end

R = zeros(D);

R(1:2,1:2) = [cos(omega) -sin(omega)
              sin(omega)  cos(omega)];

end

