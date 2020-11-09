function out = laprnd(mu,b,N)

dimy = size(mu,1);

if nargin < 3
    N = 1;
end

U = rand(dimy,N)-1/2;
out = mu - b*sign(U).*log(1-2*abs(U));