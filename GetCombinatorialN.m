function [ output ] = GetCombinatorialN( small_x, big_y )
% C_y^x
output = (prod(1:big_y)/(prod(1:small_x)*prod(1:big_y-small_x)));

end

