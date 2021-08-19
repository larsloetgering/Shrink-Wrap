function y = rect( x )
% function y = rect( x )
% note that discontinuities have "sharp edges"

x = abs( x );
y = double( x < 1/2 );

end