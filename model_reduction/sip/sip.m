%
function [G_ext, C_ext, C_neg] = sip( Gii, Cii, Gic, Cic, Gcc, Ccc )
% G_ext = Gcc - Gii \ (Gic * Gic')
% C_ext = Ccc + Gii \ (Cic * Gic' + Gic * Cic') - (Gii \ Cii) * Gic * (Gii \ Gic')

G_ext = Gcc - Gic' * (Gii \ Gic);
C_ext = Ccc -Gii \ (Cic * Gic' + Gic * Cic');
C_neg = Gic' * (Gii \ Cii) * (Gii \ Gic);

end
