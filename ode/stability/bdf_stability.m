figure(1);
hold off;

f_bdf1=@(w) (exp(i*w) - 1.0) ./ -(exp(i*w));
w=linspace(-pi,pi,101);
bdf1_stability=f_bdf1(w);
plot(bdf1_stability,'g-+'); % - means line mark
hold on;

f_bdf2=@(w) -1.5 .* (exp(i*2*w) - 4.0./3.*exp(i*w) + 1.0./3) ./ exp(i*2*w);
w=linspace(-pi,pi,101);
bdf2_stability=f_bdf2(w);
plot(bdf2_stability,'b-+'); % - means line mark
hold on;

f_bdf3=@(w) -11/6 .* (exp(i*3*w) - 18.0./11.*exp(i*2*w) + 9.0./11.*exp(i*w) - 2.0./11) ./ exp(i*3*w);
w=linspace(-pi,pi,101);
bdf3_stability=f_bdf3(w);
plot(bdf3_stability,'k-+'); 
hold on;

f_bdf4=@(w) -25/12 .* (exp(i*4*w) - 48.0./25.*exp(i*3*w) + 36.0./25.*exp(i*2*w) - 16.0./25.*exp(i*w) + 3/25) ./ exp(i*4*w);
w=linspace(-pi,pi,101);
bdf4_stability=f_bdf4(w);
plot(bdf4_stability,'r-+'); 
hold on;

saveas(1,'figure/bdf.jpg');
