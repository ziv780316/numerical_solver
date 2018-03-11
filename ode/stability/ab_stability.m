figure(1);
hold off;

f_ab1=@(w) (exp(i*w) - 1.0) ./ -(1.0);
w=linspace(-pi,pi,101);
ab1_stability=f_ab1(w);
plot(ab1_stability,'b-+'); % - means line mark
hold on;

f_ab2=@(w) (exp(i*2*w) - exp(i*w)) ./ ((-3.0/2.0)*exp(i*w) + (1.0/2.0));
w=linspace(-pi,pi,101);
ab2_stability=f_ab2(w);
plot(ab2_stability,'r-+'); % - means line mark
hold on;

f_ab3=@(w) (exp(i*3*w) - exp(i*2*w)) ./ ((-23.0/12.0)*exp(i*2*w) + (4.0/3.0)*exp(i*w) - (5.0/12.0));
w=linspace(-pi,pi,101);
ab3_stability=f_ab3(w);
plot(ab3_stability,'k-+'); % - means line mark
hold on;

axis equal;

saveas(1,'figure/ab.jpg');
