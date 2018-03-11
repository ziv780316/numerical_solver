figure(1);
hold off;

f_am1=@(w) (exp(i*w) - 1.0) ./ -exp(i*w);
w=linspace(-pi,pi,101);
am1_stamility=f_am1(w);
plot(am1_stamility,'b-+'); % - means line mark
hold on;

%f_am2=@(w) (exp(i*w) - 1.0) ./ -((1.0/2.0).*exp(i*w) + (1.0/2.0));
%w=linspace(-pi,pi,101);
%am2_stamility=f_am2(w);
%plot(am2_stamility,'r-+'); % - means line mark
%hold on;

f_am3=@(w) (exp(i*2*w) - exp(i*w)) ./ -((5.0/12.0).*exp(i*2*w) + (2.0/3.0).*exp(i*w) - (1.0/12.0));
w=linspace(-pi,pi,101);
am3_stamility=f_am3(w);
plot(am3_stamility,'k-+'); % - means line mark
hold on;

f_am4=@(w) (exp(i*3*w) - exp(i*2*w)) ./ -((3.0/8.0).*exp(i*3*w) + (19.0/24.0).*exp(i*2*w) - (5.0/24.0).*exp(i*w) + (1.0/24.0));
w=linspace(-pi,pi,101);
am4_stamility=f_am4(w);
plot(am4_stamility,'r-+'); % - means line mark
hold on;

axis equal;

saveas(1,'figure/am.jpg');
