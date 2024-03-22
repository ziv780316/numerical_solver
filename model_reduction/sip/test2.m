%
format long e;

% r12 1 2 1
% c11 2 0 c11
% r22 2 0 r22
% r23 2 3 1
% c22 3 0 c22
% r34 3 4 1
% c23 2 3 c23

g22 = 0;
g23 = 1;
G  = [...
1 -1 0 0;...
-1 1+g22+g23 -g23 0;...
0 -g23 g23+1 -1;...
0 0 -1 1;...
];

c12 = 0;
c22 = 1;
c23 = 0;
c33 = 2;
C  = [...
c12 -c12 0 0;...
-c12 c22+c23+c12 -c23 0;...
0 -c23 c33+c23 0;...
0 0 0 0;...
];

% xi = [2]
% xc = [1 3 4]
P = [...
0 1 0 0;...
1 0 0 0;...
0 0 1 0;...
0 0 0 1;...
];

% diagonal pivoting
Gp = P * G * P';
Cp = P * C * P';
fprintf('Gp=\n');
Gp
fprintf('Cp=\n');
Cp
Gii = Gp(1:1,1:1);
Cii = Cp(1:1,1:1);
Gic = Gp(1:1,2:4);
Cic = Cp(1:1,2:4);
Gcc = Gp(2:4,2:4);
Ccc = Cp(2:4,2:4);
fprintf('Gii =\n');
Gii
fprintf('Cii=\n');
Cii
fprintf('Gic=\n');
Gic
fprintf('Cic=\n');
Cic
fprintf('Gcc =\n');
Gcc
fprintf('Ccc=\n');
Ccc

[G_ext, C_ext, C_neg] = sip( Gii, Cii, Gic, Cic, Gcc, Ccc );
fprintf('G_ext=\n');
G_ext
fprintf('C_ext=\n');
C_ext
fprintf('C_neg=\n');
C_neg
fprintf('C_all=\n');
C_all = C_ext + C_neg

M = [...
-(Gii \ Gic);...
eye(3,3);...
];

fprintf('M''*Gp*M=\n');
G_sip = M'*Gp*M
fprintf('M''*Cp*M=\n');
C_sip = M'*Cp*M

fprintf('|G_sip-G_ext|=%.10e\n', max(max(G_sip-G_ext)));
fprintf('|C_sip-C_all|=%.10e\n', max(max(C_sip-C_all)));
