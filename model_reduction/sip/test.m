%
format long e;

% r1 1 2 1
% c1 2 0 1
% r2 2 3 1
% c2 3 0 1
% r3 3 4 1
G  = [...
1 -1 0 0;...
-1 2 -1 0;...
0 -1 2 -1;...
0 0 -1 1;...
];
C  = [...
0 0 0 0;...
0 2 -1 0;...
0 -1 2 0;...
0 0 0 0;...
];

% xi = [2 3]
% xc = [1 4]
P = [...
0 1 0 0;...
0 0 1 0;...
1 0 0 0;...
0 0 0 1;...
];

% diagonal pivoting
Gp = P * G * P';
Cp = P * C * P';
fprintf('Gp=\n');
Gp
fprintf('Cp=\n');
Cp
Gii = Gp(1:2,1:2);
Cii = Cp(1:2,1:2);
Gic = Gp(1:2,3:4);
Cic = Cp(1:2,3:4);
Gcc = Gp(3:4,3:4);
Ccc = Cp(3:4,3:4);
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
C_ext + C_neg

M = [...
-(Gii \ Gic);...
eye(2,2);...
];

M'*Gp*M
M'*Cp*M
