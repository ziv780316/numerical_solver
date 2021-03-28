*
.param cap='1' z0=1 pi=3.1415926535 f='1/(2*pi)'
p1 1 0 port=1 z0='1'
c1 1 2 'cap'
p2 2 0 port=2 z0='1'
.ac lin 1 'f' 'f'
.lin sparcalc=1 format=touchstone filename=output dataformat=ri
.end
