*
.param res='1' z0=1 f='1'
p1 1 0 port=1 z0='2'
r1 1 2 'res'
p2 2 0 port=2 z0='3'
.ac lin 1 'f' 'f'
.lin sparcalc=1 format=touchstone filename=output dataformat=ri
.end
