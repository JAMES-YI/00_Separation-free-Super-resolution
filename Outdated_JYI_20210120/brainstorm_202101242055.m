%% circuit example
% - from circuitexample.m
A=[1, 0
   -2, 2+1j+4   ];

y=[4 
   0 ]

x=inv(A)*y

voc=x(2,1)*4

angle(voc)*360/(2*pi)

%% circuit homework
% - from circuithomework
A=[1, 0, 0
   -10, 20+20*1j, -10
   0, -10, -5*1j+10
   ]

y=[10*exp(1j* 30/360*2*pi)
   0
   0]
    
x=inv(A)*y

