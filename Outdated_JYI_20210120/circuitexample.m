A=[1, 0
   -2, 2+1j+4   ];

y=[4 
   0 ]

x=inv(A)*y

voc=x(2,1)*4

angle(voc)*360/(2*pi)