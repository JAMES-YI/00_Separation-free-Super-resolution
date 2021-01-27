format long
A=[-1,0;0,0]; 
B=[1,1;1,1];
t=-1;
for i=1:1
sum(svd(A+t*B))
end
