clear
syms l Izz Iyy A
k = [		A*l*l,       0,         0,       0,         0,  -A*l*l,        0,         0,        0,         0;
0,  12*Izz,   6*Izz*l,       0,         0,       0,  -12*Izz,   6*Izz*l,        0,         0;
0, 6*Izz*l, 4*Izz*l*l,       0,         0,       0, -6*Izz*l, 2*Izz*l*l,        0,         0;
0,       0,         0,  12*Iyy,   6*Iyy*l,       0,        0,         0,  -12*Iyy,   6*Iyy*l;
0,       0,         0, 6*Iyy*l, 4*Iyy*l*l,       0,        0,         0, -6*Iyy*l, 2*Iyy*l*l;
-A*l*l,       0,         0,       0,         0,   A*l*l,        0,         0,        0,         0;
0, -12*Izz,  -6*Izz*l,       0,         0,       0,   12*Izz,  -6*Izz*l,        0,         0;
0, 6*Izz*l, 2*Izz*l*l,       0,         0,       0, -6*Izz*l, 4*Izz*l*l,        0,         0;
0,       0,         0, -12*Iyy,  -6*Iyy*l,       0,        0,         0,   12*Iyy,  -6*Iyy*l;
0,       0,         0, 6*Iyy*l, 2*Iyy*l*l,       0,        0,         0, -6*Iyy*l, 4*Iyy*l*l;];

K=k;
K(10:15,10:15)=0;
K(6:15,6:15) = K(6:15,6:15)+k;
display(k);
display(K);