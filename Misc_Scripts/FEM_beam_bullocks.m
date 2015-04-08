%%
% * * ** * * %****************************************************************************
%****************************************************************************
%************************ BEAMANALYSIS ***********************************
%****************************************************************************
%************************* Written by ******************************************
%****************************************************************************
%******************* TARTIBU KWANDA 2008 **********************************
%****************************************************************************
%****************************************************************************
%****************************************************************************
%Uniform beam and Stepped beam finite element program. Selectable number of step and 
%number of elements. Solves for eigenvalues and eigenvectors of a beam with user defined 
%dimensions. This program is able to calculate the natural frequencies of different uniform and 
%stepped beam geometries and material’s properties.
%default values are included in the program for the purpose of showing how to input data.
    echo off
%     clf;
    clear all;
    inp = input('Input "1" to enter beam dimensions, "Enter" to use default ... ');
    if (isempty(inp))
        inp = 0;
    else
    end
    if inp == 0
%   input size of beam and material
%   xl(i) = length of element (step)i
%   w(i) = width of element (step)i
%   t(i) = thickness of element (step)i
%   e = Young's modulus
%   bj = global degree of freedom number corresponding to the jth local degree
%   of freedom of element i
%   a(i) = area of cross section of element i
%   ne = number of elements
%   n = total number of degree of freedom
%   no = number of nodes
    format short
    xl = [40 32 24];
    xi = [1.333333 6.75 0.08333];
    w = [4 6 2];
    t = [2 3 1];
    e = 206e+6;
    rho = 7.85e-6;
    bj = [1 2 3 4;3 4 5 6;5 6 7 8];
    ne = 3;
    n = 8;
    else
       

   
    ne = input ('Input number of elements, default 3 ... ');
    if (isempty(ne))
        ne = 3;
    else
    end
    
    xl = input ('Input lengths of stepped beam, default [40 32 24], ... ');
    if (isempty(xl))
        xl = [40 32 24];
    else
    end
    
    w = input ('Input widths of stepped beam, default [4 6 2], ... ');
    if (isempty(w))
        w = [4 6 2];
    else
    end
    
    t = input ('Input thickness of stepped beam, default [2 3 1], ... ');
    if (isempty(t))
        t = [2 3 1];
    else
    end
    
    e = input('Input modulus of material, mN/mm^2, default mild steel 206e+6 ... ');
    if (isempty(e))
        e=206e+6;
    else
    end
    
    rho = input('Input density of material,kg/mm^3 , default mild steel 7.85e-6 ... ');
    if (isempty(rho))
        rho = 7.85e-6;
    else
    end
    
    bj = input('Input global degree of freedom, default global degree of freedom [1 2 3 4;3 4 5 6;5 6 7 8] ... ');
    if (isempty(bj))
        bj = [1 2 3 4;3 4 5 6;5 6 7 8];
    else
    end
    end

 

%   Calculate area (a), area moment of Inertia (xi) and mass per unit of length (xmas) of
%   the stepped beam
    a = w.*t;
%   Define area moment of inertia according to flap-wise or edge-wise
%   vibration of the stepped beam.
    vibrationdirection = input('enter "1" for edge-wise vibration, "enter" for flap-wise vibration ... ');
    if (isempty(vibrationdirection))
        vibrationdirection = 0;
    else
    end
    if vibrationdirection == 0
        xi=w.*t.^3/12;
    else
        xi=t.*w.^3/12;
    end
    
    for i=1:ne
       xmas(i)=a(i)*rho;
    end
   


%   Size the stiffness and mass matrices
    no = ne+1;
    n = 2*no;  
bigm = zeros(n,n);
bigk = zeros(n,n); 
%   Now build up the global stiffness and consistent mass matrices, element
%   by element
for ii =1:ne
        ai=zeros(4,n);   
       i1=bj(ii,1);
       i2=bj(ii,2);
       i3=bj(ii,3);
       i4=bj(ii,4);
       ai(1,i1)=1;
       ai(2,i2)=1;
       ai(3,i3)=1;
       ai(4,i4)=1;
       xm(1,1)=156;
       xm(1,2)=22*xl(ii);
       xm(1,3)=54;
       xm(1,4)=-13*xl(ii);
       xm(2,2)=4*xl(ii)^2;
       xm(2,3)=13*xl(ii);
       xm(2,4)=-3*xl(ii)^2;
       xm(3,3)=156;
       xm(3,4)=-22*xl(ii);
       xm(4,4)=4*xl(ii)^2;
       xk(1,1)=12;
       xk(1,2)=6*xl(ii);
       xk(1,3)=-12;
       xk(1,4)=6*xl(ii);
       xk(2,2)=4*xl(ii)^2;
       xk(2,3)=-6*xl(ii);
       xk(2,4)=2*xl(ii)^2;
       xk(3,3)=12;
       xk(3,4)=-6*xl(ii);
       xk(4,4)=4*xl(ii)^2;
       for i=1:4
          for j=1:4
             xm(j,i)=xm(i,j);
             xk(j,i)=xk(i,j);
          end
       end
       for i=1:4
          for j=1:4
             xm(i,j)=(((xmas(ii)*xl(ii))/420))*xm(i,j);
             xk(i,j)=((e*xi(ii))/(xl(ii)^3))*xk(i,j);
          end
       end
       for i=1:n
          for j=1:4
             ait(i,j)=ai(j,i);
          end
       end
       xka=xk*ai;
       xma=xm*ai;
       aka=ait*xka;
       ama=ait*xma;
      for i=1:n
          for j=1:n
             bigm(i,j)=bigm(i,j)+ama(i,j);
             bigk(i,j)=bigk(i,j)+aka(i,j);
          end
       end
end    
%   Application of boundary conditions
%   Rows and columns corresponding to zero displacements are deleted
    
     bigk(1:2,:) = [];
     bigk(:,1:2) = [];
        
     bigm(1:2,:) = [];
     bigm(:,1:2) = [];
   


%   Calculation of eigenvector and eigenvalue
    
    [L, V] = eig (bigk,bigm)
    
%   Natural frequency
    V1 = V.^(1/2)
    W = diag(V1)
     f=W/(2*pi)   


