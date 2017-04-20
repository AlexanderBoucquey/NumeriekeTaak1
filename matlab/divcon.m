% function [lam,U] =divcon(A)
%
% Verdeel- en heersmethode om de eigenwaarden en eigenvectoren van een
% symmetrische, tridiagonale matrix te berekenen
%
% Invoer:
% A             een symmetrisch, tridiagonale matrix
%
% Uitvoer:
% lam           eigenwaarden
% U             eigenvectoren

function [lam,U] = divcon(A)

  % Initialize
    n= size(A,1);

    if n==0
      lam = []; U=[];
      return
    end
    
    
    if n==1
      U=1.0;
      lam=A(1,1);
      return
    end

    if n==2
      [U,lam]=eig(A);
      lam=diag(lam);
      return
    end
    if n>2
      s=floor(n/2);
      U = zeros(n);
      % Start the recursive procedure 
      beta = A(s+1,s);
      T1 = A(1:s,1:s); T1(s,s) = T1(s,s)-beta;
      T2 = A(s+1:n,s+1:n); T2(1,1) = T2(1,1)-beta;
      [lam1,U1]=divcon(T1);
      [lam2,U2]=divcon(T2);
      U(1:s,1:s)=U1;
      U(s+1:n,s+1:n)=U2;
      z = [U1(end,:) U2(1,:)]';
      if(abs(beta)<1e-15)
          disp('small beta');
          lam = [lam1;lam2];
          return;
      end
      d = [lam1;lam2];
     
      [l1,Q1]=eig_r1diag(d,z,beta);
      lam =l1;
      U = U*Q1;
    
    end

function [lam,U] = eig_r1diag(a,b,ro)
% a,b columns!! ro number
%
% [lam U] = eig_r1diag(a,b,ro):
%
% The elements of lam are the EIGENVALUES and the columns of U are the
% orthonormal EIGENVECTORS, with last elements nonnegative, of the real
% symmetric matrix A = [diag(a) +ro*b* b'], so a diagonal matrix plus a rank-one modification.
%
% This code is an addapted version of the program eva written by Gragg which starts its iterations at the (iterates
% from the) poles closest to the zeros [3]. 

     n = length(a);   U = eye(n); V=[];

   if ro ==0
        lam=a; return 
   end

   if n < 2
        lam = a(n)+ro*b(n)^2;  return
   end

if ro < 0
     a=-a;
     rho=-ro;
else rho=ro;
end

%    Order A = [diag(a) b; b' c] and scale it so b is nonnegative.
     [a p] = sort(-a);   a = - a;   b = b(p);   U = U(:,[p; n]);

     p = find(b < 0);
     if length(p) > 0
        b(p) = - b(p);   U(:,p) = - U(:,p);
     end

%    b-deflation.
     m = 3*(n + 6);   tolb = m*eps*norm(b);

     p = find(b <= tolb);   l = length(p);
     if l > 0
        b(p) = zeros(l,1);
     end

%    Combo-deflation
    p = find(b > 0);   l = length(p);   k = 1;
     while k < l
        i = p(k);      j = p(k+1);   q = [i j]';
        a0 = a(i);     a1 = a(j);    b0 = b(i);     b1 = b(j);
        s = [a0 a1];   s = abs(s);   t = [b0 b1];   t = max([s t]);
        if b0 < b1
           r = b0/b1;
        else
           r = b1/b0;
        end
        d = a0 - a1;   s = d/(r + 1/r);   tolc = 2*eps*t;
        if s <= tolc
           s = r*r;   t = 1 + s;   s = s*d/t;   d = sqrt(t);
           if b0 < b1
              a(q) = [a0-s; a1+s];   b(q) = [0; b1*d];
              cs = 1/d;              sn = r/d;
           else
              a(q) = [a1+s; a0-s];   b(q) = [0; b0*d];
              cs = r/d;              sn = 1/d;
           end
           U(:,q) = U(:,q)*[-cs sn; sn cs];   p = p([1:k-1 k+1:l]');
           l = l - 1;                         k = max(1,k-1);
        else
           k = k + 1;
        end
     end

     lam = [a];

%    Solve the eigenproblem for the ordered and reduced arrow matrix, again
%    called A.  The index k refers to the kth eigenpair lam(p(k)),
%    U(:,p(k))
     a = a(p);   b = b(p);   p = [p];   n = length(p);   m = 2*(n + 6);

     if n < 2
        lam(p) = a(n)+rho*b(n)^2;   return
     end

     for k = 1:n

%       The index vector q omits reference to the neighboring finite poles.
        if k > 1  &  k < n
           q = [1:k-2 k+1:n];
        else
           if k == 1
              q = 2:n;
           else
              q = 1:n-2;
           end
        end

%       Determine a closest pole, shift, to lam(p(k)).  Shift A.  Compute
%       the initial approximation x to lam(p(k)) - shift.    
        if k == 1
           shift = a(1);          as = a - shift;             cs = -1/rho;
           t = b(q)./as(q);       u = sqrt(sum(t.*t));   t = b(q).*t;
           s = sum([-t; cs])/2;   t = u*b(1);
           if s > 0
              if s < t
                 r = s/t;   x = b(1)*(r + sqrt(1 + r*r))/u;
              else
                 r = t/s;   x = b(1)*(1 + sqrt(1 + r*r))/r/u;
              end
           else
              if abs(s) < t
                 r = s/t;   x = - b(1)/(r - sqrt(1 + r*r))/u;
              else
                 r = t/s;   x = - b(1)*r/(1 + sqrt(1 + r*r))/u;
              end
           end
%x1=x+shift
        end

        if 1 < k  &  k <= n
           x = sum(a(k-1:k))/2;   t = b.*(b./(a - x));   f = sum([t;1/rho]);
           if f > 0
              shift = a(k);     as = a - shift;   pole = as(k-1);
           else
              shift = a(k-1);   as = a - shift;   pole = as(k);
           end
           t = b(q)./as(q);     dg = sum(t.*t);
           t = b(q).*t;       g = sum([t; 1/rho]);   alpha = pole*dg - g;
           s = b(k-1:k).^2;   t = pole*g;          beta = sum([s; t]);
           s = flipud(s);     t = s.*as(k-1:k);    f = - sum(t);
           beta = beta/2;     d = alpha*f;         d = sqrt(beta^2 - d);
           if beta > 0
              x = - f/(beta + d);
           else
              x = (d - beta)/alpha;
           end

        end



%       Iterate to compute lam(p(k)).
        for j = 0:20

%          Compute g, f and the termination tolerance tolt.
           v = as - x;   s = b./v;   t = b.*s;   g = sum([t(q); 1/rho]);

           if k > 1  &  k <= n
              f = g + sum(t(k-1:k));
           else
              %if k == 1
                 f = g + t(1);
              %else
              %   f = g + t(n-1);
              %end
           end


          
%          Compute dg := g' and ddg := g''/2.
           t = s.*s;   dg = sum(t(q));   t = t./v;   ddg = sum(t(q));

%          Compute alpha and beta.
           if k > 1  &  k <= n
              s = sum(b(k-1:k).^2);         d = pole - x;
              t = d - x;                    u = x*d;
              alpha = ddg + (t*dg - g)/u;   beta = dg + (t*g + s)/u;
           else
              alpha = ddg + dg/x;           beta = dg + g/x;
           end


%          Compute the new iterate.
           beta = beta/2;   d = sqrt(beta^2 - alpha*f);
           if beta > 0
              x = x - f/(beta + d);
           else
              x = x + (d - beta)/alpha;
           end

        end

%       Accept the eigenvalue.  Save a - lam(p(k)) = as - x for later use in
%       computing the eigenvectors.
         lam(p(k)) = shift + x;
	 V = [V v];
        
     end

%    Compute the modified vector b.
     e = ones(n,1);   A = e*a' - a*e' + eye(n);   A = [A];
     A = V'./A;         b = sqrt(-(prod(A)))';

%    Compute numerically orthonormal eigenvectors of the modified reduced
%    matrix.  Compute U.
    A =  b*e';           V = - A./V;  
    A = sort(V.*V);      t = sqrt(sum(A))';   A = e*t';     V = V./A;
    U(:,p) = U(:,p)*V;

%    Sort (actually merge) and quit.
     [lam p] = sort(-lam);    lam = - lam;
      if ro < 0  
         lam = - lam; 
      end   
      U = U(:,p);



    
    
    
