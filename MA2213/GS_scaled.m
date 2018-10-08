function x       =  GS_scaled(A,b,digit,type)
    %%  function to perform the Gauss elimination with scaled paritial
%        pivoting and rounding or chopping arithmatics
%   input: Ax    = b
%          digit, how many digits to keep when using rounding or chopping
%          type, be 'round' or 'chop', which is a string
if ~exist('digit','var')
    digit        = 4;
end
if ~exist('type','var')
    type         = 'round';    % using rounding arithmatic
end
n                = length(A);  % dimension of A 
global d t
d                = digit;
t                = type;
%% perform GS elimination 
for i            = 1:n-1
    %% partial pivoting
                               % largest value of A(i,j);
                               % largest value of each row
    v            = zeros(n-i+1,1);
    for j        =i:n
       v(j-i+1)  = max(abs(A(j,i:n)));
    end
    for j        =i:n
        v(j-i+1) = abs(A(j,i))/v(j-i+1);
    end
    fprintf('The relative pivoting is:\n');
    disp(v');
    [~,idx]      = max(v);
    idx          = idx + i-1;
    
    fprintf('  change %d-row and %d-th row\n',i,idx);
    v            = A(i,:);
    A(i,:)       = A(idx,:);
    A(idx,:)     = v;
    fprintf('New A after row exchange\n');
    disp(A);
                               % exchange b also
    v            = b(i);
    b(i)         = b(idx);
    b(idx)       = v;
    fprintf('New b after row exchange\n');
    disp(b);
   for j         = i+1:n
      mji        = tdiv(A(j,i),A(i,i));
      fprintf('m(%d%d) is %g\n',j,i,mji);
                               % -mji row i + row j 
      
      for k      = i+1:n
         A(j,k)  = tsub(A(j,k),   tmulti(mji,A(i,k)));
         
      end
      b(j)       = tsub(b(j),tmulti(mji,b(i)));
   end
   A(i+1:n,i)    =0;
   fprintf('The new A is:\n');
   disp(A);
   fprintf('The new b is:\n');
   disp(b);
end
%% call backward substitution
x                = backsub(A,b);
end


function y       = tadd(x,y)   % 
x                = tround(x);
y                = tround(y);
y                = tround(x+y);
end

function y       = tmulti(x,y)
x                = tround(x);
y                = tround(y);
y                = tround(x*y);
end
function y       = tdiv(x,y)
x                = tround(x);
y                = tround(y);
y                = tround(x/y);
end

function y       = tsub(x,y)
x                = tround(x);
y                = tround(y);
y                = tround(x-y);
end

function y       = tround(x)
global d t 
switch t
    case 'round'
        y        = round(x,d,'significant');
    otherwise
        ex       = ceil(log10(abs(x)));
        x        = x/(10^ex);
        y        = fix(x*(10^d))/(10^d);
        x        = x*(10^ex);
end

end

function x       = backsub(A,b)
%% A is a upper triangular matrix
n                = length(A);
x                = zeros(n,1);
for i            =n:-1:1
    x(i)         = tround(b(i)/A(i,i));
    fprintf('x(%d) is %g ',i,x(i));
    for j        = 1:i-1
       b(j)      = tround(b(j) -  tmulti(A(j,i),x(i)));
    end
    fprintf('new b is %g\n',b);
    
end
end