%%                                                            Developed by Kife Intasar Bin Iqbal (BUET)                                         %%

function P = customExpm(H,contol,nonZtol,varargin)


% Convergence Point Value
if nargin>3
    convergenceCriteria=varargin{1}; 
else
    convergenceCriteria=contol;
end

% NonZero Tolerance Value
if nargin>4
    nonZeroTol=varargin{2}; 
else
    nonZeroTol=nonZtol;
end

mat_norm=norm(H,1);
 delta=1;
n_squarings=max([0 ceil(log2(mat_norm))]); 
scaling_factor=2^n_squarings;
H=H*scaling_factor^-1;
H=nonZeroTol*round((1/nonZeroTol)*H);

%% Run Taylor series 
P=speye(size(H)); 
nextTerm=P; 
n=1; 

% Sparsity ans size check
if nnz(H)/numel(H)>0.25 && numel(H)^0.5<=64
    H=full(H);
else
    H=sparse(H);
end

while delta>convergenceCriteria
    % Compute the next term
    if issparse(nextTerm)
        nextTerm=(1/n)*H*nextTerm; % order matters
        nextTerm=nonZeroTol*round((1/nonZeroTol)*nextTerm); % Eliminate small quantity
        if nnz(nextTerm)/numel(nextTerm)>0.25
            nextTerm=full(nextTerm);
        end
    else
        nextTerm=(1/n)*nextTerm*H;
    end
    delta=norm(nextTerm,1); % check residual norm
    P=P+nextTerm; 
    n=n+1; %Add to the total and increment the counter
end
 
H=[];
nextTerm=[];

%% Squaring
P=nonZeroTol*round((1/nonZeroTol)*P);
for n=1:n_squarings
    P=P*P;
    if issparse(P)
        if nnz(P)/numel(P)<0.25
            P=nonZeroTol*round((1/nonZeroTol)*P); % Eliminate small quantity
        else
            P=full(P);
        end
    end
end
end
