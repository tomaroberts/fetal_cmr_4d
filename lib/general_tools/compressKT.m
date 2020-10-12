function [ktAcq, ktSmp] = compressKT( ktAcq )

% Compress k-t SENSE ktAcq to remove zeros and reduce file size
%
% - Based on code from Lucilio
% - Useful for compressing *_kspace.mat data
%
% Input:
% - ktAcq - expect in Josh order: [nX, nY, nF, nC, nZ]
%
% see also: uncompressKT
%
% Tom Roberts (t.roberts@kcl.ac.uk)


%% k-t Sampling pattern
dimC  = 4;
dimX  = 1;
ktSmp = single( sum( sum( ktAcq, dimC ), dimX ) ~= 0 );


%% Rearrange data order
% Josh order:    [nX, nY, nF, nC, nZ]
% Lucilio order: [nX, nY, nZ, nC, nF] 

ktSmp = permute(ktSmp, [1,2,5,4,3]);
ktSmp = ktSmp(:,:,1,:,:); % equivalent to A in Lucilio code

ktAcq = permute(ktAcq, [1,2,5,4,3]);

%% Compress K-Space

A = ktSmp;
y = ktAcq;

NY = size(y);
% NY(end+1:5)=1;
NA = size(A);
% NA(end+1:5)=1;
NR = NY./NA;     
y = y(repmat(A,NR)==1);
NY(2) = sum(dynInd(A,1,5));
y = reshape(y,NY);


%% Outputs
ktSmp = permute(A, [1,2,5,4,3]);
ktAcq = permute(y, [1,2,5,4,3]); % output in Josh order


% compressKT(...)
end


%% Subfunctions
function x=dynInd(x,ind,dim,y)

%DYNIND   Performs dynamic indexing over multidimensional arrays
%
% Lucilio code.

Ndim=length(dim);
%if iscell(ind);assert(length(ind)==Ndim,'Lenght of cell of indexes (%d) has to match length of dimensions to index (%d)',length(ind),Ndim);else assert(~(Ndim>1 && length(ind)~=Ndim),'Length of vector of indexes (%d) has to match the length of dimensions to index (%d). You may want to consider operating with cells of indexes for vectorial indexing',length(ind),Ndim);end%Takes time and may simply work without error control
%assert(length(unique(dim))==length(dim),'Dimensions to index do not allow repetitions');%Takes time and may simply work without error control

ndx=max(ndims(x),max(dim));

if ndx==1;subs={':'};%Calling repmat takes time
elseif ndx==2;subs=[{':'} {':'}];
elseif ndx==3;subs=[{':'} {':'} {':'}];
elseif ndx==4;subs=[{':'} {':'} {':'} {':'}];
elseif ndx==5;subs=[{':'} {':'} {':'} {':'} {':'}];
elseif ndx==6;subs=[{':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==7;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==8;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==9;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==10;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==11;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==12;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==13;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
elseif ndx==14;subs=[{':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'} {':'}];
else subs=repmat({':'}, [1 ndx]);
end

%Note error control is not implemented, if the user tries to index outside 
%the size of x along a given dimension, matlab should throw an error
if ~iscell(ind)
    if Ndim~=1;subs(dim)=num2cell(ind);else subs{dim}=ind;end
else
    subs(dim)=ind;%If this fails it may be because the length of the cell does not match the length of dimensions
end
if nargin<4;x=x(subs{:});else x(subs{:})=y;end

% end dynIn(...)
end