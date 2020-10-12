function ktAcq = uncompressKT( ktAcq, ktSmp )

% Uncompress k-t SENSE ktAcq to restore size according to original sampling 
% pattern
%
% - Based on code from Lucilio
%
% Input:
% - ktAcq - expect in Josh order: [nX, nY, nF, nC, nZ]
% - ktSmp - original k-t sampling pattern in shape: [1, nY, nF, 1, 1]
%
% see also: compressKT
%
% Tom Roberts (t.roberts@kcl.ac.uk)


%% Rearrange data order
% Josh order:    [nX, nY, nF, nC, nZ]
% Lucilio order: [nX, nY, nZ, nC, nF] 

ktSmp = permute(ktSmp, [1,2,5,4,3]);
ktSmp = ktSmp(:,:,1,:,:); % equivalent to A in Lucilio code

ktAcq = permute(ktAcq, [1,2,5,4,3]);

%% Uncompress K-Space

A = ktSmp;
y = ktAcq;


NY = size( y );
% NY( end+1:5 ) = 1;
NA = size( A );
% NA( end+1:5 ) = 1;
NR = NY./NA;
NR(2) = 1;
A = repmat( A, NR ); % repeat ktSmp
z = A; z(:) = 0; % seems to just be z = zeros(size(A)); ???
z( A==1 ) = y;
y = z;


%% Outputs
% restore Josh order
ktAcq = permute(y, [1,2,5,4,3]);


% compressKT(...)
end