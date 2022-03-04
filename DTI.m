function [ADCv,FAv,EigenVectors,EigenValues,DiffusionTensor]=DTI(DTIdata)
% This function will perform DTI calculations on a certain
% DTI dataset.

% ADC: A 3D matrix with the Apparent Diffuse Coefficient (ADC)
% FA: A 3D matrix with the fractional anistropy (FA)
% VectorF: A 4D matrix with the eigenvector matrix in each pixel
% Y:the  A 4D matrix with eigen values of every voxel
% DifT: A 3D matrix with all  Diffusion tensors [Dxx,Dxy,Dxz,Dyx,Dyy,Dyz,Dzx,Dzy,Dzz]

S=zeros([1 1 1 length(DTIdata)],'single');
% Make a 3D matrix to store the zero gradient voxel volume(s)
S0=1;% Make a matrix to store the gradients
H=zeros(length(DTIdata),3,'single');
% Make a vector to store the different B-values (timing)
Bvalue=zeros(length(DTIdata),1,'single');
% Read the input data (DTIdata), and seperate the zero gradient
% and other gradients into different matrices.
 voxelg=0;
for i=1:length(DTIdata)
    
        voxelg=voxelg+1;
        S(:,:,:,voxelg)=single(DTIdata(i).VoxelData);
        H(voxelg,:)=single(DTIdata(i).Gradient);
        Bvalue(voxelg) = single(DTIdata(i).Bvalue);
    
end
%create b matrix
b=zeros([3 3 size(H,1)]);
for i=1:size(H,1)
    b(:,:,i)=Bvalue(i)*H(i,:)'*H(i,:);
end

% Convert measurement intentensity into absorption (log)
Slog=zeros(size(S),'single');
for i=1:size(H,1)
    Slog(:,:,:,i)=log((S(:,:,:,i)./S0)+eps);
end

Bv=squeeze([b(1,1,:),2*b(1,2,:),2*b(1,3,:),b(2,2,:),2*b(2,3,:),b(3,3,:)])';
% Create a matrix to store the Diffusion tensor of every voxel
% [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz]
% DifT=zeros([1 1 1 9],'single');
% Create a matrix to store the eigen values of every voxel
% Y=zeros([1 1 1 3],'single');
% Create a matrix to store the fractional anistropy (FA)
% FA=zeros([1 1 1],'single');
% Create a matrix to store the Apparent Diffuse Coefficient (ADC)
% ADC=zeros([1 1 1],'single');
%create A 4D matrix with the eigenvector matrix in each pixel
% VectorF=zeros([1 1 1 9],'single');
% Loop through all voxel coordinates
for x=1:size(S,1)
    for y=1:size(S,2)
        for z=1:size(S,3)
            
            % Only process a pixel if it isn't background
           
                
                % Calculate the Diffusion tensor [Dxx,Dxy,Dxz,Dyy,Dyz,Dzz],
                % with a simple matrix inverse.
                Z=-squeeze(Slog(x,y,z,:));
                M=Bv\Z;
               
                % The DiffusionTensor (Remember it is a symetric matrix,
                % thus for instance Dxy == Dyx)
                DiffusionTensor=[M(1) M(2) M(3); M(2) M(4) M(5); M(3) M(5) M(6)];
                % Calculate the eigenvalues and vectors, and sort the 
                % eigenvalues from small to large
                [EigenVectors,D]=eig(DiffusionTensor); EigenValues=diag(D);
              
                
                % Regulating of the eigen values (negative eigenvalues are
                % due to noise and other non-idealities of MRI)
                if((EigenValues(1)<0)&&(EigenValues(2)<0)&&(EigenValues(3)<0)), EigenValues=abs(EigenValues);end
                if(EigenValues(1)<=0), EigenValues(1)=eps; end
                if(EigenValues(2)<=0), EigenValues(2)=eps; end
                
                % Apparent Diffuse Coefficient
                ADCv=(EigenValues(1)+EigenValues(2)+EigenValues(3))/3;
                
                % Fractional Anistropy (2 different definitions exist)
                % First FA definition:
                %FAv=(1/sqrt(2))*( sqrt((EigenValues(1)-EigenValues(2)).^2+(EigenValues(2)-EigenValues(3)).^2+(EigenValues(1)-EigenValues(3)).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
                % Second FA definition:
                FAv=sqrt(1.5)*( sqrt((EigenValues(1)-ADCv).^2+(EigenValues(2)-ADCv).^2+(EigenValues(3)-ADCv).^2)./sqrt(EigenValues(1).^2+EigenValues(2).^2+EigenValues(3).^2) );
                
                % Store the results of this pixel in the volume matrices
%                 ADC(x,y,z)=ADCv;
%                 Y(x,y,z,:)=EigenValues;
%                 DifT(x,y,z,:)=reshape(DiffusionTensor,1,9);
%                 
%                 FA(x,y,z)=FAv;
%                 VectorF(x,y,z,:)=reshape(EigenVectors,1,9);

            
        end
    end
end
