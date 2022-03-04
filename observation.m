clc
clear
%%

%creat signal in 6/15/21 direction
% u6=[1/sqrt(2) 0 1/sqrt(2);
%     -1/sqrt(2) 0 1/sqrt(2);
%     0 1/sqrt(2) 1/sqrt(2);
%     0 1/sqrt(2) -1/sqrt(2);
%     1/sqrt(2) 1/sqrt(2) 0;
%     -1/sqrt(2) 1/sqrt(2) 0];
s=2*cos(pi/5);
A=1/sqrt(1+s^2);
B=s/sqrt(1+s^2);
C=(s-1)/2;
D=s/2;
u6=[0 A -B;0 -A -B;A -B 0;A B 0;-B 0 A; -B 0 -A];
u15=[0 0 1;0 1 0;C -D -0.5;C D -0.5;-C -D -0.5;C -D 0.5;-0.5 C D;0.5 -C D;0.5 C -D;0.5 C D;-D 0.5 C;D -0.5 C;D 0.5 -C;-D -0.5 -C;-1 0 0];
u21=[0 0 1;0 A -B;0 -A -B;0 1 0;C -D -0.5;C D -0.5;-C -D -0.5;C -D 0.5;-0.5 C D;0.5 -C D;0.5 C -D;0.5 C D;A -B 0;A B 0;-D 0.5 C;D -0.5 C;D 0.5 -C;-D -0.5 -C;-B 0 A; -B 0 -A;-1 0 0];
% [u162,fcs] = icosahedron(2);
%%
B1=[100 300 500 1000 2000 3000 4000];%s/mm^2
thegma1=[2 5 10 20 30];%%ms

Delta=[20 40 60 80 100];%ms
for t1=1:length(B1)
    for t2=1:length(thegma1)
        for t3=1:length(Delta)
qq1(t1,t2,t3)=sqrt(B1(t1)*10^3./(Delta(t3)-thegma1(t2)/3));
q1(t1,t2,t3)=qq1(t1,t2,t3)/(2*pi);
        end
    end
end

%%
 MD=cell(length(B1),length(thegma1),length(Delta));
FA=cell(length(B1),length(thegma1),length(Delta));
sim_azi=cell(length(B1),length(thegma1),length(Delta));
sim_ele=cell(length(B1),length(thegma1),length(Delta));
EigenValues1=cell(length(B1),length(thegma1),length(Delta));

EigenValues2=cell(length(B1),length(thegma1),length(Delta));
EigenValues3=cell(length(B1),length(thegma1),length(Delta));
Dif=cell(length(B1),length(thegma1),length(Delta));
sim_vector1=cell(length(B1),length(thegma1),length(Delta));
u=u21;
n_molecule=50000;
n_step=131;%%
displacement=zeros(1,n_molecule);
signal=zeros(1,length(u));
step_location=cell(1,n_molecule);


for k=2:5
for kkk=1:35
    a=k
    aa=kkk
fileID=fopen(['D:\doctorat\programme\C\inmyocyte_outmyocyte_permeability_DiaOrienHeter\inmyocyte_outmyocyte_permeability_DiaOrienHeter\startpoint_mean_r=0.01_std_r=0.001_step=1_DM=50_0.3_0.3_0.3_azi=0_ele=-85-5-85\',num2str(k),'_',num2str(kkk),'_startpoint_mean_r=0.01_std_r=0.001_step=1_DM=50_0.3_0.3_0.3_per=0.002_molecules=50000_DELTA=0.001_tau=0.1_D=1_2.5_cycle=130_newnew.txt']);
% fileID=fopen('D:\doctorat\programme\C\inmyocyte_outmyocyte_permeability_heterogeneity2\inmyocyte_outmyocyte_permeability_heterogeneity2\r_mean=10_std=1_varied angle_0.7_0.7_0.7_new\8_r_mean=0.01_r_std=0.001_DM=50_0.7_0.7_0.7_per=0.002_molecules=50000_DELTA=0.001_tau=0.1_D=1_2.5_cycle=120_newnew.txt');
 aaa=textscan(fileID,'%f');
  fclose(fileID);  %%


for i=1:n_molecule
    for j=1:n_step

       step_location{i}(j,:)=aaa{1}(3*n_step*(i-1)+3*j-2:3*n_step*(i-1)+3*j);
    end
end
for t1=1:length(B1)
for t2=1:length(thegma1)
for t3=1:length(Delta)
   

for j=1:length(u)
        for i=1:n_molecule

        displacement(i)=dot((sum(step_location{i}(Delta(t3)+1:Delta(t3)+thegma1(t2),:))-sum(step_location{i}(1:thegma1(t2),:)))/(thegma1(t2)),u(j,:));
%         displacement(i)=dot(step_location{i}(t+1,:)-step_location{i}(1,:),u(j,:));
        end
        signal(j)=(1/n_molecule)*sum(cos(qq1(t1,t2,t3)*displacement));
% signal(t,j)=(1/n_molecule)*sum(cos(qq1*displacement));
end


DTIdata=struct();
n_voxel=1;



for j=1:length(u)

    DTIdata(j).VoxelData = abs(signal(j));
    DTIdata(j).Gradient = u(j,:);
    DTIdata(j).Bvalue=B1(t1);%%%(s/mm^2)
%   DTIdata(j).Bvalue=B1(t);%%%(s/mm^2)
end

 [MDv,FAv,VectorF,Y,DifT]=DTI(DTIdata);
 Dif{t1,t2,t3}(k,kkk,:,:)=DifT;
MD{t1,t2,t3}(k,kkk)=MDv;
FA{t1,t2,t3}(k,kkk)=FAv;
[EigenVectors,D]=eig(DifT);
EigenValues=diag(D);
EigenValues1{t1,t2,t3}(k,kkk)=EigenValues(1);
EigenValues2{t1,t2,t3}(k,kkk)=EigenValues(2);
EigenValues3{t1,t2,t3}(k,kkk)=EigenValues(3);
VectorF1=reshape(VectorF,3,3);
Fiber_dir=VectorF1(:,end);%%fiber direction
sim_vector=squeeze(Fiber_dir)';
sim_vector1{t1,t2,t3}(k,kkk,:,:)=sim_vector;
%ele~[-1/2*pi,1/2*pi]
sim_ele{t1,t2,t3}(k,kkk)=asin(sim_vector(3))*180/pi;

%azi~[0,181]
if sim_vector(1)>=0
  sim_azi{t1,t2,t3}(k,kkk)=atan(sim_vector(2)/sim_vector(1))*180/pi;
else
  sim_azi{t1,t2,t3}(k,kkk)=(atan(sim_vector(2)/sim_vector(1))+pi)*180/pi;
end
if  sim_azi{t1,t2,t3}(k,kkk)>90
    sim_azi{t1,t2,t3}(k,kkk)=sim_azi{t1,t2,t3}(k,kkk)-180;
end
% if  sim_azi{t1,t2,t3}(k,kkk)<0
%     sim_azi{t1,t2,t3}(k,kkk)=sim_azi{t1,t2,t3}(k,kkk)+180;
% end
%  if sim_vector(max_location)*real_vector(max_location)<0
%     sim_vector1=-sim_vector;
% else
%     sim_vector1=sim_vector;
% end
% %ele~[-1/2*pi,1/2*pi]
% sim_ele(k)=asin(sim_vector1(3))*180/pi;
%
% %azi~[0,181]
% if sim_vector(1)>=0
%   sim_azi(k)=atan(sim_vector1(2)/sim_vector1(1))*180/pi;
% else

%   sim_azi(k)=(atan(sim_vector1(2)/sim_vector1(1))+pi)*180/pi;
% end
end
end
end
end
end