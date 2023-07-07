function [Fs,Fs2,Fs3,Fs4,Fs5,Fs6,Fs7,Fs8,Fs9,Fs10,Fs11,Fs12,Fs13,Fs14,Fs15,Fs16,Fs17,Fs18,Fs19,Fs20]...
    = Main_RF(z,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z20,ind)

%% 
% if size(z,2) ~= 1
%     z=z.';
% end
% function [Fs,Fs2,Fs3,Fs4,Fs5,Fs6,Fs7,Fs8,Fs9,Fs10,Fs11,Fs12,Fs13,Fs14,Fs15,Fs16,Fs17,Fs18,Fs19,Fs20]...
%     = Main_RF(z,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z20,zz,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10,...
%      zz11,zz12,zz13,zz14,zz15,zz16,zz17,zz18,zz19,zz20,ind)
%% c-phi correlation
% cor=[1,0;0,1];
% L = chol(cor,'lower');
% 
% for i=1:9
%     A = [z(i) zz(i)]' ;        A =  L * A;   z(i) = A(1,1);     zz(i) = A(2,1);
%     A = [z2(i) zz2(i)]' ;      A =  L * A;  z2(i) = A(1,1);     zz2(i) = A(2,1);
%     A = [z3(i) zz3(i)]' ;      A =  L * A;  z3(i) = A(1,1);     zz3(i) = A(2,1);
%     A = [z4(i) zz4(i)]' ;      A =  L * A;  z4(i) = A(1,1);     zz4(i) = A(2,1);
%     A = [z5(i) zz5(i)]' ;      A =  L * A;  z5(i) = A(1,1);     zz5(i) = A(2,1);
%     A = [z6(i) zz6(i)]' ;	   A =  L * A;	z6(i) = A(1,1);	    zz6(i) = A(2,1);
%     A = [z7(i) zz7(i)]' ;	   A =  L * A;	z7(i) = A(1,1);	    zz7(i) = A(2,1);
%     A = [z8(i) zz8(i)]' ;	   A =  L * A;	z8(i) = A(1,1);	    zz8(i) = A(2,1);
%     A = [z9(i) zz9(i)]' ;	   A =  L * A;	z9(i) = A(1,1);	    zz9(i) = A(2,1);
%     A = [z10(i) zz10(i)]' ;	   A =  L * A;	z10(i) = A(1,1);	zz10(i) = A(2,1);
%     A = [z11(i) zz11(i)]' ;	   A =  L * A;	z11(i) = A(1,1);	zz11(i) = A(2,1);
%     A = [z12(i) zz12(i)]' ;	   A =  L * A;	z12(i) = A(1,1);	zz12(i) = A(2,1);
%     A = [z13(i) zz13(i)]' ;	   A =  L * A;	z13(i) = A(1,1);	zz13(i) = A(2,1);
%     A = [z14(i) zz14(i)]' ;	   A =  L * A;	z14(i) = A(1,1);	zz14(i) = A(2,1);
%     A = [z15(i) zz15(i)]' ;	   A =  L * A;	z15(i) = A(1,1);	zz15(i) = A(2,1);
%     A = [z16(i) zz16(i)]' ;	   A =  L * A;	z16(i) = A(1,1);	zz16(i) = A(2,1);
%     A = [z17(i) zz17(i)]' ;	   A =  L * A;	z17(i) = A(1,1);	zz17(i) = A(2,1);
%     A = [z18(i) zz18(i)]' ;	   A =  L * A;	z18(i) = A(1,1);	zz18(i) = A(2,1);
%     A = [z19(i) zz19(i)]' ;	   A =  L * A;	z19(i) = A(1,1);	zz19(i) = A(2,1);
%     A = [z20(i) zz20(i)]' ;	   A =  L * A;	z20(i) = A(1,1);	zz20(i) = A(2,1);
% end
%%    cohesion Interval Field 
load matlab.mat  %vector that contains coordinates of centroids of elements
Lh       = 20;   %correlation length of interval field - [m]
Lv       = 2;   %correlation length of interval field - [m]
muRF    = 23; %mean value 
sigmaRF = 3.45;  %radius
% corrRF  = @(d1,d2) exp(-d1.^2/Lh^2-d2.^2/Lv^2);   %correlation function
% corrRF  = @(d1,d2)  exp(-d1/Lh-d2/Lv) .* (1+d1/Lh) .* (1+d2/Lv);   %correlation function
corrRF  = @(d1,d2)  exp(-d1/Lh-d2/Lv);
Perc    = 0.99; %percentage  
RF  = CreateRF(...
    'Permeability of Soil',...     %name of interval field
    muRF,...     %mean value 
    sigmaRF,...    %radius
    @(d1,d2)corrRF(d1,d2),...  %correlation function
    Centroids,...   %centroids of finite elements
    Perc...         %percentage of represented variability
    );
RFc = RF.Std2Phy(z);
RFc2 = RF.Std2Phy(z2);
RFc3 = RF.Std2Phy(z3);
RFc4 = RF.Std2Phy(z4);
RFc5 = RF.Std2Phy(z5);
RFc6 = RF.Std2Phy(z6);
RFc7 = RF.Std2Phy(z7);
RFc8 = RF.Std2Phy(z8);
RFc9 = RF.Std2Phy(z9);
RFc10 = RF.Std2Phy(z10);
RFc11 = RF.Std2Phy(z11);
RFc12 = RF.Std2Phy(z12);
RFc13 = RF.Std2Phy(z13);
RFc14 = RF.Std2Phy(z14);
RFc15 = RF.Std2Phy(z15);
RFc16 = RF.Std2Phy(z16);
RFc17 = RF.Std2Phy(z17);
RFc18 = RF.Std2Phy(z18);
RFc19 = RF.Std2Phy(z19);
RFc20 = RF.Std2Phy(z20);
%%   friction angle Interval Field (489)
% Lh       = 25;   %correlation length of interval field - [m]
% Lv       = 2.5;   %correlation length of interval field - [m]
% muRF    = 20; %mean value 
% sigmaRF = 4;  %radius
% % corrRF  = @(d1,d2) exp(-d1.^2/Lh^2-d2.^2/Lv^2);   %correlation function
% corrRF  = @(d1,d2) exp(-d1.^2/Lh^2-d2.^2/Lv^2);   %correlation function
% Perc    = 0.95; %percentage  
% RF  = CreateRF(...
%     'Permeability of Soil',...     %name of interval field
%     muRF,...     %mean value 
%     sigmaRF,...    %radius
%     @(d1,d2)corrRF(d1,d2),...  %correlation function
%     Centroids,...   %centroids of finite elements
%     Perc...         %percentage of represented variability
%     );
% RFphi = RF.Std2Phy(zz);
% RFphi2 = RF.Std2Phy(zz2);
% RFphi3 = RF.Std2Phy(zz3);
% RFphi4 = RF.Std2Phy(zz4);
% RFphi5 = RF.Std2Phy(zz5);
% RFphi6 = RF.Std2Phy(zz6);
% RFphi7 = RF.Std2Phy(zz7);
% RFphi8 = RF.Std2Phy(zz8);
% RFphi9 = RF.Std2Phy(zz9);
% RFphi10 = RF.Std2Phy(zz10);
% RFphi11 = RF.Std2Phy(zz11);
% RFphi12 = RF.Std2Phy(zz12);
% RFphi13 = RF.Std2Phy(zz13);
% RFphi14 = RF.Std2Phy(zz14);
% RFphi15 = RF.Std2Phy(zz15);
% RFphi16 = RF.Std2Phy(zz16);
% RFphi17 = RF.Std2Phy(zz17);
% RFphi18 = RF.Std2Phy(zz18);
% RFphi19 = RF.Std2Phy(zz19);
% RFphi20 = RF.Std2Phy(zz20);

%%  input
Pref.Str2Num = 'never';
x_doc=xml_read('D:\Cheng\FS\initial\Slope0.xml',Pref);
%%  edit 
for i=1:1:1033
    x_doc.Functions.Func3Ds.Fn3D(1).Points.Points_(i).ATTRIBUTE.Z= RFc(i);
    x_doc.Functions.Func3Ds.Fn3D(3).Points.Points_(i).ATTRIBUTE.Z= RFc2(i);  
    x_doc.Functions.Func3Ds.Fn3D(5).Points.Points_(i).ATTRIBUTE.Z= RFc3(i);
    x_doc.Functions.Func3Ds.Fn3D(7).Points.Points_(i).ATTRIBUTE.Z= RFc4(i);
    x_doc.Functions.Func3Ds.Fn3D(9).Points.Points_(i).ATTRIBUTE.Z= RFc5(i);
    x_doc.Functions.Func3Ds.Fn3D(11).Points.Points_(i).ATTRIBUTE.Z= RFc6(i);
    x_doc.Functions.Func3Ds.Fn3D(13).Points.Points_(i).ATTRIBUTE.Z= RFc7(i);  
    x_doc.Functions.Func3Ds.Fn3D(15).Points.Points_(i).ATTRIBUTE.Z= RFc8(i);
    x_doc.Functions.Func3Ds.Fn3D(17).Points.Points_(i).ATTRIBUTE.Z= RFc9(i);
    x_doc.Functions.Func3Ds.Fn3D(19).Points.Points_(i).ATTRIBUTE.Z= RFc10(i);
    x_doc.Functions.Func3Ds.Fn3D(21).Points.Points_(i).ATTRIBUTE.Z= RFc11(i);
    x_doc.Functions.Func3Ds.Fn3D(23).Points.Points_(i).ATTRIBUTE.Z= RFc12(i);  
    x_doc.Functions.Func3Ds.Fn3D(25).Points.Points_(i).ATTRIBUTE.Z= RFc13(i);
    x_doc.Functions.Func3Ds.Fn3D(27).Points.Points_(i).ATTRIBUTE.Z= RFc14(i);
    x_doc.Functions.Func3Ds.Fn3D(29).Points.Points_(i).ATTRIBUTE.Z= RFc15(i);
    x_doc.Functions.Func3Ds.Fn3D(31).Points.Points_(i).ATTRIBUTE.Z= RFc16(i);
    x_doc.Functions.Func3Ds.Fn3D(33).Points.Points_(i).ATTRIBUTE.Z= RFc17(i);  
    x_doc.Functions.Func3Ds.Fn3D(35).Points.Points_(i).ATTRIBUTE.Z= RFc18(i);
    x_doc.Functions.Func3Ds.Fn3D(37).Points.Points_(i).ATTRIBUTE.Z= RFc19(i);
    x_doc.Functions.Func3Ds.Fn3D(39).Points.Points_(i).ATTRIBUTE.Z= RFc20(i);
end


% for j=1:1:1033
%     x_doc.Functions.Func3Ds.Fn3D(2).Points.Points_(j).ATTRIBUTE.Z=   RFphi(i);
%     x_doc.Functions.Func3Ds.Fn3D(4).Points.Points_(j).ATTRIBUTE.Z=   RFphi2(i);
%     x_doc.Functions.Func3Ds.Fn3D(6).Points.Points_(j).ATTRIBUTE.Z=   RFphi3(i);
%     x_doc.Functions.Func3Ds.Fn3D(8).Points.Points_(j).ATTRIBUTE.Z=   RFphi4(i);
%     x_doc.Functions.Func3Ds.Fn3D(10).Points.Points_(j).ATTRIBUTE.Z=  RFphi5(i);
%     x_doc.Functions.Func3Ds.Fn3D(12).Points.Points_(j).ATTRIBUTE.Z=  RFphi6(i);
%     x_doc.Functions.Func3Ds.Fn3D(14).Points.Points_(j).ATTRIBUTE.Z=  RFphi7(i);
%     x_doc.Functions.Func3Ds.Fn3D(16).Points.Points_(j).ATTRIBUTE.Z=  RFphi8(i);
%     x_doc.Functions.Func3Ds.Fn3D(18).Points.Points_(j).ATTRIBUTE.Z=  RFphi9(i);
%     x_doc.Functions.Func3Ds.Fn3D(20).Points.Points_(j).ATTRIBUTE.Z=  RFphi10(i);
%     x_doc.Functions.Func3Ds.Fn3D(22).Points.Points_(j).ATTRIBUTE.Z=  RFphi11(i);
%     x_doc.Functions.Func3Ds.Fn3D(24).Points.Points_(j).ATTRIBUTE.Z=  RFphi12(i);
%     x_doc.Functions.Func3Ds.Fn3D(26).Points.Points_(j).ATTRIBUTE.Z=  RFphi13(i);
%     x_doc.Functions.Func3Ds.Fn3D(28).Points.Points_(j).ATTRIBUTE.Z=  RFphi14(i);
%     x_doc.Functions.Func3Ds.Fn3D(30).Points.Points_(j).ATTRIBUTE.Z=  RFphi15(i);
%     x_doc.Functions.Func3Ds.Fn3D(32).Points.Points_(j).ATTRIBUTE.Z=  RFphi16(i);
%     x_doc.Functions.Func3Ds.Fn3D(34).Points.Points_(j).ATTRIBUTE.Z=  RFphi17(i);
%     x_doc.Functions.Func3Ds.Fn3D(36).Points.Points_(j).ATTRIBUTE.Z=  RFphi18(i);
%     x_doc.Functions.Func3Ds.Fn3D(38).Points.Points_(j).ATTRIBUTE.Z=  RFphi19(i);
%     x_doc.Functions.Func3Ds.Fn3D(40).Points.Points_(j).ATTRIBUTE.Z=  RFphi20(i);
% end
%%  
Pref.CellItem = false;
Pref.StructItem = false;
path = ['D:\Cheng\FS\', num2str(ind),'\Slope.xml'];
xml_write(path,x_doc,'GSIData',Pref);

%%   GeoStudio
xml_fn = ['D:\Cheng\FS\',num2str(ind), '\Slope.xml'];
cmdstr = ['SolveServer.exe ','"', xml_fn, '"'];
system(cmdstr)

%%   FS
pathout0 = ['D:\Cheng\FS\', num2str(ind)];
pathout1 = [pathout0, '\????????????????????????\001\slip_surface.csv'];
pathout2 = [pathout0, '\???????????????????????? (2)\001\slip_surface.csv'];
pathout3 = [pathout0, '\???????????????????????? (3)\001\slip_surface.csv'];
pathout4 = [pathout0, '\???????????????????????? (4)\001\slip_surface.csv'];
pathout5 = [pathout0, '\???????????????????????? (5)\001\slip_surface.csv'];
pathout6 = [pathout0, '\???????????????????????? (6)\001\slip_surface.csv'];
pathout7 = [pathout0, '\???????????????????????? (7)\001\slip_surface.csv'];
pathout8 = [pathout0, '\???????????????????????? (8)\001\slip_surface.csv'];
pathout9 = [pathout0, '\???????????????????????? (9)\001\slip_surface.csv'];
pathout10 = [pathout0, '\???????????????????????? (10)\001\slip_surface.csv'];
pathout11 = [pathout0, '\???????????????????????? (11)\001\slip_surface.csv'];
pathout12 = [pathout0, '\???????????????????????? (12)\001\slip_surface.csv'];
pathout13 = [pathout0, '\???????????????????????? (13)\001\slip_surface.csv'];
pathout14 = [pathout0, '\???????????????????????? (14)\001\slip_surface.csv'];
pathout15 = [pathout0, '\???????????????????????? (15)\001\slip_surface.csv'];
pathout16 = [pathout0, '\???????????????????????? (16)\001\slip_surface.csv'];
pathout17 = [pathout0, '\???????????????????????? (17)\001\slip_surface.csv'];
pathout18 = [pathout0, '\???????????????????????? (18)\001\slip_surface.csv'];
pathout19 = [pathout0, '\???????????????????????? (19)\001\slip_surface.csv'];
pathout20 = [pathout0, '\???????????????????????? (20)\001\slip_surface.csv'];
Output_Fs = csvread(pathout1,1,0); % read table skipping the first line
Fs = min(Output_Fs(:,2)) ;   %get the value from the table (it is saved in the 3rd column)
Output_Fs2 = csvread(pathout2,1,0); % read table skipping the first line
Fs2 = min(Output_Fs2(:,2)) ;
Output_Fs3 = csvread(pathout3,1,0); % read table skipping the first line
Fs3 = min(Output_Fs3(:,2)) ;
Output_Fs4 = csvread(pathout4,1,0); % read table skipping the first line
Fs4 = min(Output_Fs4(:,2)) ;
Output_Fs5 = csvread(pathout5,1,0); % read table skipping the first line
Fs5 = min(Output_Fs5(:,2)) ;
Output_Fs6 = csvread(pathout6,1,0); % read table skipping the first line
Fs6 = min(Output_Fs6(:,2)) ;   %get the value from the table (it is saved in the 3rd column)
Output_Fs7 = csvread(pathout7,1,0); % read table skipping the first line
Fs7 = min(Output_Fs7(:,2)) ;
Output_Fs8 = csvread(pathout8,1,0); % read table skipping the first line
Fs8 = min(Output_Fs8(:,2)) ;
Output_Fs9 = csvread(pathout9,1,0); % read table skipping the first line
Fs9 = min(Output_Fs9(:,2)) ;
Output_Fs10 = csvread(pathout10,1,0); % read table skipping the first line
Fs10 = min(Output_Fs10(:,2)) ;
Output_Fs11 = csvread(pathout11,1,0); % read table skipping the first line
Fs11 = min(Output_Fs11(:,2)) ;   %get the value from the table (it is saved in the 3rd column)
Output_Fs12 = csvread(pathout12,1,0); % read table skipping the first line
Fs12 = min(Output_Fs12(:,2)) ;
Output_Fs13 = csvread(pathout13,1,0); % read table skipping the first line
Fs13 = min(Output_Fs13(:,2)) ;
Output_Fs14 = csvread(pathout14,1,0); % read table skipping the first line
Fs14 = min(Output_Fs14(:,2)) ;
Output_Fs15 = csvread(pathout15,1,0); % read table skipping the first line
Fs15 = min(Output_Fs15(:,2)) ;
Output_Fs16 = csvread(pathout16,1,0); % read table skipping the first line
Fs16 = min(Output_Fs16(:,2)) ;   %get the value from the table (it is saved in the 3rd column)
Output_Fs17 = csvread(pathout17,1,0); % read table skipping the first line
Fs17 = min(Output_Fs17(:,2)) ;
Output_Fs18 = csvread(pathout18,1,0); % read table skipping the first line
Fs18 = min(Output_Fs18(:,2)) ;
Output_Fs19 = csvread(pathout19,1,0); % read table skipping the first line
Fs19 = min(Output_Fs19(:,2)) ;
Output_Fs20 = csvread(pathout20,1,0); % read table skipping the first line
Fs20 = min(Output_Fs20(:,2)) ;
% rmdir(pathout0);
end