%n = 124;
b = zeros(180,20);
k = 0;
tic
%%
% x=randn(10000,80);  %mc
%y=randn(1000,9);
%% lhs
p=900;
q=80;
x = lhsdesign(p,q); 
for i=1:q
    x(1:p,i)=sqrt(2)*1*erfinv(2*x(:,i)-1);
end

%%
parfor i=0:179
    fprintf('i = %d\n',i);
    ii = i*20;
    z=x(ii+1,:).';
    z2=x(ii+2,:).';
    z3=x(ii+3,:).';
    z4=x(ii+4,:).';
    z5=x(ii+5,:).';
    z6=x(ii+6,:).';
    z7=x(ii+7,:).';
    z8=x(ii+8,:).';
    z9=x(ii+9,:).';
    z10=x(ii+10,:).';
    z11=x(ii+11,:).';
    z12=x(ii+12,:).';
    z13=x(ii+13,:).';
    z14=x(ii+14,:).';
    z15=x(ii+15,:).';
    z16=x(ii+16,:).';
    z17=x(ii+17,:).';
    z18=x(ii+18,:).';
    z19=x(ii+19,:).';
    z20=x(ii+20,:).';
%     zz=y(ii+1,:).';
%     zz2=y(ii+2,:).';
%     zz3=y(ii+3,:).';
%     zz4=y(ii+4,:).';
%     zz5=y(ii+5,:).';
%     zz6=y(ii+6,:).';
%     zz7=y(ii+7,:).';
%     zz8=y(ii+8,:).';
%     zz9=y(ii+9,:).';
%     zz10=y(ii+10,:).';
%     zz11=y(ii+11,:).';
%     zz12=y(ii+12,:).';
%     zz13=y(ii+13,:).';
%     zz14=y(ii+14,:).';
%     zz15=y(ii+15,:).';
%     zz16=y(ii+16,:).';
%     zz17=y(ii+17,:).';
%     zz18=y(ii+18,:).';
%     zz19=y(ii+19,:).';
%     zz20=y(ii+20,:).';
%     [Fs,Fs2,Fs3,Fs4,Fs5,Fs6,Fs7,Fs8,Fs9,Fs10,Fs11,Fs12,Fs13,Fs14,Fs15,Fs16,Fs17,Fs18,Fs19,Fs20]...
%      =Main_RF(z,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z20,zz,zz2,zz3,zz4,zz5,zz6,zz7,zz8,zz9,zz10,...
%      zz11,zz12,zz13,zz14,zz15,zz16,zz17,zz18,zz19,zz20,i);
 [Fs,Fs2,Fs3,Fs4,Fs5,Fs6,Fs7,Fs8,Fs9,Fs10,Fs11,Fs12,Fs13,Fs14,Fs15,Fs16,Fs17,Fs18,Fs19,Fs20]...
     =Main_RF(z,z2,z3,z4,z5,z6,z7,z8,z9,z10,z11,z12,z13,z14,z15,z16,z17,z18,z19,z20,i);
    b(i+1,:) = [Fs,Fs2,Fs3,Fs4,Fs5,Fs6,Fs7,Fs8,Fs9,Fs10,Fs11,Fs12,Fs13,Fs14,Fs15,Fs16,Fs17,Fs18,Fs19,Fs20];
end
toc