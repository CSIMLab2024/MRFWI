clear;
clc; 
close all
load model.mat
%% model paramaters
dt1=1.0e-11;
T=3e-8;
f=400e6;
srcway=1; %1 ricker  |0 blackharrispulse
dxway=1; %1 noStability  |0 Stability
dx1=0.01;
[dx,dz,dt,srcpulse,nt_wf,dt_wf,z,x,nz,nx,t,vpmin,vpmax,vpmine, vpmaxe,nt] = Para_PRE(ep,sig,mu,T,dt1,srcway,dxway,f,dx1);
%% model antenna
% offset=30;
% ds=20;
% dg=20;
offset=20;
ds=5;
dg=5;
air_d=20;

gprmode=1; %1 common offset  |0 multiple offset
lm=1;
bj1=0;
bj2=0;
[recloc,srcloc,ng,nsrc] = bj_Aat_PRE(nx,dx,ds,dg,offset,air_d,gprmode,lm,bj1,bj2);
%% Forward
FWImode=11; %10 MR-FWI  |00 CFWI | 11 parallel MR-FWI  |01 parallel CFWI
nnoise=0; % noise  
parx=59;
%multiple offset   Forward_PRE_mu 
%common offset  Forward_PRE_co 
% [gatherobs,wavefield_gradient] = bj_Forward_PRE_mu(recloc,srcloc,ng,nsrc,srcpulse,dx,ep,sig,mu,T,dt,dz,x,z,nnoise,parx,nt_wf,dt_wf,t,FWImode);

[gatherobs] = bj_Forward_PRE_co1(recloc,srcloc,ng,nsrc,srcpulse,dx,ep,sig,mu,T,dt,dz,x,z,nnoise,parx,nt_wf,dt_wf,t,FWImode,bj1,bj2);
cd(:,:)=gatherobs(:,1,:);


%% Inversion
 ff=1e6*[100 200 250 300 400 500];
iteration=250;
FWImode=11; %10 MR-FWI  |00 CFWI | 11 parallel MR-FWI  |01 parallel CFWI
k =1;
kk=0;
n=1;
err=0.05;
N=find(ep(1,:)>1);
f1=0.5;
szz=5;
sxx=10;
pur=0.05;
%
lam=1.0e-6;
k =1;kk=0;n=1;err=0.05;%k
N=find(ep(1,:)>1);
iteration=250;
 pur=0.05;
  outstep = [1, 1, 1]; % 1st for time sampling, 2nd for spatial sampling along x, 3rd for spatial sampling along z;
 plotopt = [1 1 0.001]; 
 plotopt_back = [plotopt(1:2), plotopt(3)^2 * 10];
 dt_wf=dt*2; nt_wf=floor((nt-1)*dt/dt_wf)+1;
%% ++++++++Start iteration 
while (k<=iteration)    
srcway=1; %1 ricker  |0 blackharrispulse
[gatherobs1,srcpulse1] = Mulitpe_frequency(ff,nsrc,dt,nt,ng,gatherobs,srcpulse,srcway,n);
%%
parfor isrc = 1:nsrc
    [ep_sse,mu_sse, sig_sse,xse,zse,zz] = Mr1(bj1,bj2,srcloc,recloc,dx,ep_ss,mu_ss,sig_ss,isrc,lm);
     [forwavefield,gatherpre]=TMGenerateData(srcloc(zz,:),recloc(zz,:),  srcpulse1,T,isrc,ep_sse,mu_sse,sig_sse,dt,dx,dz,xse,zse,outstep,plotopt,dt_wf,nt_wf);
      [joint_source,delta_T] = residual_wt(gatherpre,gatherobs1(:,:,isrc));%%%%求差residual_wt
      res(isrc)=sum(delta_T(:).*delta_T(:));
      joint_source=joint_source';
        bpwavefield=TMRunBackwardnew1_hr(recloc(zz,:),srcloc(zz,:),joint_source,T,isrc,ep_sse,mu_sse,sig_sse,dt,dx,dz,xse,zse,outstep,plotopt,t,dt_wf);
       [a_dk11(:,:,isrc),a_dk1e1(:,:,isrc),~,a_illum1(:,:,isrc)]=image_conditionnew1(forwavefield,bpwavefield,sig_sse,ep_sse,dx,dt);  
end
   
%% calcuate the gradient 
[a_dk1e,a_dk1] = Recon1(nx,nz,nsrc,srcloc,recloc,bj1,bj2,dx,lm,a_dk11,a_dk1e1,a_illum1);

% Use the penality method to build the objective function
    res=sum(res); 
    residual(k)=res;
    res0=residual(k); 
    display(['residual = ',num2str( residual(k) ),' k=',num2str(k)]);    

    a_dk1 = smooth2a(a_dk1,3,4);
   a_dk1e = smooth2a(a_dk1e,3,4);

    clear a_wavefield bp_wavefield

 a_dk1=norm_shot(a_dk1);
   a_dk1e=norm_shot(a_dk1e);

sigg=sig_ss;
     epp=ep_ss;
     sig_ss=norm_shot(sig_ss);

    ep_ss=norm_shot(ep_ss);
     a_s_mean=(sum(sig_ss(:).*sig_ss(:)))^0.5;
     a_g_mean=(sum(a_dk1(:).*a_dk1(:)))^0.5;
     a_alpha=a_s_mean/a_g_mean*pur;
    
    a_s_meane=(sum(ep_ss(:).*ep_ss(:)))^0.5;
    a_g_meane=(sum(a_dk1e(:).*a_dk1e(:)))^0.5;
    a_alphae=a_s_meane/a_g_meane*pur;

 % Back tracking to find the numerical step length
    if k==1
    f1=0.5;
    end
    %% total_variation
%     am=0;
%  %a_dk1=-a_dk1e;
% % am=sum(a_dk1(:))./sum(a_dk1e(:));
% L =   5.0;
% g =   5.0;
% it =    20;
% S =     1;
%     a_dk1v = total_variation(a_dk1,g,L,it,S);
%     %[a_dk1v]=Split_TV(a_dk1);
%     a_dk1ev = total_variation(a_dk1e,g,L,it,S);
%     lam1=norm(a_dk1)/norm(a_dk1v);
%     a_dk1v=a_dk1v*lam1;
%      lam1=norm(a_dk1e)/norm(a_dk1ev);
%      a_dk1ev=a_dk1ev*lam1;
     
     a_dk1(:,1:N(1))=0; a_dk1e(:,1:N(1))=0;
     a_dk1v(:,1:N(1))=0; a_dk1ev(:,1:N(1))=0;
    sig11=sig_ss+a_alpha*f1*a_dk1;
    
   
     ep11=ep_ss+a_alphae*f1*a_dk1e;
   

      sig11=invnorm_shot(sig11,sigg);
      ep11=invnorm_shot(ep11,epp);
sig11(sig11<vpmin)=vpmin;sig11(sig11>vpmax)=vpmax;
  ep11(ep11<vpmine)=vpmine;ep11(ep11>vpmaxe)=vpmaxe;


parfor isrc = 1:nsrc
      [ep11e,mu_sse,sig11e,xse,zse,zz] = Mr1(bj1,bj2,srcloc,recloc,dx,ep11,mu,sig11,isrc,lm);
    [~,gatherpre]=TMGenerateData(srcloc(zz,:),recloc(zz,:), srcpulse1,T,isrc,ep11e,mu_sse,sig11e,dt,dx,dz,xse,zse,outstep,plotopt,dt_wf,nt_wf);
[joint_source,delta_T] = residual_wt(gatherpre,gatherobs1(:,:,isrc));
      res1(isrc)=sum(delta_T(:).*delta_T(:));
end

res1=sum(res1);
  display(['f1= ',num2str(f1),' res1= ',num2str(res1)]);
    if res1>res0
        while res1>res0 && f1>0.0001
            f2=f1; res2=res1;
            f1=f1*0.5;
    sig11=sig_ss+a_alpha*f1*a_dk1;
    
   
    ep11=ep_ss+a_alphae*f1*a_dk1e;
      sig11=invnorm_shot(sig11,sigg);
    ep11=invnorm_shot(ep11,epp);
  sig11(sig11<vpmin)=vpmin;sig11(sig11>vpmax)=vpmax;
 ep11(ep11<vpmine)=vpmine;ep11(ep11>vpmaxe)=vpmaxe;


parfor isrc = 1:nsrc
         [ep11e,mu_sse,sig11e,xse,zse,zz] = Mr1(bj1,bj2,srcloc,recloc,dx,ep11,mu,sig11,isrc,lm);
[~,gatherpre]=TMGenerateData(srcloc(zz,:),recloc(zz,:),srcpulse1,T,isrc,ep11e,mu_sse,sig11e,dt,dx,dz,xse,zse,outstep,plotopt,dt_wf,nt_wf);
    [joint_source,delta_T] = residual_wt(gatherpre,gatherobs1(:,:,isrc));
      res1(isrc)=sum(delta_T(:).*delta_T(:));
end
 
res1=sum(res1);

    display(['f1= ',num2str(f1),' res1= ',num2str(res1)]);
        end
    else
        f2=f1*2;

    sig11=sig_ss+a_alpha*f2*a_dk1;
    
   
   ep11=ep_ss+a_alphae*f2*a_dk1e;
     sig11=invnorm_shot(sig11,sigg);
    ep11=invnorm_shot(ep11,epp);
sig11(sig11<vpmin)=vpmin;sig11(sig11>vpmax)=vpmax;
 ep11(ep11<vpmine)=vpmine;ep11(ep11>vpmaxe)=vpmaxe;

parfor isrc = 1:nsrc
        [ep11e,mu_sse,sig11e,xse,zse,zz] = Mr1(bj1,bj2,srcloc,recloc,dx,ep11,mu,sig11,isrc,lm);
 [~,gatherpre]=TMGenerateData(srcloc(zz,:),recloc(zz,:),srcpulse1,T,isrc,ep11e,mu_sse,sig11e,dt,dx,dz,xse,zse,outstep,plotopt,dt_wf,nt_wf);
[~,delta_T] = residual_wt(gatherpre,gatherobs1(:,:,isrc));
      res2(isrc)=sum(delta_T(:).*delta_T(:));
end
 
res2=sum(res2);
        display(['f2= ',num2str(f2),' res2= ',num2str(res2)]);

end
    gama=(f1^2*(res0-res2)+f2^2*(res1-res0))/(2*res0*(f1-f2)+2*res1*f2-2*res2*f1);
    display(['gama= ',num2str(gama),' numerical step_length= ',num2str(gama*a_alphae)]);
    
sig11=sig_ss+a_alpha*gama*a_dk1;
       ep11=ep_ss+a_alphae*gama*a_dk1e;
      sig11=invnorm_shot(sig11,sigg);
     ep11=invnorm_shot(ep11,epp);
sig11(sig11<vpmin)=vpmin;sig11(sig11>vpmax)=vpmax;
 ep11(ep11<vpmine)=vpmine;ep11(ep11>vpmaxe)=vpmaxe;

parfor isrc = 1:nsrc
       [ep11e,mu_sse,sig11e,xse,zse,zz] = Mr1(bj1,bj2,srcloc,recloc,dx,ep11,mu,sig11,isrc,lm);

      [~,gatherpre]=TMGenerateData(srcloc(zz,:),recloc(zz,:), srcpulse1,T,isrc,ep11e,mu_sse,sig11e,dt,dx,dz,xse,zse,outstep,plotopt,dt_wf,nt_wf);
[~,delta_T] = residual_wt(gatherpre,gatherobs1(:,:,isrc));
      res3(isrc)=sum(delta_T(:).*delta_T(:));
end
 
res3=sum(res3);
 

    display(['res3= ',num2str(res3)]);
    if (res3>res1 || res3>res2)
        if res1>res2
            res0=res2;

            gama=f2; 
        else
            res0=res1;
            gama=f1; 
        end
        sig11=sig_ss+a_alpha*gama*a_dk1;     
       ep11=ep_ss+a_alphae*gama*a_dk1e;
      sig11=invnorm_shot(sig11,sigg);
     ep11=invnorm_shot(ep11,epp);
sig11(sig11<vpmin)=vpmin;sig11(sig11>vpmax)=vpmax;
 ep11(ep11<vpmine)=vpmine;ep11(ep11>vpmaxe)=vpmaxe;

    else
        res0=res3;
    end
    sig_ss=sig11;
    siggg(:,:,k)= sig_ss;
    ep_ss=ep11;
       epppp(:,:,k)= ep_ss;
  
   if k>1 && ((residual(k) - res0)/residual(k)) < err 
   n = n + 1
   f1=0.5;
   end

   k = k +1
   clear res a_dk11 a_dk1e1 a_illum1
 end






