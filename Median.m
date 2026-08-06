clear,clc


disp('Step 2: IFMs')

mu(:,:,2)=[0.75 0.85 0.50
    0.65 0.50  0.65
    0.50 0.85  0.75
    0.75 0.50  0.65];

nu(:,:,2)=[0.15 0.10 0.40
     0.25  0.40 0.25
     0.40  0.10 0.15
     0.15  0.40 0.25];
          
mu(:,:,1)=[0.35 0.50 0.25
    0.25 0.05  0.50
    0.75 0.65  0.35
    0.50 0.50  0.65];

nu(:,:,1)=[0.55 0.40 0.65
    0.65 0.95  0.40
    0.15 0.25  0.55
    0.40 0.40  0.25];

mu(:,:,4)=[0.85 0.75 0.75
    0.95 0.65  0.75 
    0.75 0.85  0.65
    0.65 0.65  0.75];

nu(:,:,4)=[0.10 0.15 0.15
    0.05 0.25  0.15
    0.15 0.10  0.25
    0.25 0.25  0.15]; 

          
mu(:,:,3)=[0.25  0.35  0.50
    0.50  0.65  0.35
    0.15  0.50  0.75
    0.75 0.50  0.65];
  
nu(:,:,3)=[0.65  0.55  0.40
    0.40  0.25  0.55
    0.80  0.40  0.15
    0.15  0.40  0.25];

[m,n]=size(mu(:,:,1));
t=4;

for k=1:t
    temp=[reshape(mu(:,:,k)',1,[]); reshape(nu(:,:,k)',1,[])];
   eval(['Y',num2str(k),'=  sprintf(''(%4.2f, %4.2f) (%4.2f, %4.2f) (%4.2f, %4.2f) \n'', temp)'])
end
clear Y1 Y2 Y3  Y4

mm=median(mu,3);

mn=median(nu,3);

disp('步骤3，F的mode matrix：')

% temp=[reshape(mode_m',1,[]); reshape(mode_n',1,[])];

mod=[];
mod=[mod,mm(:,1),mn(:,1),mm(:,2),mn(:,2), mm(:,3),mn(:,3)];

Median_m=sprintf('(%4.2f, %4.2f) (%4.2f, %4.2f) (%4.2f, %4.2f)  \n', mod')

%%% 第3步，求Fk的逆序数 

for k=1:t
     Zm(:,:,k)=reshape(mu(:,:,k)',1,12);   %将x转成1*16（行接行）矩阵
     Tm{k}=Zm(:,:,k);
 %  temp=[reshape(Zm(:,:,k)',1,[])];
%   eval(['Zm',num2str(k),'=  sprintf(''%4.2f  %4.2f  %4.2f \n'', temp)'])
%    temp=[reshape(Z(:,:,k)',1,[])];
%   eval(['Z',num2str(k),'=  sprintf(''%4.2f,  %4.2f %4.2f,  %4.2f \n'', temp)'])
end
% disp('步骤4，mu(k)排成一行为：')
Zm;
%% 步骤5：Calculate the number of inversions of Y(k).

% [m,n]=size(R(:,:,k));

disp('muk的每个元素所产生的逆序数为：')  
for k=1:t
%t{k}=R(:,:,k);      
sm{k}=0;
     for j=2:length(Tm{k})     
     um{k}=sum((Tm{k}(1:j-1)>Tm{k}(j)));
     sm{k}=[sm{k},um{k}];   %?
     end
 
%  temp=[reshape(sm{k}',1,[])];
%  eval(['sm',num2str(k),'=  sprintf(''%4.0f  %4.0f  %4.0f  \n'', temp)'])
 %s{k}    %各元素的逆序数向量；
inv_num_m(k)=sum(sm{k});   %向量的逆序数
end
inv_num_m   % ?

for k=1:t
    Zn(:,:,k)=reshape(nu(:,:,k)',1,12);   %将x转成1*16（行接行）矩阵
    Tn{k}=Zn(:,:,k);
%    temp=[reshape(Tn{k}',1,[])];
%    eval(['Tn',num2str(k),'=  sprintf(''%4.2f  %4.2f  %4.2f \n'', temp)'])
%    temp=[reshape(Z(:,:,k)',1,[])];
%   eval(['Z',num2str(k),'=  sprintf(''%4.2f,  %4.2f %4.2f,  %4.2f \n'', temp)'])
end
% disp('步骤4，nu(k)排成一行为：')
Zn;


%% 步骤5：Calculate the number of inversions of Y(k).

% [m,n]=size(R(:,:,k));
disp('nuk的每个元素所产生的逆序数为：')  
for k=1:t
%t{k}=R(:,:,k);      
sn{k}=0;
     for j=2:length(Tn{k})     
     un{k}=sum((Tn{k}(1:j-1)>Tn{k}(j)));
     sn{k}=[sn{k},un{k}];   %?
     end
 
%  temp=[reshape(sn{k}',1,[])];
%  eval(['sn',num2str(k),'=  sprintf(''%4.0f  %4.0f  %4.0f  \n'', temp)'])
 %s{k}    %各元素的逆序数向量；
 
inv_num_n(k)=sum(sn{k});   %向量的逆序数
end

inv_num_n   % ?

disp('总的的逆序数为：')  

inv_num=inv_num_m+inv_num_n


% 步骤6：Calculate the number of inversions of F(median).

%先转成1*16（行接行）矩阵

% disp('步骤5，Y(*)排成一行：')

Y_m=reshape(mm',1,12);

smo=0;
for j=2:length(Y_m)     
     um=sum((Y_m(1:j-1)>Y_m(j)));
     smo=[smo,um];   %?
end
     disp('Ymode的mu每个元素所产生的逆序数为：') 
% smo   %各元素的逆序数向量；
% Reshape smo into a 4x3 matrix
 smo_matrix = reshape(smo', 4, 3);
 disp('The smo matrix is:');
smo_matrix
inv_num_mo=sum(smo)  %平均矩阵的逆序数

Y_n=reshape(mn',1,12);

sno=0;
for j=2:length(Y_n)     
     un=sum((Y_n(1:j-1)>Y_n(j)));
     sno=[sno,un];   %?
end
     disp('Ymode的nu每个元素所产生的逆序数为：') 
 % sno   %各元素的逆序数向量；
   sno_matrix = reshape(sno', 4, 3);
 disp('The sno matrix is:');
sno_matrix
inv_num_no=sum(sno)  %平均矩阵的逆序数

% temp=[reshape(smo',1,[]); reshape(sno',1,[])];
%    eval(['Y_modeInv',num2str(k),'=  sprintf(''(%1.0f, %1.0f) (%1.0f, %1.0f) (%1.0f, %1.0f) \n'', temp)'])


inv_num_mode=inv_num_mo+inv_num_no


disp('步骤3，Yk分熵值为：')
 
 for k=1:t 
    Pm(:,:,k)=mu(:,:,k)./sum(sum(mu(:,:,k)));
    Em(k)=-sum(sum((Pm(:,:,k).*log(Pm(:,:,k)))))./log(m*n);   
    
    Pn(:,:,k)=nu(:,:,k)./sum(sum(nu(:,:,k)));
    En(k)=-sum(sum((Pn(:,:,k).*log(Pn(:,:,k)))))./log(m*n); 
 end
 
 Em
 
 En
 
 
disp('步骤3，Yk总熵值为：')
 
 
 entropyk=Em+En
 

 Pme=mm./sum(sum(mm));
 Eme=-sum(sum((Pme.*log(Pme))))./log(m*n)   
    
 Pne=mn./sum(sum(mn));
 Ene=-sum(sum((Pne.*log(Pne))))./log(m*n) 
  
  disp('E*的总熵值为：')
  
  entropy_median=Eme+Ene
  
   disp('基于逆序的CCs为：')
   
    IC=inv_num_mode./(inv_num_mode+abs(inv_num-inv_num_mode))
  
  
 disp('基于熵的CCs为：')
 
 EC=entropy_median./(entropy_median+abs(entropyk-entropy_median))


%% 步骤7： Calculate the relative closeness of inversions.

disp('Step 6,integrated CC：')

   
    CC=(IC+EC)/2
    
%% 步骤8： Calculate the inversion-based weights of DMs

disp('Step 7, weights of DMs 为：')

 lambda=CC/sum(CC)
          

%%%% 第9步，求加权决策矩阵

disp('Step 9，加权到决策者')

for k=1:4 
    % weight(:,:,k)=ones(m,3)*lambda(k)
     tau(:,:,k)=1-(1-mu(:,:,k)).^lambda(k);
   % tau(:,:,k)=mu(:,:,k).^weight(:,:,k);
   % upsilon(:,:,k)=sqrt(1-(1-nu(:,:,k).^2).^weight(:,:,k));  
     upsilon(:,:,k)=nu(:,:,k).^lambda(k);
   tmp=[reshape(tau(:,:,k)',1,[]);  reshape(upsilon(:,:,k)',1,[])];
   eval(['F',num2str(k),'=  sprintf(''(%4.2f, %4.2f)  (%4.2f, %4.2f)  (%4.2f, %4.2f)\n'', tmp)']);
end



%%%% 转换成方案矩阵

disp('第10步，群决策矩阵')

xii=permute(tau,[3,2,1]);
oi=permute(upsilon,[3,2,1]);
for k=1:4    
   temp=[reshape(xii(:,:,k)',1,[]);  reshape(oi(:,:,k)',1,[])];
   eval(['H',num2str(k),'=  sprintf(''(%4.2f, %4.2f)  (%4.2f, %4.2f)  (%4.2f, %4.2f)  \n'', temp)']);
end

disp('第11步，加权到属性')

weights=[0.4,0.4,0.2];

for k=1:4      
    % 正确：对每个属性j使用对应的权重weights(j)
    for j=1:3
        xi(:,j,k)=1-(1-xii(:,j,k)).^weights(j);
        o(:,j,k)=oi(:,j,k).^weights(j);
    end
    temp=[reshape(xi(:,:,k)',1,[]);  reshape(o(:,:,k)',1,[])];
   eval(['G',num2str(k),'=  sprintf(''(%4.2f, %4.2f) & (%4.2f, %4.2f) & (%4.2f, %4.2f)  \n'', temp)']);
end
  

%%%% 第3步，求理想解

xi_nega=min(xi,[],3);%%%负理想解
xi_posi=max(xi,[],3);%%%正理想解


o_nega=max(o,[],3);%%%负理想解
o_posi=min(o,[],3);%%%正理想解


xi_nega;

xi_posi;

o_nega;

o_posi;

posi = [];
nega = [];

posi=[posi, xi_posi(:,1), o_posi(:,1),xi_posi(:,2), o_posi(:,2),xi_posi(:,3), o_posi(:,3) ];

nega=[nega, xi_nega(:,1), o_nega(:,1),xi_nega(:,2), o_nega(:,2),xi_nega(:,3), o_nega(:,3)];

disp('第12步，理想解')


%%% 输出正负理想解 

disp('正理想解为：')

G_posi=sprintf('(%4.2f, %4.2f) & (%4.2f, %4.2f) & (%4.2f, %4.2f)\n', posi')

disp('负理想解为：')

G_nega=sprintf('(%4.2f, %4.2f) & (%4.2f, %4.2f) & (%4.2f, %4.2f)\n', nega')

%%% 参照矩阵的模平方

G_quare_posi=sum(sum(xi_posi.^2+o_posi.^2+(1-xi_posi-o_posi).^2 )); %%G+模平方

G_quare_nega=sum(sum(xi_nega.^2+o_nega.^2+(1-xi_nega-o_nega).^2 ));  %%G-模平方
   
for k=1:t
   
    G_G_posi(k)=sum(sum(xi(:,:,k).*xi_posi+o(:,:,k).*o_posi+(1-xi(:,:,k)-o(:,:,k)).*(1-xi_posi-o_posi)));  %%G与G_{+}内积
    G_G_nega(k)=sum(sum(xi(:,:,k).*xi_nega+o(:,:,k).*o_nega+(1-xi(:,:,k)-o(:,:,k)).*(1-xi_nega-o_nega)));  %%G与G_{-}内积 
    G_square(k)=sum(sum(xi(:,:,k).^2+o(:,:,k).^2+(1-xi(:,:,k)-o(:,:,k)).^2 ));%%G模平方
    
    a(k)=abs(1-G_G_posi(k)./G_square(k)-G_G_posi(k)./G_quare_posi);
    b(k)=a(k)+abs(1-G_G_posi(k)./G_quare_posi);
    
    
    an(k)=abs(1-G_G_nega(k)./G_square(k)-G_G_nega(k)./G_quare_nega);
    %an(k)=abs(1-G_G_nega(k)./G_G_nega(k)-G_G_nega(k)./G_quare_nega);
    bn(k)=an(k)+abs(1-G_G_nega(k)./G_quare_nega);
        

%    NP_G_G_posi(k)=abs(1-G_G_posi(k)./G_square(k)-G_G_posi(k)./G_quare_posi)./(abs(1-G_G_posi(k)./G_square(k)-G_G_posi(k)./G_quare_posi)+abs(1-G_G_posi(k)./G_quare_posi));   %%G在G_{+}上标准化投影
%    NP_G_G_nega(k)=abs(1-G_G_nega(k)./G_square(k)-G_G_nega(k)./G_quare_nega)./(abs(1-G_G_nega(k)./G_square(k)-G_G_nega(k)./G_quare_nega)+abs(1-G_G_nega(k)./G_quare_nega));   %%G在G_{-}上标准化投影

NP_G_G_posi(k)=a(k)./b(k);
   
%     ap(k)=abs(1-G_G_posi(k)./G_square(k)-G_G_posi(k)./G_quare_posi);
%     bp(k)=ap(k)+abs(1-G_G_posi(k)./G_quare_posi);
%    cp=ap./bp
   
%     a(k)=abs(1-G_G_nega(k)./G_G_nega(k)-G_G_nega(k)./G_quare_nega);
%    b(k)=a(k)+abs(1-G_G_nega(k)/G_quare_nega);
%     c=a./b
   
 %  d(k)=abs(1-G_G_nega(k)./G_G_nega(k)-G_G_nega(k)./G_quare_nega)+abs(1-G_G_nega(k)./G_quare_nega);
 
 NP_G_G_nega(k)=an(k)./bn(k);
   
   
   URC=NP_G_G_posi./(NP_G_G_posi+NP_G_G_nega);
 % vpa(URC,5)
   
end

 NP_G_G_nega;

% a, b, c, cp
% 
% ap, bp, cp

% NP_nega=abs(1-G_G_nega./G_G_nega-G_G_nega./G_quare_nega)./(abs(1-G_G_nega./G_square-G_G_nega./G_quare_nega)+abs(1-G_G_nega./G_quare_nega))


disp('Step 17, group utility measurement')



 NP_G_G_posi
 
 NP_G_G_nega
 
 %NP_G_G_ne

 URC
