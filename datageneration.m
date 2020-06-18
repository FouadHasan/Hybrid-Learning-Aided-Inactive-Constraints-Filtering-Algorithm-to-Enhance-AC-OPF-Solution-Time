clc;
clear all;
close all;
define_constants;
scenario=1;%number of scenario

%percentage(+/-) of load variation
pc=30;
pc=pc/100; 
%loading the system
% file='pglib_opf_case39_epri';
% file='pglib_opf_case118_ieee';
% file='pglib_opf_case300_ieee';
% file='pglib_opf_case500_tamu'; 
file='pglib_opf_case1354_pegase';

mpc=loadcase(file);

nb=length(mpc.bus(:,1));% number of buses
nl=length(mpc.branch(:,1));% number of lines
ng=length(mpc.gen(:,1));

mpc = ext2int(mpc);

% mpc.branch(:,6:8)=mpc.branch(:,6:8)*0.75;
%extracting load buses
j=1;
for i=1:nb
 if (mpc.bus(i,PD)~=0)
    lb(j)=mpc.bus(i,1);
    j=j+1;
 end
end

%extracting generator buses
 gb=find(mpc.bus(:,2)~=1);
% gb=mpc.gen(:,1);
% ng=length(gb);% number of generators

%storing the base values of load
p=mpc.bus(lb,PD);%real power
q=mpc.bus(lb,QD);%reactive power
pf=cos(atan(q./p));%power factor

%Lower limit of load
pl=p*0.70;
ql=q*0.70;
%upper limit of load
pu=p*1.1;
qu=q*1.1;
%-------------------------------------------------------------------------------------------------------------
%%generating scebnarios(random)
busdata=zeros(nb*(scenario),6);
gendata=zeros(ng*(scenario),7);
branchdata=zeros(nl*(scenario),4);
data=zeros(nb,4,(scenario)); 

mlinputbus=zeros(nb*(scenario),nb*2);
mlinputbranch=zeros(nl*(scenario),nb*2);
mlinputgen=zeros(ng*(scenario),nb*2);


for i=1:scenario
    pd=(pl+(pu-pl).*rand(length(lb),1));
%     qd=pd.*tan(acos(pf));
    qd=(ql+(qu-ql).*rand(length(lb),1));
    mpc.bus(lb,PD)=pd;
    mpc.bus(lb,QD)=qd;
%   disp(num2str([p q pd qd pf cos(atan(qd./pd))],'%.2f'))
    mpopt = mpoption('out.all', 0);
    results=runopf(mpc,mpopt);
    time(i)=results.et;
    busdata((i-1)*nb+1:i*nb,:)=results.bus(:,[1,3,4,8,16,17]);%[bus no,p,q,v,MU_vmax,MU_vmin]
    gendata((i-1)*ng+1:i*ng,:)=results.gen(:,[1,2,3,22:25]);%[gen no,pg,qg,MU_PMAX,MU_PMIN,MU_QMAX,MU_QMAX]
    branchdata((i-1)*nl+1:i*nl,:)=results.branch(:,[1,2,18,19]);%[from bus,to bus,MU_SF, MU ST]
    
    %when some buses have multiple generators; use this line
    %"pq_g=[mpc.gen(:,1) results.gen(:,[PG,QG])];" and follwing for loop.
    % Comment out the "data(gb,1:2,i)=results.gen(:,[PG,QG]);"
 %%   
%     pq_g=[mpc.gen(:,1) results.gen(:,[PG,QG])];% storing Pg Qg 
%     
%     for m=1:length(gb)
%     a=find(pq_g(:,1)==gb(m));
%     data(gb(m),1,i)=sum(pq_g(a,2));
%     data(gb(m),2,i)=sum(pq_g(a,3));
%     end
 %% 
    %when each bus has only one generator; use the following line  "data(gb,1:2,i)=results.gen(:,[PG,QG]);" and 
    %comment out the previous for loop and line "pq_g=[mpc.gen(:,1) results.gen(:,[PG,QG])];"   
 %%
    data(gb,1:2,i)=results.gen(:,[PG,QG]);
   %% 
     
    data(:,3:4,i)=results.bus(:,[PD,QD]);%storing Pd Qd
    
    Pgd=zeros(nb,nb);
    Qgd=zeros(nb,nb);
    for k=1:nb
    Pgd(k,:)=((data(:,1,i)-data(:,3,i))');
    Qgd(k,:)=((data(:,2,i)-data(:,4,i))');
    end
    mlinputbus((i-1)*nb+1:i*nb,1:nb)=Pgd;
    mlinputbus((i-1)*nb+1:i*nb,(nb+1):2*nb)=Qgd;
    
    Pgd=[];
    Qgd=[];
    Pgd=zeros(nl,nb);
    Qgd=zeros(nl,nb);
    for k=1:nl
    Pgd(k,:)=(data(:,1,i)-data(:,3,i))';
    Qgd(k,:)=(data(:,2,i)-data(:,4,i))';
    end
    mlinputbranch((i-1)*nl+1:i*nl,1:nb)=Pgd;
    mlinputbranch((i-1)*nl+1:i*nl,(nb+1):2*nb)=Qgd;
    
    Pgd=[];
    Qgd=[];
    Pgd=zeros(ng,nb);
    Qgd=zeros(ng,nb);
    for k=1:ng
    Pgd(k,:)=(data(:,3,i)');
    Qgd(k,:)=(data(:,4,i)');
    end
    mlinputgen((i-1)*ng+1:i*ng,1:nb)=Pgd;
    mlinputgen((i-1)*ng+1:i*ng,(nb+1):2*nb)=Qgd;
    
end
    

%getting the violated bus
activeb=(busdata(:,5)~=0|busdata(:,6)~=0);
mlinputbus(:,end+1)=activeb;


%getting the violated branch
activel=(branchdata(:,3)~=0|branchdata(:,4)~=0);
mlinputbranch(:,end+1)=activel;

mlinputgen(:,end+1:end+2)=gendata(:,2:3);

%%data formatting for Neural network

%generator
i=1:ng:ng*scenario;
inputgen=(mlinputgen(i,1:end-2));
outputgen=(mlinputgen(:,end-1:end));
outputgenp=reshape(outputgen(:,1),ng,scenario);
outputgenq=reshape(outputgen(:,2),ng,scenario);

a=find(all(inputgen==0));% find(all()) gives the indices of those columns; which columns are all zero.
inputgen_unfiltered=inputgen';
inputgen(:,a)=[];
inputgen=inputgen';


% bus
i=1:nb:nb*scenario;
inputbus=(mlinputbus(i,1:end-1));
outputbus=(mlinputbus(:,end));
outputbus=reshape(outputbus,nb,scenario);
% a=find(all(inputbus==0));% find(all()) gives the indices of those columns; which columns are all zero.
% inputbus(:,a)=[];
inputbus=inputbus';


%branch
i=1:nl:nl*scenario;
inputbranch=(mlinputbranch(i,1:end-1));
outputbranch=(mlinputbranch(:,end));
outputbranch=reshape(outputbranch,nl,scenario);
% a=find(all(inputbranch==0));% find(all()) gives the indices of those columns; which columns are all zero.
% inputbranch(:,a)=[];
inputbranch=inputbranch';


%removing the scenarios which are numerically failed to converge.
j=1;
x=[];
for i=1:scenario
    if outputbus(:,i)==1
        x(j)=i;
        j=j+1;
    end
end

if x
    inputbus(:,x)=[];
    outputbus(:,x)=[];
    inputbranch(:,x)=[];
    outputbranch(:,x)=[];
    inputgen(:,x)=[];
    inputgen_unfiltered(:,x)=[];
    outputgenp(:,x)=[];
    outputgenq(:,x)=[];
    time(x)=[];
end
    % rounding to nearest integer
    inputbus=round(inputbus);
    outputbus=round(outputbus);
    inputbranch=round(inputbranch);
    outputbranch=round(outputbranch);
    inputgen=round( inputgen);
    inputgen_unfiltered=round(inputgen_unfiltered);
    outputgenp=round(outputgenp);
    outputgenq=round(outputgenq);
    outputgen=[outputgenp;outputgenq];
%saving data


gen_index=mpc.gen(:,1);

% save('trainingdata1354_3','inputgen','inputgen_unfiltered','outputgen','outputbus','outputbranch','time','gen_index','-v7.3')
disp('done')
 
% xlswrite('dataset.xlsx',inputgen_unfiltered',1);
% xlswrite('dataset.xlsx',outputgen',2);
% xlswrite('dataset.xlsx',outputbus',3);
% xlswrite('dataset.xlsx',outputbranch',4);
% xlswrite('dataset.xlsx',gen_index,5);

