function [result]=anv9(sequence)
tesseq=sequence;
n=length(tesseq);
na=0;
nc=0;
ng=0;
nt=0;
ua=0;
uc=0;
ug=0;
ut=0;
daup=0;
dcup=0;
dgup=0;
dtup=0;
covacup=0;
covagup=0;
covatup=0;
covcgup=0;
covctup=0;
covgtup=0;
numerical=zeros(1,n);
for k=1:n
    if tesseq(k)=='A'|tesseq(k)=='a'
        numerical(k)=1;
        na=na+1;
    else if tesseq(k)=='C'|tesseq(k)=='c'
            numerical(k)=2;
            nc=nc+1;
        else if tesseq(k)=='G'|tesseq(k)=='g'
                numerical(k)=3;
                ng=ng+1;
            else if tesseq(k)=='T'|tesseq(k)=='t'
                    numerical(k)=4;
                    nt=nt+1;
                end
            end
        end
    end
end

supa=0;
supc=0;
supg=0;
supt=0;
posa=zeros(1,na);
posc=zeros(1,nc);
posg=zeros(1,ng);
post=zeros(1,nt);
a=1;
c=1;
g=1;
t=1;
indicatora=zeros(1,n);
indicatorc=zeros(1,n);
indicatorg=zeros(1,n);
indicatort=zeros(1,n);
accuindia=zeros(1,n);
accuindic=zeros(1,n);
accuindig=zeros(1,n);
accuindit=zeros(1,n);
if numerical(1)==1
     posa(a)=1;
     supa=supa+1;
     a=a+1;
     indicatora(1)=1;
 else if numerical(1)==2
            posc(c)=1;
            supc=supc+1;
            c=c+1;
           indicatorc(1)=1;
        else if numerical(1)==3
                posg(g)=1;
                supg=supg+1;
                g=g+1;
                indicatorg(1)=1;
            else if numerical(1)==4
                    post(t)=k;
                    supt=supt+k;
                    t=t+1;
                    indicatort(1)=1;
                end
            end
      end
end
  
for k=2:n
    if numerical(k)==1
        posa(a)=k;
        supa=supa+k;
        a=a+1;
        indicatora(k)=1;
    else if numerical(k)==2
            posc(c)=k;
            supc=supc+k;
            c=c+1;
           indicatorc(k)=1;
        else if numerical(k)==3
                posg(g)=k;
                supg=supg+k;
                g=g+1;
                indicatorg(k)=1;
            else if numerical(k)==4
                    post(t)=k;
                    supt=supt+k;
                    t=t+1;
                    indicatort(k)=1;
                end
            end
        end
    end
end



for k=1:n
    for j=1:k
        accuindia(k)=accuindia(k)+indicatora(j);
        accuindic(k)=accuindic(k)+indicatorc(j);
        accuindig(k)=accuindig(k)+indicatorg(j);
        accuindit(k)=accuindit(k)+indicatort(j);
    end
end


%u2
ua2=sum(accuindia)/na;
uc2=sum(accuindic)/nc;
ug2=sum(accuindig)/ng;
ut2=sum(accuindit)/nt;
%u3
ua3=mean(accuindia);
uc3=mean(accuindic);
ug3=mean(accuindig);
ut3=mean(accuindit);



for k=1:n
    daup=daup+(accuindia(k)-ua3)*(accuindia(k)-ua3);
    dcup=dcup+(accuindic(k)-uc3)*(accuindic(k)-uc3);
    dgup=dgup+(accuindig(k)-ug3)*(accuindig(k)-ug3);
    dtup=dtup+(accuindit(k)-ut3)*(accuindit(k)-ut3);
end

for k=1:n
    covacup=covacup+(accuindia(k)-ua3)*(accuindic(k)-uc3);
    covagup=covagup+(accuindia(k)-ua3)*(accuindig(k)-ug3);
    covatup=covatup+(accuindia(k)-ua3)*(accuindit(k)-ut3);
    covcgup=covcgup+(accuindic(k)-uc3)*(accuindig(k)-ug3);
    covctup=covctup+(accuindic(k)-uc3)*(accuindit(k)-ut3);
    covgtup=covgtup+(accuindig(k)-ug3)*(accuindit(k)-ut3);
end

%D0618%
da0618=daup/(na*na);
dc0618=dcup/(nc*nc);
dg0618=dgup/(ng*ng);
dt0618=dtup/(nt*nt);

cov0618ac=covacup/(na*nc);
cov0618ag=covagup/(na*ng);
cov0618at=covatup/(na*nt);
cov0618cg=covcgup/(nc*ng);
cov0618ct=covctup/(nc*nt);
cov0618gt=covgtup/(ng*nt);

result=[na,nc,ng,nt,ua2,uc2,ug2,ut2,da0618,dc0618,dg0618,dt0618,cov0618ac,cov0618ag,cov0618at,cov0618cg,cov0618ct,cov0618gt];