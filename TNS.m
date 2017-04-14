function [frameFout, TNScoeffs] = TNS(frameFin, frameType)
%Function for Temporal Noise Shaping
load('TableB219.mat');
if (strcmp(frameType,'ESH'))
    frameFout=zeros(128,8);
    TNScoeffs=zeros(4,8);
    for i=1:8
          Sw=zeros(128,1);
          P=zeros(42,1);
          for j=0:41
             P(j+1)=sum(frameFin(B219b(j+1,2)+1:B219b(j+1,3)+1,i).^2);
             for k=B219b(j+1,2):B219b(j+1,3)
                 Sw(k+1)=sqrt(P(j+1));
             end
          end
          for k = 127:-1:1
              Sw(k)=(Sw(k)+Sw(k+1))/2;
          end
          for k=2:128
              Sw(k)=(Sw(k)+Sw(k-1))/2;
          end
          Xw=frameFin(:,i)./Sw;
          cor=autocorr(Xw,4);
          r=cor(2:5);
          R=toeplitz(cor(1:4));
          a=R\r;
          for n=1:4
             a(n)=(round((a(n)*10)+0.5)-0.5)/10;
             if a(n)<-0.75
                 a(n)=-0.75;
             elseif a(n)>0.75
                a(n)=0.75;
             end
          end
          TNScoeffs(:,i)=a;
          fil=[1;-a];
          frameFout(:,i)=filter(fil,1,frameFin(:,i));     
    end        
else
    P=zeros(69,1);
    Sw=zeros(1024,1);
    for j=0:68
        P(j+1)=sum(frameFin(B219a(j+1,2)+1:B219a(j+1,3)+1).^2);
        for k=B219a(j+1,2):B219a(j+1,3)
            Sw(k+1)=sqrt(P(j+1));
        end
    end
    for k = 1023:-1:1
           Sw(k)=(Sw(k)+Sw(k+1))/2;
    end
    for k=2:1024
            Sw(k)=(Sw(k)+Sw(k-1))/2;
    end
    Xw=frameFin./Sw;
    cor=autocorr(Xw,4);
    r=cor(2:5);
    R=toeplitz(cor(1:4));
    a=R\r;
    for n=1:4
             a(n)=(round((a(n)*10)+0.5)-0.5)/10;
             if a(n)<-0.75
                 a(n)=-0.75;
             elseif a(n)>0.75
                a(n)=0.75;
             end
    end
    TNScoeffs=a;
    fil=[1;-a];
    frameFout=filter(fil,1,frameFin);
end
end

