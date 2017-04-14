function [S, sfc, G] = AACquantizer(frameF, frameType, SMR)
%Function to quantize a frameF
%Outputs: symbols S
%sfc: DPCM of global gains
%G: a(0)
load('TableB219.mat');
if (strcmp(frameType,'ESH'))
    S=zeros(128,8);
    sfc=zeros(42,8);
    G=zeros(1,8);
    for j=1:8
        P=zeros(42,1);
        T=zeros(42,1);
        for b=1:42
            P(b)=sum(frameF(B219b(b,2)+1:B219b(b,3)+1,1,j));
            T(b)=P(b)/SMR(b,j);
        end
        a=zeros(42,1);
        MQ=8191;
        a(:)=(16/3)*log2((max(frameF(:,1,j))^(3/4))/MQ);
        done=zeros(42,1);
        magicNumber=0.4054;
        X=zeros(128,1);
        while (~all(done))
            for b=1:42
                if (~done(b))
                    for k=B219b(b,2)+1:B219b(b,3)+1
                        S(k,j)=sign(frameF(k,1,j))*fix((abs(frameF(k,1,j))*(2^(-0.25*a(b))))^0.75+magicNumber);
                        X(k)=sign(S(k,j))*(abs(S(k,j))^(4/3))*2^(0.25*a(b));
                    end
                    Pe=0;
                    for k=B219b(b,2)+1:B219b(b,3)+1
                        Pe=Pe+(frameF(k,1,j)-X(k))^2;
                    end
                    if (Pe>T(b))
                        done(b)=1;
                    else
                        a(b)=a(b)+1;
                    end
                end
            end
            adif=a(41,1);
            for i=1:41
                adif(i)=a(i+1)-a(i);
            end
            if (max(adif)>60)
                done(:)=1;
                for b=1:41
                    if (a(b+1)-a(b)>60)
                        a(b+1)=a(b)+60;
                    elseif (a(b+1)-a(b)<-60)
                        a(b+1)=a(b)-60;
                    end
                    for k=B219b(b,2)+1:B219b(b,3)+1
                        S(k,j)=sign(frameF(k,1,j))*fix((abs(frameF(k,1,j))*(2^(-0.25*a(b))))^0.75+magicNumber);
                    end
                end
            end
        end
        G(1,j)=a(1);
        sfc(1,j)=round(a(1));
        for b=2:42
            sfc(b,j)=a(b)-a(b-1);
        end
    end        
else
    P=zeros(69,1);
    T=zeros(69,1);
    for b=1:size(SMR,1)
        P(b)=0;
        for k=B219a(b,2)+1:B219a(b,3)+1
            P(b)=P(b)+frameF(k)^2;
        end
        T(b)=P(b)/SMR(b);
    end
    a=zeros(69,1);
    MQ=8191;
    a(:)=(16/3)*log2((max(frameF)^(3/4))/MQ);
    done=zeros(69,1);
    magicNumber=0.4054;
    S=zeros(1024,1);
    X=zeros(1024,1);
    while (~all(done))
        for b=1:69
            if (~done(b))
                for k=B219a(b,2)+1:B219a(b,3)+1
                    S(k)=sign(frameF(k))*fix((abs(frameF(k))*(2^(-0.25*a(b))))^0.75+magicNumber);
                    X(k)=sign(S(k))*(abs(S(k))^(4/3))*2^(0.25*a(b));
                end
                Pe=0;
                for k=B219a(b,2)+1:B219a(b,3)+1
                    Pe=Pe+(frameF(k)-X(k))^2;
                end
                if (Pe>T(b))
                    done(b)=1;
                else
                    a(b)=a(b)+1;
                end
            end
        end
        adif=zeros(68,1);
        for i=1:68
            adif(i)=abs(a(i+1)-a(i));
        end
        if (max(adif)>60)
            done(:)=1;
            for b=1:68
                if (a(b+1)-a(b)>60)
                    a(b+1)=a(b)+60;
                elseif (a(b+1)-a(b)<-60)
                    a(b+1)=a(b)-60;
                end
                for k=B219a(b,2)+1:B219a(b,3)+1
                    S(k)=sign(frameF(k))*fix((abs(frameF(k))*(2^(-0.25*a(b))))^0.75+magicNumber);
                end
        end
    end
    G=a(1);
    sfc=zeros(69,1);
    sfc(1)=round(a(1));
    for b=2:69
        sfc(b)=a(b)-a(b-1);
    end
end
end

