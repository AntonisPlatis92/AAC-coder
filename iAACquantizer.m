function frameF = iAACquantizer(S, sfc, G, frameType)
%Inverse quantization
%Given the S, sfc, G and frameType it reconstructs the frameF
load('TableB219.mat');
if (strcmp(frameType,'ESH'))
    for i=1:8
        a=zeros(42,1);
        a(1)=G(1,i);
        for b=2:42
            a(b)=sfc(b,i)+a(b-1);
        end
        for b=1:42
            for k=B219b(b,2)+1:B219b(b,3)+1
               frameF(k,i)=sign(S(k,i))*(abs(S(k,i))^(4/3))*(2^(0.25*a(b)));
            end
        end
    end
else   
    a=zeros(69,1);
    a(1)=G;
    frameF=zeros(1024,1);
    for b=2:69
        a(b)=sfc(b)+a(b-1);
    end
    for b=1:69
        for k=B219a(b,2)+1:B219a(b,3)+1
           frameF(k)=sign(S(k))*(abs(S(k))^(4/3))*(2^(0.25*a(b)));
        end
    end

end

