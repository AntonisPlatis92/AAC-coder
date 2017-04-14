function SMR = psycho(frameT, frameType, frameTprev1, frameTprev2)
%function to calculate the SMR of this frame
load('TableB219.mat');
sfl=zeros(69);
sfs=zeros(42);
if (strcmp(frameType,'ESH'))
    for i=1:42
        for j=1:42
            if (j>=i)
                tmpx=3*(B219b(j,5)-B219b(i,5));
            else
                tmpx=1.5*(B219b(j,5)-B219b(i,5));
            end
            tmpz=8*min((tmpx-0.5)^2-2*(tmpx-0.5),0);
            tmpy=15.811389+7.5*(tmpx+0.474)-17.5*((1.0+(tmpx+0.474)^2)^0.5);
            if (tmpy<-100)
                x=0;
            else
                x=10^((tmpz+tmpy)/10);
            end
            sfs(i,j)=x;
        end
    end
    sw=zeros(256,1);
    sw1=zeros(256,1);
    sw2=zeros(256,1);
    subFrameT1=zeros(256,1);
    index=449;
    SMR=zeros(42,8);
    for j=1:8
        subFrameT=frameT(index:index+255);
        if j==1
            subFrameT=frameTprev1(1345:1600);
            subFrameT2=frameTprev1(1217:1472);
        elseif j==2
            subFrameT1=frameT(index-128:index+127);
            subFrameT2=frameTprev1(1345:1600);
        else
            subFrameT1=frameT(index-128:index+127);
            subFrameT2=frameT(index-256:index-1);
        end
        for n=1:256
            sw(n)=subFrameT(n)*(0.5-0.5*cos((pi*(n-1+0.5))/256));
            sw1(n)=subFrameT1(n)*(0.5-0.5*cos((pi*(n-1+0.5))/256));
            sw2(n)=subFrameT2(n)*(0.5-0.5*cos((pi*(n-1+0.5))/256));
        end
        fftFrame=fft(sw);
        fftFrame=fftFrame(1:128);
        fftFrame1=fft(sw1);
        fftFrame1=fftFrame1(1:128);
        fftFrame2=fft(sw2);
        fftFrame2=fftFrame2(1:128);
        r=abs(fftFrame);
        r1=abs(fftFrame1);
        r2=abs(fftFrame2);
        f=angle(fftFrame);
        f1=angle(fftFrame1);
        f2=angle(fftFrame2);
        rPred=2*r1-r2;
        fPred=2*f1-f2;
        c=zeros(size(r,1),1);
        for w=1:size(r,1)
            c(w)=sqrt((r(w)*cos(f(w))-rPred(w)*(cos(fPred(w))))^2+(r(w)*sin(f(w))-rPred(w)*sin(fPred(w)))^2)/(r(w)+abs(rPred(w)));
        end
        e=zeros(42,1);
        c1=zeros(42,1);
        for b=1:42
            e(b)=0;
            c1(b)=0;
            for k=B219b(b,2)+1:B219b(b,3)+1
                e(b)=e(b)+r(k)^2;
                c1(b)=c1(b)+c(k)*(r(k)^2);
            end
        end
        ecb=zeros(42,1);
        ct=zeros(42,1);
        cb=zeros(42,1);
        en=zeros(42,1);
        tb=zeros(42,1);
        SNR=zeros(42,1);
        bc=zeros(42,1);
        nb=zeros(42,1);
        npart=zeros(42,1);
        for b=1:42
            ecb(b)=0;
            ct(b)=0;
            for k=1:42
                ecb(b)=ecb(b)+e(k)*sfs(k,b);
                ct(b)=ct(b)+c1(k)*sfs(k,b);
            end
            cb(b)=ct(b)/ecb(b);
            tb(b)=-0.299-0.43*log(cb(b));
            if (tb(b)<0)
                tb(b)=0;
            end
            sumSfs=0;
            for k=1:42
                sumSfs=sumSfs+sfs(k,b);
            end
            en(b)=ecb(b)/sumSfs;
            SNR(b)=tb(b)*6+(1-tb(b))*18;
            bc(b)=10.^(-SNR(b)/10);
            nb(b)=en(b)*bc(b);
            qthr=eps*256*10^(B219b(b,6)/10);
            npart(b)=max(nb(b),qthr);
            SMR(b,j)=e(b)/npart(b);
        end
    index=index+128;
    end
else
    for i=1:69
        for j=1:69
            if (i>=j)
                tmpx=3*(B219a(j,5)-B219a(i,5));
            else
                tmpx=1.5*(B219a(j,5)-B219a(i,5));
            end
            tmpz=8*min(((tmpx-0.5)^2)-2*(tmpx-0.5),0);
            tmpy=15.811389+7.5*(tmpx+0.474)-17.5*((1+(tmpx+0.474)^2)^0.5);
            if (tmpy<-100)
                x=0;
            else
                x=10^((tmpz+tmpy)/10);
            end
            sfl(i,j)=x;
        end
    end
    sw=zeros(size(frameT,1),1);
    sw1=zeros(size(frameT,1),1);
    sw2=zeros(size(frameT,1),1);
    for n=1:2048
        sw(n)=frameT(n)*(0.5-0.5*cos((pi*(n-1+0.5))/2048));
        sw1(n)=frameTprev1(n)*(0.5-0.5*cos((pi*(n-1+0.5))/size(frameT,1)));
        sw2(n)=frameTprev2(n)*(0.5-0.5*cos((pi*(n-1+0.5))/size(frameT,1)));
    end
    save('frameT', 'frameT');
    fftFrame=fft(sw);
    fftFrame=fftFrame(1:1024);
    fftFrame1=fft(sw1);
    fftFrame1=fftFrame1(1:1024);
    fftFrame2=fft(sw2);
    fftFrame2=fftFrame2(1:1024);
    r=abs(fftFrame);
    r1=abs(fftFrame1);
    r2=abs(fftFrame2);
    f=angle(fftFrame);
    f1=angle(fftFrame1);
    f2=angle(fftFrame2);
    rPred=2*r1-r2;
    fPred=2*f1-f2;
    c=zeros(size(r,1),1);
    for w=1:size(r,1)
        c(w)=sqrt((r(w)*cos(f(w))-rPred(w)*(cos(fPred(w))))^2+(r(w)*sin(f(w))-rPred(w)*sin(fPred(w)))^2)/(r(w)+abs(rPred(w)));
    end
    rSqrt=r.^2;
    e=zeros(69,1);
    c1=zeros(69,1);
    for b=1:69
        e(b)=0;
        c1(b)=0;
        for k=B219a(b,2)+1:B219a(b,3)+1
            e(b)=e(b)+r(k)^2;
            c1(b)=c1(b)+c(k)*(r(k)^2);
        end
    end
    ecb=zeros(69,1);
    ct=zeros(69,1);
    cb=zeros(69,1);
    en=zeros(69,1);
    tb=zeros(69,1);
    SNR=zeros(69,1);
    bc=zeros(69,1);
    nb=zeros(69,1);
    npart=zeros(69,1);
    SMR=zeros(69,1);
    for b=1:69
        ecb(b)=0;
        ct(b)=0;
        for k=1:69
            ecb(b)=ecb(b)+e(k)*sfl(k,b);
            ct(b)=ct(b)+c1(k)*sfl(k,b);
        end
        cb(b)=ct(b)/ecb(b);
        tb(b)=-0.299-0.43*log(cb(b));
        if (tb(b)<0)
            tb(b)=0;
        end
        sumSfl=0;
        for k=1:69
            sumSfl=sumSfl+sfl(k,b);
        end
        en(b)=ecb(b)/sumSfl;%sum(sfl(1:69,b));

        SNR(b)=tb(b)*6+(1-tb(b))*18;
        bc(b)=10.^(-SNR(b)/10);
        nb(b)=en(b)*bc(b);
        qthr=eps*1024*10^(B219a(b,6)/10);
        npart(b)=max(nb(b),qthr);
        SMR(b)=e(b)/npart(b);
    end
end
end