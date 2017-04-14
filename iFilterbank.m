function frameT = iFilterbank(frameF, frameType, winType)
% Inverse Function of filterbank. Transform MDCT Coeffs to frameT
%First we call the imdct4()
if (strcmp(frameType,'ESH'))
    preWindowFrameT=zeros(256,2,8);
    for i=1:8
        preWindowFrameT(:,1,i)=imdct4(frameF(:,1,i));
        preWindowFrameT(:,2,i)=imdct4(frameF(:,2,i));
    end
else
    preWindowFrameT=zeros(2048,2);
    preWindowFrameT(:,1)=imdct4(frameF(:,1));
    preWindowFrameT(:,2)=imdct4(frameF(:,2));
end
%Then we create the window depending on the frameType
if (strcmp(frameType,'OLS'))
    winSize=2048;
    winHalf=1024;
    a=6;
    window=zeros(2048,2);
    if (strcmp(winType,'KDB'))
      num = kaiser(1025, a);
      denum = sum(num(1:1025));
      window = zeros(2048, 2);
      for i = 1:1024
         window(i, :) = sqrt(sum(num(1:i))/denum);                 
      end
      window(1025:2048, :)= window(1024:-1:1,:);
     
    elseif (strcmp(winType,'SIN'))
      for n=1:winSize
          window(n,:)=sin((pi/winSize)*(n-1+1/2));
      end
   end

elseif (strcmp(frameType,'LSS'))
    index=1;
    winSize=2048;
    winHalf=1024;
    window=zeros(2048,2);
    if (strcmp(winType,'KDB'))
        a=6;
        num = kaiser(1025, a);
        denum = sum(num(1:1025));
        window = zeros(2048, 2);
        for i = 1:1024
           window(i, :) = sqrt(sum(num(1:i))/denum);                 
        end
        for n=1025:1472
            window(n,:)=1;
        end
        index=index+448;
        a=4;
        num = kaiser(129, a);
        denum = sum(num(1:129));
        winTest=zeros(256,2);
        for i = 1:128
           winTest(i, :) = sqrt(sum(num(1:i))/denum);                 
        end
        window(1473:1600,:)=winTest(128:-1:1,:);
        for n=1601:2048
            window(n,:)=0;
        end
    
    elseif (strcmp(winType,'SIN'))
        for n=index:winHalf
            window(n,:)=sin((pi/2048)*(n-1+1/2));
        end
        index=index+winHalf;
        for n=index:index+448-1
            window(n,:)=1;
        end
        index=index+448;
        for n=index:index+128-1
            window(n,:)=sin((pi/256)*(n+128-1473+1/2));
        end
        index=index+128;
        for n=index:winSize
            window(n,:)=0;
        end
        window1=window;
        save('Window1.mat','window1');
    end
elseif (strcmp(frameType,'LPS'))
    window=zeros(2048,2);
    if (strcmp(winType,'KDB'))
        for n=1:448
            window(n,:)=0;
        end
        a=4;
        num = kaiser(129, a);
        denum = sum(num(1:129));
        winTest=zeros(256,2);
        for i = 1:128
           winTest(i, :) = sqrt(sum(num(1:i))/denum);                 
        end
        window(449:576,:)=winTest(1:128,:);
        for n=577:1024
            window(n,:)=1;
        end
        a=6;
        num = kaiser(1025, a);
        denum = sum(num(1:1025));
        winTest = zeros(2048, 2);
        for i = 1:1024
           winTest(i, :) = sqrt(sum(num(1:i))/denum);                 
        end
        window(1025:2048, :)= winTest(1024:-1:1,:);
        
    elseif (strcmp(winType,'SIN'))
        for n=1:448
            window(n,:)=0;
        end
        for n=449:576
            window(n,:)=sin((pi/256)*(n-449+1/2));
        end
        for n=577:1024
            window(n,:)=1;
        end
        for n=1025:2048
            window(n,:)=sin((pi/2048)*(n-1+1/2));
        end
     end
elseif (strcmp(frameType,'ESH'))
    window=zeros(256,2);
    if (strcmp(winType,'KDB'))
      a=4;
      num = kaiser(129, a);
      denum = sum(num(1:129));
      window = zeros(256, 2);
      for i = 1:128
         window(i, :) = sqrt(sum(num(1:i))/denum);                 
      end
      window(129:256, :)= window(128:-1:1,:);
     
    elseif (strcmp(winType,'SIN'))
      for n=1:256
          window(n,:)=sin((pi/256)*(n-1+1/2));
      end
   end
end

%We apply the window to the pre-window frame to reconstruct frameT
if (strcmp(frameType,'ESH'))
    frameT=zeros(2048,2);
    index=449;
    windowedFrameT=zeros(256,2,8);
    for i=1:8
        windowedFrameT(:,:,i)=preWindowFrameT(:,:,i).*window;
        for n=0:255
            frameT(index+n,1)=frameT(index+n,1)+windowedFrameT(n+1,1,i);
            frameT(index+n,2)=frameT(index+n,2)+windowedFrameT(n+1,2,i);
        end
        index=index+128;
    end
    
else
    size(window);
    frameT=preWindowFrameT.*window;
end


end


