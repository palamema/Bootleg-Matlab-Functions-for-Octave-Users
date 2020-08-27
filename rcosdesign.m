function b = rcosdesign(beta,span,sps,shape)
% rcosdesign, for poor bozos like me who can't afford MATLAB

    % If no input argument is given for "shape", declare as empty string
    if nargin==3
        shape='';
    end
    
    n=-(sps*span)/2:(sps*span)/2; %Sample index
    t=n/sps; %Time index

    %First case: normal raised cosine filter
    if strcmpi(shape, 'normal')
        
        %Divide into three parts for readability:
        x=sinc(t);
        y=cos(pi*beta*t);
        z=1-(2*beta*t).^2;
        
        a=x.*y./z;

        %Fill in discontinuities:
        if beta ~= 0
            i1 = t==-1/(2*beta); %Index of first discontinuity
            i2 = t==1/(2*beta); %Index of second discontinuity
            d=pi/4*sinc(1/(2*beta)); %Value of discontinuities
            a(i1)=d;
            a(i2)=d;
        end
        
    %Second case: root raised cosine filter
    else
        
        %Divide into three parts for readability:
        x=sin(pi*t*(1-beta));
        y=4*beta*t.*cos(pi*t*(1+beta));
        z=pi*t.*(1-(4*beta*t).^2);

        a=(x+y)./z;
        
        %Fill in discontinuities:
        i0 = t==0; %Index of discontinuity at t=0
        a(i0)=1+beta*(4/pi-1); %Value of discontinuity at t=0        
        if beta ~= 0
            i1 = t==-1/(4*beta); %Index of first discontinuity
            i2 = t==1/(4*beta); %Index of last discontinuity
            %Divide into three parts for readability:
            da=beta/sqrt(2);
            db=(1+2/pi)*sin(pi/(4*beta));
            dc=(1-2/pi)*cos(pi/(4*beta));
            d=da*(db+dc); %Value of discontinuity (all parts combined)
            a(i1)=d;
            a(i2)=d;
        end
    end
    
    %Normalize filter to have a total energy of 1:
    E=sqrt(sum(a.^2)); %Energy of filter
    b=a/E;
    
end
